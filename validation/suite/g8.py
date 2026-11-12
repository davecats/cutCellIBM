"""G8 decomposition independence: converged steady runs serial and on 4 and 6 ranks must give
the same forces. cylinderRe1 (scotch: cut cells across processor faces; and 5 ranks cut along
x, whose processor face at the left of the cylinder splits linked groups: 6 links cross it)
and the SA channel with immersed walls (SA wall law; decomposed across the channel, so the
IBM wall distance and the wall samples cross processors)."""
import re
import subprocess

import numpy as np

import common
from common import write_dict

ID, KIND, TIER = "G8", "gate", 0
TITLE = "decomposition: serial vs 4 and 6 ranks"
STRESSES, RELIES = ["K14"], ["G7"]
REGRESSION = None
RTOL = 1e-8


def decompose(d, n, method, split="1 %d 1"):
    coeffs = "coeffs { n (%s); }" % (split % n) if method == "simple" else ""
    write_dict(d / "system" / "decomposeParDict", "dictionary", "decomposeParDict",
               "numberOfSubdomains %d;\nmethod %s;\n%s\n" % (n, method, coeffs))


def cylinder(ctx, n, method="scotch"):
    d = common.copy_case(common.TUTORIALS / "cylinderRe1", ctx.root / ("cylinderRe1_np%d" % n))
    common.sh("blockMesh", d, "log.blockMesh")
    # the solver builds the geometry itself; ibmSetGeometry only writes fields (its 0/Cf with
    # sliced patches cannot be decomposed)
    if n > 1:
        decompose(d, n, method, "%d 1 1")
        common.sh("decomposePar", d, "log.decomposePar")
    common.sh("ibmSimpleFoam", d, "log.ibmSimpleFoam", np_=n)
    t = (d / "log.ibmSimpleFoam").read_text()
    m = re.search(r"cell linking: .*\((\d+) over coupled faces\)", t)
    return np.append(common.ibm_force(t, "cylinder")["total"], int(m.group(1)) if m else 0)


def cylinder_split(ctx, n):
    """slabs along x: the processor face next to the leftmost point of the cylinder separates
    small cut cells from their masters"""
    return cylinder(ctx, n, "simple")


def channel(ctx, n):
    d = ctx.root / ("channelSA_np%d" % n)
    subprocess.run([common.sys.executable, str(common.TUTORIALS / "turbulence" / "channelSA" / "makeCase.py"),
                    str(d), "--kind", "ibm", "--mesh", "wr"], check=True, capture_output=True)
    ctl = (d / "system" / "controlDict").read_text()
    app = re.search(r"application\s+(\w+)", ctl).group(1)
    if (d / "system" / "blockMeshDict").exists():
        common.sh("blockMesh", d, "log.blockMesh")
    if n > 1:
        decompose(d, n, "simple")
        common.sh("decomposePar", d, "log.decomposePar")
    common.sh(app, d, "log." + app, np_=n)
    t = (d / ("log." + app)).read_text()
    # the driving gradient (u_tau^2); residualControl cannot trigger in this 1-D problem (Uy and p
    # residuals normalised by ~0), so the iteration must be stationary (with nuTilda relaxed by 0.9
    # it ended in a two-iteration cycle whose end state depended on the decomposition)
    g, change = common.driving_gradient(t)
    return np.array([g[-1], change()])


def run(ctx, res):
    jobs = [(f, n) for f in (cylinder, channel) for n in (1, 4, 6)] + [(cylinder_split, 5)]
    out = common.pmap(lambda j: (j[0].__name__, j[1], j[0](ctx, j[1])), jobs, max(1, ctx.jobs // 4))
    F = {(k, n): v for k, n, v in out}
    res.data["forces"] = {"%s_np%d" % k: v for k, v in F.items()}
    for n in (4, 6):
        c, s = F[("cylinder", n)], F[("cylinder", 1)]
        res.le("cylinder: |drag_np%d/drag_serial - 1|" % n, abs(c[0]/s[0] - 1), RTOL)
        res.report("cylinder: lift np%d vs serial (zero by symmetry: solver noise)" % n, float(c[1]), float(s[1]))
        g, gs = F[("channel", n)][0], F[("channel", 1)][0]
        res.le("channel: |G_np%d/G_serial - 1| (driving gradient, u_tau^2)" % n, abs(g/gs - 1), RTOL)
    c, s = F[("cylinder_split", 5)], F[("cylinder", 1)]
    res.ge("cylinder, 5 slabs along x: links over processor faces", int(c[3]), 1)
    res.le("cylinder, 5 slabs along x: |drag/drag_serial - 1|", abs(c[0]/s[0] - 1), RTOL)
    for n in (1, 4, 6):
        res.le("channel np%d: driving gradient stationary (largest change over 100 iterations)" % n,
               F[("channel", n)][1], 1e-10)
