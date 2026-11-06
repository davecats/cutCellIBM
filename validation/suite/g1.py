"""G1 geometry: fluid volume and wall area of a cylinder (2-D, outside and inside) and a sphere
(3-D, outside and inside) against the exact values, closure of the wall segments, uniform field."""
import math

import numpy as np

import common
import geometry

ID, KIND, TIER = "G1", "gate", 0
TITLE = "geometry: volume, wall area, closure (cylinder 2-D, sphere 3-D, inverted)"
STRESSES, RELIES = ["K1"], []
REGRESSION = (1e-9,)

R = 1.0
C2 = (1.5123, 1.4871, 0.0)    # off the mesh lines
C3 = (1.5123, 1.4871, 1.5033)
SCHEMES = common.TUTORIALS / "order" / "laplaceDisk" / "system"


def make(ctx, dim, N, invert, ns):
    d = ctx.root / ("%dd_N%d_inv%d_ns%d" % (dim, N, int(invert), ns))
    L = 3.0
    if dim == 2:
        geometry.box(d, (0, 0, 0), (L, L, L/N), (N, N, 1))
        body = {"type": "cylinder", "centre": C2, "axis": (0, 0, 1), "radius": R}
    else:
        geometry.box(d, (0, 0, 0), (L, L, L), (N, N, N))
        body = {"type": "sphere", "centre": C3, "radius": R}
    body["invert"] = invert
    geometry.ibm_properties(d, {"body": body}, nSubSamples=ns, writeFields=False, writeForces=False)
    geometry.control_dict(d, "ibmSetGeometry")
    for f in ("fvSchemes", "fvSolution"):
        (d / "system" / f).write_text((SCHEMES / f).read_text())
    common.sh("blockMesh", d, "log.blockMesh")
    common.sh("ibmSetGeometry", d, "log.ibmSetGeometry")
    t = (d / "log.ibmSetGeometry").read_text()
    V = common.last_match(t, r"fluid volume ([-+0-9.eE]+)")
    A = common.last_match(t, r"total wall area ([-+0-9.eE]+)")
    net = np.linalg.norm(common.vec(common.last_match(t, r"net wall area vector (\([^)]*\))", cast=str)))
    g = common.last_match(t, r"uniform field, max magnitude: ([-+0-9.eE]+)")
    nabs = common.last_match(t, r"absorbed (\d+) sliver", cast=int)
    dz = L/N if dim == 2 else 1.0
    if dim == 2:
        Vb, Ab = math.pi*R*R*dz, 2*math.pi*R*dz
        Vd = L*L*dz
    else:
        Vb, Ab = 4/3*math.pi*R**3, 4*math.pi*R*R
        Vd = L**3
    Vex = Vb if invert else Vd - Vb
    return {"dim": dim, "N": N, "invert": invert, "ns": ns, "h": L/N,
            "eV": abs(V - Vex)/Vb, "eA": abs(A - Ab)/Ab, "net": net/A, "grad": g, "absorbed": nabs}


def run(ctx, res):
    N2 = [20, 40] if ctx.quick else [20, 40, 80, 160]
    N3 = [16, 32] if ctx.quick else [16, 32, 64]
    jobs = [(2, N, inv, 32) for inv in (False, True) for N in N2]
    jobs += [(2, N, False, 8) for N in N2]
    jobs += [(3, N, inv, 16) for inv in (False, True) for N in N3]
    rows = common.pmap(lambda j: make(ctx, *j), jobs, ctx.jobs)
    res.data["rows"] = rows
    res.data["metrics"] = {"%dd_N%d_inv%d_ns%d_%s" % (r["dim"], r["N"], r["invert"], r["ns"], k): r[k]
                           for r in rows for k in ("eV", "eA")}
    res.le("closure: max |net wall area|/area (closed bodies)", max(r["net"] for r in rows), 1e-10)
    res.le("closure: max |cut-cell gradient of a uniform field|", max(r["grad"] for r in rows), 1e-10)
    for dim, inv, ns in [(2, False, 32), (2, True, 32), (3, False, 16), (3, True, 16), (2, False, 8)]:
        s = [r for r in rows if (r["dim"], r["invert"], r["ns"]) == (dim, inv, ns)]
        lab = "%s %s, nSubSamples %d" % ("cylinder" if dim == 2 else "sphere",
                                        "inside" if inv else "outside", ns)
        pa = common.order([r["h"] for r in s], [r["eA"] for r in s])
        note = "rel. volume error %s; absorbed slivers %s" % (
            " ".join("%.2e" % r["eV"] for r in s), " ".join(str(r["absorbed"]) for r in s))
        if ns == 8:
            res.report("volume error (max rel.), " + lab, max(r["eV"] for r in s), note=note)
            res.report("wall area order, " + lab, pa,
                       note="rel. area error %s" % " ".join("%.2e" % r["eA"] for r in s))
        else:
            # sub-sampling quadrature, not a discretisation error: bounded, no order
            res.le("volume error (max rel. to the body), " + lab, max(r["eV"] for r in s), 1e-4, note)
            res.ge("wall area order, " + lab, pa, 1.8,
                   "rel. area error %s" % " ".join("%.2e" % r["eA"] for r in s))
