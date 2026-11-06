"""V7 Spalart-Allmaras plane channel, Re_b = 2e4 (Re_tau ~ 545, tutorials/turbulence/channelSA):
immersed walls against the body-fitted channel with the same model (validation/channelSA_compare.py).
Tier 1: wall-resolved (y+ ~ 0.5), the walls cutting their cells at f = 0.05, 0.3, 0.5, linked
(default) and unlinked at f = 0.3. Tier 2 (not here yet): wall functions, inclined walls."""
import re
import subprocess

import numpy as np

import common
from channelSA_compare import load

ID, KIND, TIER = "V7", "case", 1
TITLE = "SA channel, immersed vs body-fitted walls (wall-resolved)"
STRESSES, RELIES = ["K12", "K13", "K10"], ["V1", "V2", "V5"]
REGRESSION = (1e-4,)

MAKE = common.TUTORIALS / "turbulence" / "channelSA" / "makeCase.py"


def one(ctx, kind, f=None, ma=0.5):
    name = "bf_wr" if kind == "bf" else "ibm_wr_f%g_ma%g" % (f, ma)
    d = ctx.root / name
    if not (common.REUSE and d.exists()):
        args = [str(d), "--kind", kind, "--mesh", "wr"]
        if kind == "ibm":
            args += ["--offset", str(f), "--mergeAlpha", str(ma)]
        subprocess.run([common.sys.executable, str(MAKE)] + args, check=True, capture_output=True)
        common.set_keys(d / "system" / "controlDict", endTime=60000)
    app = re.search(r"application\s+(\w+)", (d / "system" / "controlDict").read_text()).group(1)
    if (d / "system" / "blockMeshDict").exists():
        common.sh("blockMesh", d, "log.blockMesh")
    (d / "0").mkdir(exist_ok=True)
    common.sh("postProcess -func writeCellCentres -time 0", d, "log.postProcess")
    if kind == "ibm":
        common.sh("ibmSetGeometry", d, "log.ibmSetGeometry")
    common.sh(app, d, "log." + app)
    log = (d / ("log." + app)).read_text()
    r = load(str(d))
    out = {"name": name, "kind": kind, "f": f, "ma": ma, "utau": float(r["utau"]),
           "converged": "SIMPLE solution converged" in log}
    for side, m in (("bottom", r["bottom"]), ("top", ~r["bottom"])):
        o = np.argsort(r["yp"][m])
        for q in (5, 30, 100, 300):
            out["%s_u+(%d)" % (side, q)] = float(np.interp(q, r["yp"][m][o], r["up"][m][o]))
    return out


def run(ctx, res):
    jobs = [("bf",), ("ibm", 0.05), ("ibm", 0.3), ("ibm", 0.5), ("ibm", 0.3, 0)]
    rows = common.pmap(lambda j: one(ctx, *j), jobs, ctx.jobs)
    res.data["rows"] = rows
    res.data["metrics"] = {r["name"] + "_utau": r["utau"] for r in rows}
    res.le("runs not converged (SIMPLE residualControl, 60000 iterations)", sum(not r["converged"] for r in rows), 0)
    bf = [r for r in rows if r["kind"] == "bf"][0]
    for r in rows:
        if r["kind"] == "bf":
            continue
        du = max(abs(r["%s_u+(%d)" % (s, q)] - bf["%s_u+(%d)" % (s, q)]) for s in ("bottom", "top")
                 for q in (5, 30, 100, 300))
        res.rel("u_tau, f %g, mergeAlpha %g, vs body-fitted" % (r["f"], r["ma"]), r["utau"], bf["utau"], 0.005,
                "max |du+| at y+ = 5, 30, 100, 300: %.3f" % du)
