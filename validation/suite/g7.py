"""G7 regression: the final wall force of every tutorial against the stored baseline,
|F - F_base| <= 1e-6 |F_base| per body. A larger difference is either a defect or an intended
change; in the latter case the baseline is updated in the same commit with the reason."""
import re

import common

ID, KIND, TIER = "G7", "gate", 0
TITLE = "regression: forces of every tutorial"
STRESSES, RELIES = ["K1", "K3", "K4", "K5", "K6", "K8", "K9"], []
REGRESSION = None   # own check below (relative to the force magnitude of each body)

CASES = ["cylinderRe1", "cylinderRe1_stl", "cylinderRe1_alphaField", "cylinderRe1_pimple",
         "cylinderRe1_pimple_freeSlip", "twoCylinders", "oscillatingCylinder",
         "galileanTranslating", "order/taylorCouette"]
RTOL = 1e-6


def one(ctx, name):
    d = common.copy_case(common.TUTORIALS / name, ctx.root / name.replace("/", "_"))
    app = re.search(r"application\s+(\w+)", (d / "system" / "controlDict").read_text()).group(1)
    for cmd in ("blockMesh", "ibmSetGeometry", app):
        common.sh(cmd, d, "log." + cmd)
    t = (d / ("log." + app)).read_text()
    bodies = sorted(set(re.findall(r"IBM force on (\S+):", t)))
    return name, {b: common.ibm_force(t, b)["total"].tolist() for b in bodies}


def run(ctx, res):
    rows = common.pmap(lambda n: one(ctx, n), CASES, ctx.jobs)
    res.data["forces"] = {n: f for n, f in rows}
    res.data["metrics"] = {"%s/%s_F%s" % (n, b, "xyz"[i]): F[i] for n, f in rows for b, F in f.items()
                           for i in range(3)}
    base = common.load_baseline(ID)
    if base is None:
        res.report("baseline", len(rows), note="no stored baseline: this run can be accepted as one")
        return
    bf = base["data"]["forces"]
    for n, f in rows:
        for b, F in f.items():
            if n not in bf or b not in bf[n]:
                res.report("%s/%s: not in the baseline" % (n, b), None)
                continue
            Fb = bf[n][b]
            mag = max(sum(x*x for x in Fb)**0.5, 1e-300)
            d = max(abs(x - y) for x, y in zip(F, Fb))/mag
            res.le("%s/%s: max |F - F_base|/|F_base|" % (n, b), d, RTOL,
                   "baseline %s" % base["commit"])
