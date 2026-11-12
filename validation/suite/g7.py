"""G7 regression: the final wall force and moment of every tutorial against the stored baseline,
per body |F - F_base| <= 1e-6 max(|F_base|, |M_base|/L) and |M - M_base| <= 1e-6 max(|M_base|, L |F_base|),
L the size of the body (largest distance of its wall cells from their area centroid). The load
scale covers bodies whose net force is a cancellation (Taylor-Couette: |F| ~ 1e-8 against a torque
of 0.28) and those whose moment is (a cylinder in a stream). A larger difference is either a
defect or an intended change; in the latter case the baseline is updated in the same commit with
the reason."""
import re

import numpy as np

import common
from foamio import read_field

ID, KIND, TIER = "G7", "gate", 0
TITLE = "regression: forces and moments of every tutorial"
STRESSES, RELIES = ["K1", "K3", "K4", "K5", "K6", "K8", "K9"], []
REGRESSION = None   # own check below (relative to the load scale of each body)

CASES = ["cylinderRe1", "cylinderRe1_stl", "cylinderRe1_alphaField", "cylinderRe1_pimple",
         "cylinderRe1_pimple_freeSlip", "twoCylinders", "oscillatingCylinder",
         "galileanTranslating", "order/taylorCouette"]
RTOL = 1e-6


def sizes(d, bodies):
    """{body: L} from the geometry written at time 0 (ibmBodyIndex, ibmSw, cell centres)"""
    common.sh("postProcess -func writeCellCentres -time 0", d, "log.writeCellCentres")
    C = read_field(d / "0" / "C")[0]
    idx = np.rint(read_field(d / "0" / "ibmBodyIndex")[0]).astype(int)
    Sw = np.linalg.norm(read_field(d / "0" / "ibmSw")[0], axis=1)
    out = {}
    for i, b in enumerate(bodies):
        w = (idx == i) & (Sw > 0)
        c = (Sw[w, None]*C[w]).sum(axis=0)/Sw[w].sum()
        out[b] = float(np.linalg.norm(C[w] - c, axis=1).max())
    return out


def one(ctx, name):
    d = common.copy_case(common.TUTORIALS / name, ctx.root / name.replace("/", "_"))
    app = re.search(r"application\s+(\w+)", (d / "system" / "controlDict").read_text()).group(1)
    for cmd in ("blockMesh", "ibmSetGeometry", app):
        common.sh(cmd, d, "log." + cmd)
    # body order of ibmProperties = ibmBodyIndex
    bodies = list(dict.fromkeys(re.findall(r"^\s*body (\S+):", (d / "log.ibmSetGeometry").read_text(), re.M)))
    t = (d / ("log." + app)).read_text()
    return name, {b: common.ibm_force(t, b)["total"].tolist() for b in bodies}, \
        {b: common.ibm_moment(t, b)["total"].tolist() for b in bodies}, sizes(d, bodies)


def run(ctx, res):
    rows = common.pmap(lambda n: one(ctx, n), CASES, ctx.jobs)
    res.data["forces"] = {n: f for n, f, _, _ in rows}
    res.data["moments"] = {n: m for n, _, m, _ in rows}
    res.data["sizes"] = {n: L for n, _, _, L in rows}
    res.data["metrics"] = {"%s/%s_%s%s" % (n, b, q, "xyz"[i]): V[i] for n, f, m, _ in rows
                           for q, fm in (("F", f), ("M", m)) for b, V in fm.items() for i in range(3)}
    base = common.load_baseline(ID)
    if base is None:
        res.report("baseline", len(rows), note="no stored baseline: this run can be accepted as one")
        return
    bf, bm = base["data"]["forces"], base["data"].get("moments", {})
    for n, f, m, L in rows:
        for b in f:
            if n not in bf or b not in bf[n]:
                res.report("%s/%s: not in the baseline" % (n, b), None)
                continue
            Fb = np.array(bf[n][b])
            # a baseline without moments (before G7 compared them): this run's moment sets the scale
            Mb = np.array(bm[n][b]) if b in bm.get(n, {}) else None
            Mscale = np.linalg.norm(Mb if Mb is not None else m[b])
            scale = max(np.linalg.norm(Fb), Mscale/L[b], 1e-300)
            res.le("%s/%s: max |F - F_base| / max(|F_base|, |M_base|/L)" % (n, b),
                   float(np.abs(np.array(f[b]) - Fb).max()/scale), RTOL,
                   "baseline %s, L %.4g" % (base["commit"], L[b]))
            if Mb is None:
                res.report("%s/%s: moment not in the baseline" % (n, b), None)
                continue
            scale = max(Mscale, L[b]*np.linalg.norm(Fb), 1e-300)
            res.le("%s/%s: max |M - M_base| / max(|M_base|, L |F_base|)" % (n, b),
                   float(np.abs(np.array(m[b]) - Mb).max()/scale), RTOL,
                   "baseline %s, L %.4g" % (base["commit"], L[b]))
