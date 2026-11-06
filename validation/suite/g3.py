"""G3 hydrostatic container (tutorials/order/hydrostatic): fluid at rest under gravity in an
immersed container around an immersed obstacle; exact U = 0, p = -y + const. The spurious
velocity must converge at second order and not grow against the baseline. Linked (default)
and unlinked."""
import os

import numpy as np

import common
from foamio import latest_time, read_field

ID, KIND, TIER = "G3", "gate", 0
TITLE = "hydrostatic container: spurious velocity"
STRESSES, RELIES = ["K5", "K6", "K2"], ["G1", "G2"]
REGRESSION = (0.05,)

SRC = common.TUTORIALS / "order" / "hydrostatic"


def one(ctx, N, ma):
    d = ctx.root / ("N%d_ma%g" % (N, ma))
    common.copy_case(SRC, d)
    common.set_keys(d / "system" / "blockMeshDict", N=N)
    common.set_keys(d / "constant" / "ibmProperties", mergeAlpha=ma)
    for cmd in ("blockMesh", "ibmSetGeometry", "postProcess -func writeCellCentres -time 0", "ibmPimpleFoam"):
        common.sh(cmd, d, "log." + cmd.split()[0])
    t = latest_time(str(d))
    C, _ = read_field(os.path.join(d, "0", "C"))
    U, _ = read_field(os.path.join(d, t, "U"))
    p, _ = read_field(os.path.join(d, t, "p"))
    geo = t if os.path.exists(os.path.join(d, t, "ibmAlpha")) else "0"
    a, _ = read_field(os.path.join(d, geo, "ibmAlpha"))
    offf = os.path.join(d, t, "ibmCentroidOffset")
    off = read_field(offf)[0] if os.path.exists(offf) else np.zeros_like(C)
    y = C[:, 1] + off[:, 1]
    fl = a > 0
    ep = p + y
    ep -= np.mean(ep[fl & (a == 1)])
    mU = np.linalg.norm(U, axis=1)
    return {"N": N, "ma": ma, "t": t, "maxU": float(mU[fl].max()),
            "rmsU": float(np.sqrt(np.mean(mU[fl]**2))), "maxEp": float(np.abs(ep[fl]).max())}


def run(ctx, res):
    Ns = [40, 80] if ctx.quick else [40, 80, 160]
    rows = common.pmap(lambda j: one(ctx, *j), [(N, ma) for N in Ns for ma in (0.5, 0)], ctx.jobs)
    res.data["rows"] = rows
    res.data["metrics"] = {"N%d_ma%g_%s" % (r["N"], r["ma"], k): r[k] for r in rows for k in ("maxU", "rmsU")}
    for ma in (0.5, 0):
        s = sorted([r for r in rows if r["ma"] == ma], key=lambda r: r["N"])
        note = "max|U| " + " ".join("%.3e" % r["maxU"] for r in s)
        res.ge("max|U| order, mergeAlpha %g" % ma, common.order([1/r["N"] for r in s], [r["maxU"] for r in s]),
               1.8, note)
        res.report("pressure error order (max), mergeAlpha %g" % ma,
                   common.order([1/r["N"] for r in s], [r["maxEp"] for r in s]),
                   note="max|p + y| " + " ".join("%.3e" % r["maxEp"] for r in s))
