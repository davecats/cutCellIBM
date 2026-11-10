"""V1 Taylor-Couette flow between two immersed cylinders (tutorials/order/taylorCouette): exact
u_theta = A r + B/r and pressure; L2 errors over all fluid cells and over the cut cells,
N = 40...320, linked (default) and unlinked; secondOrder false reported. Torque on each cylinder
(postProcessing/ibmForces/<t>/moments.dat, C1) against the exact 4 pi nu B per unit depth (the inner
one rotates: it includes the transpose stress nu (n ^ omega) of the rotating wall)."""
import glob
import os

import numpy as np

import common
from foamio import latest_time, read_field

ID, KIND, TIER = "V1", "case", 1
TITLE = "Taylor-Couette, exact solution, N = 40...320"
STRESSES, RELIES = ["K3", "K5", "K6", "K2"], ["G1", "G2"]
REGRESSION = (0.01,)

SRC = common.TUTORIALS / "order" / "taylorCouette"
RI, RO, OM, CX, CY = 0.5, 1.5, 1.0, 1.6, 1.6
NU, L = 1.0, 3.2   # constant/transportProperties, system/blockMeshDict (square of side L, depth L/N)


def torques(d, N):
    """relative errors of the final torque on the inner and the outer cylinder (bodies 0 and 1)"""
    f = sorted(glob.glob(os.path.join(d, "postProcessing", "ibmForces", "*", "moments.dat")),
               key=lambda q: float(os.path.basename(os.path.dirname(q))))[-1]
    m = np.loadtxt(f, ndmin=2)[-1]
    dz = L/N
    exact = 4*np.pi*NU*OM*RI**2*RO**2/(RO**2 - RI**2)
    return {"Min": float((-(m[3] + m[6])/dz - exact)/exact), "Mout": float(((m[9] + m[12])/dz - exact)/exact)}


def errors(d):
    A = -OM*RI**2/(RO**2 - RI**2)
    B = OM*RI**2*RO**2/(RO**2 - RI**2)
    t = latest_time(str(d))
    C, _ = read_field(os.path.join(d, "0", "C"))
    U, _ = read_field(os.path.join(d, t, "U"))
    p, _ = read_field(os.path.join(d, t, "p"))
    a, _ = read_field(os.path.join(d, t, "ibmAlpha"))
    xv = common.value_points(d, t, C)
    x, y = xv[:, 0] - CX, xv[:, 1] - CY
    r = np.hypot(x, y)
    fl = a > 0
    ut = A*r + B/r
    eU = np.hypot(U[:, 0] + ut*y/r, U[:, 1] - ut*x/r)
    pex = A*A*r*r/2 + 2*A*B*np.log(np.maximum(r, 1e-300)) - B*B/(2*r*r)
    ep = p - pex
    ep -= np.sum(ep[fl]*a[fl])/np.sum(a[fl])

    def l2(e, m):
        return float(np.sqrt(np.sum(e[m]**2*a[m])/np.sum(a[m])))
    cut = fl & (a < 1)
    return {"uL2": l2(eU, fl), "uL2cut": l2(eU, cut), "uLinf": float(eU[fl].max()),
            "pL2": l2(ep, fl), "pL2cut": l2(ep, cut)}


def one(ctx, N, ma, so):
    d = ctx.root / ("N%d_ma%g_so%d" % (N, ma, so))
    common.copy_case(SRC, d)
    common.set_keys(d / "system" / "blockMeshDict", N=N)
    common.set_keys(d / "constant" / "ibmProperties", mergeAlpha=ma, secondOrder=bool(so))
    # pseudo-time to the steady state: the viscous time R^2/nu = 1 over 6 units at dt 0.02
    for cmd in ("blockMesh", "postProcess -func writeCellCentres -time 0", "ibmSetGeometry", "ibmPimpleFoam"):
        common.sh(cmd, d, "log." + cmd.split()[0])
    return {"N": N, "ma": ma, "so": so, **errors(d), **torques(d, N)}


def run(ctx, res):
    Ns = [40, 80] if ctx.quick else [40, 80, 160, 320]
    jobs = [(N, ma, 1) for N in Ns for ma in (0.5, 0)] + [(N, 0, 0) for N in Ns]
    rows = common.pmap(lambda j: one(ctx, *j), jobs, ctx.jobs)
    res.data["rows"] = rows
    res.data["metrics"] = {"N%d_ma%g_so%d_%s" % (r["N"], r["ma"], r["so"], k): r[k]
                           for r in rows for k in ("uL2", "pL2", "uL2cut", "Min", "Mout")}
    for ma, so in ((0.5, 1), (0, 1), (0, 0)):
        s = sorted([r for r in rows if (r["ma"], r["so"]) == (ma, so)], key=lambda r: r["N"])
        h = [1/r["N"] for r in s]
        lab = "mergeAlpha %g%s" % (ma, "" if so else ", secondOrder false")
        # p: the cut-cell pressure converges below second order (a continuity error in the cut
        # cells, doc/secondOrderStudy.md): reported, open
        for k, name, tol in (("uL2", "u L2 order (all fluid cells)", 1.8),
                             ("pL2", "p L2 order (all fluid cells)", None),
                             ("uL2cut", "u L2 order (cut cells)", None)):
            p = common.order(h, [r[k] for r in s])
            note = " ".join("%.3e" % r[k] for r in s) + ("; open: doc/secondOrderStudy.md" if k == "pL2" else "")
            if so and tol is not None:
                res.ge("%s, %s" % (name, lab), p, tol, note)
            else:
                res.report("%s, %s" % (name, lab), p, note=note)
        for k, name in (("Min", "inner (rotating)"), ("Mout", "outer (at rest)")):
            e = [abs(r[k]) for r in s]
            p = common.order(h, e)
            note = "rel. error " + " ".join("%+.2e" % r[k] for r in s)
            if so and ma > 0:
                res.ge("torque order, %s cylinder, %s" % (name, lab), p, 1.8, note)
            else:
                res.report("torque order, %s cylinder, %s" % (name, lab), p, note=note)
        if so and ma > 0 and s[-1]["N"] >= 320:
            res.le("torque rel. error, both cylinders, N = %d, %s" % (s[-1]["N"], lab),
                   max(abs(s[-1]["Min"]), abs(s[-1]["Mout"])), 1e-3)
