#!/usr/bin/env python3
"""Taylor-Couette errors by cell class and wall torque, for first-order (and second-order) runs
with or without linking (tutorials/order/taylorCouette geometry: Ri 0.5, Ro 1.5, centre (1.6 1.6),
inner cylinder at Omega = 1). Needs the geometry fields of ibmSetGeometry in 0/ (ibmAlpha,
ibmLinkRole, ibmNoSlipCoeff, ibmSw, C) and U, p, ibmUb at the latest time.

Classes, per wall (inner r < 1, outer r > 1):
  slave    linked cell without the master role (role 1)
  master   group master (role 2)
  cut      unlinked cut cell (0 < alpha < 1, no role)
  layer1   full cell (alpha = 1, no role) with a face neighbour that is cut or solid
  bulk     every other full cell
Errors at the cell centres (the first-order value points), velocity |U - Uex|, pressure with
the alpha-weighted fluid mean of the error removed; rms over the class (alpha-weighted).
Torque: the torque the solution transmits, 4 pi nu b from a least-squares fit u_theta = a r + b/r over
the bulk cells, against the exact 4 pi nu B (steady: the same through every circle).
Usage: taylorCouette_classes.py case...
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from foamio import read_field, latest_time

Ri, Ro, Om, cx, cy, L = 0.5, 1.5, 1.0, 1.6, 1.6, 3.2
A = -Om*Ri**2/(Ro**2 - Ri**2); B = Om*Ri**2*Ro**2/(Ro**2 - Ri**2)
CLASSES = ("slave", "master", "cut", "layer1", "bulk")


def exact(x, y):
    r = np.hypot(x, y); ut = A*r + B/r
    return np.stack([-ut*y/r, ut*x/r], 1), A*A*r*r/2 + 2*A*B*np.log(r) - B*B/(2*r*r)


def analyse(case):
    t = latest_time(case)
    rd0 = lambda f: read_field(os.path.join(case, "0", f))[0]
    C = rd0("C"); a = rd0("ibmAlpha"); c = rd0("ibmNoSlipCoeff")
    role = rd0("ibmLinkRole") if os.path.exists(os.path.join(case, "0", "ibmLinkRole")) else 0*a
    U = read_field(os.path.join(case, t, "U"))[0]; p = read_field(os.path.join(case, t, "p"))[0]
    Ub = read_field(os.path.join(case, t, "ibmUb"))[0]
    n = int(round(np.sqrt(len(a)))); h = L/n
    I = np.arange(len(a)); i = I % n; j = I // n
    nbmin = np.ones(len(a))
    for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        nbmin = np.minimum(nbmin, a[((j + dj) % n)*n + (i + di) % n])
    x = C[:, 0] - cx; y = C[:, 1] - cy; r = np.hypot(x, y)
    Uex, pex = exact(x, y)
    fl = a > 0
    eU = np.hypot(U[:, 0] - Uex[:, 0], U[:, 1] - Uex[:, 1])
    ep = p - pex; ep -= np.sum(ep[fl]*a[fl])/np.sum(a[fl])
    cls = np.full(len(a), "", dtype=object)
    cls[fl & (a >= 1)] = "bulk"
    cls[fl & (a >= 1) & (nbmin < 1)] = "layer1"
    cls[fl & (a < 1)] = "cut"
    cls[role == 2] = "master"
    cls[role == 1] = "slave"
    out = {}
    for wall, m in (("inner", r < 1.0), ("outer", r >= 1.0)):
        for k in CLASSES:
            s = m & (cls == k)
            if k == "bulk":
                s = s & (np.abs(r - 1.0) < 0.5 - 2*h)
            if not s.any():
                continue
            w = a[s]
            out[(wall, k)] = (s.sum(), np.sqrt(np.sum(w*eU[s]**2)/np.sum(w)),
                              np.sqrt(np.sum(w*ep[s]**2)/np.sum(w)), np.abs(ep[s]).max())
    blk = fl & (a >= 1) & (np.abs(r - 1.0) < 0.5 - 2*h)
    ut = (-y*U[:, 0] + x*U[:, 1])/r
    M = np.stack([r[blk], 1/r[blk]], 1)
    b = np.linalg.lstsq(M, ut[blk], rcond=None)[0][1]
    return n, out, b/B - 1


def main():
    for case in sys.argv[1:]:
        n, out, eT = analyse(case)
        print("%s  N=%d  torque error %+.3e" % (case, n, eT))
        print("    %-6s %-7s %6s %10s %10s %10s" % ("wall", "class", "cells", "U rms", "p rms", "p max"))
        for (wall, k), (cnt, u, prms, pmax) in out.items():
            print("    %-6s %-7s %6d %10.3e %10.3e %10.3e" % (wall, k, cnt, u, prms, pmax))


if __name__ == "__main__":
    main()
