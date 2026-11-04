#!/usr/bin/env python
"""Viscous drag of NACA 0012 runs by chordwise bin and, for IBM runs, by wall-cell class.

    naca0012_wallClasses.py [--time T] [--ref bfcase] case1 case2 ...

IBM: the wall shear the solver applies, <time>/ibmWallShear (nu_w |u_t|/dWall), on the wall
area |Sw| of each wall cell, along its tangential velocity; the drag part of it in 1e-4 C_d.
Classes from 0/ibmLinkRole (if present): slave (1), master (2), cut (unlinked cut cell),
full (alpha = 1 wall cell behind a closed face). Also the area-weighted means per class of
the wall-law weight (<time>/ibmWallLawWeight), of d+ = u_tau dWall/nu (u_tau from the
applied shear, dWall = 0/ibmWallDistance) and of h+ = u_tau h/nu (h from the cell volume).
Body-fitted (--ref): -wallShearStress on the airfoil patch, face area from the arc length.
"""
import argparse, math, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from foamio import read_field, latest_time
from naca0012_surface import airfoil, owner_of_patch

NU, DZ = 1e-6, 0.01
BINS = [0, .002, .01, .05, .2, .5, .8, .95, .99, 1.0001]


def drag_dir(case):
    txt = open(os.path.join(case, "0", "U")).read()
    m = re.search(r"internalField\s+uniform\s+\(([^)]*)\)", txt)
    if not m:
        txt = open(os.path.join(case, "0.orig", "U")).read()
        m = re.search(r"internalField\s+uniform\s+\(([^)]*)\)", txt)
    u = np.array(m.group(1).split(), float)[:2]
    return u/np.linalg.norm(u), np.linalg.norm(u)


def ibm_cells(case, t):
    r0 = lambda f: read_field(os.path.join(case, "0", f))[0]
    C, a, Sw, dw = r0("C"), r0("ibmAlpha"), r0("ibmSw"), r0("ibmWallDistance")
    U = read_field(os.path.join(case, t, "U"))[0]
    tau = read_field(os.path.join(case, t, "ibmWallShear"))[0]
    fw = os.path.join(case, t, "ibmWallLawWeight")
    w = read_field(fw)[0] if os.path.exists(fw) else 0.0
    w = np.broadcast_to(np.asarray(w, float), len(a))      # uniform without a wall law
    role = r0("ibmLinkRole") if os.path.exists(os.path.join(case, "0", "ibmLinkRole")) else 0*a
    Aw = np.linalg.norm(Sw, axis=1)
    wall = np.where((a > 0) & (Aw > 1e-14) & (tau > 0))[0]
    n = Sw[wall]/Aw[wall, None]
    u = U[wall]
    ut = u - np.sum(u*n, 1)[:, None]*n
    that = ut/np.maximum(np.linalg.norm(ut, axis=1), 1e-30)[:, None]
    cls = np.where(a[wall] < 1, "cut", "full").astype(object)
    cls[role[wall] == 2] = "master"
    cls[role[wall] == 1] = "slave"
    # cell size from the cut-cell candidates: h = sqrt(V/DZ) with V from alpha and the
    # fraction field is not written; use the distance between neighbouring wall centroids
    return dict(C=C[wall], Aw=Aw[wall], tau=tau[wall], t=that[:, :2], w=w[wall], d=dw[wall],
                cls=cls, a=a[wall])


def chord_x(XY):
    from scipy.spatial import cKDTree
    P, _, _ = airfoil()
    _, k = cKDTree(P).query(XY)
    return P[k, 0], np.sign(P[k, 1] + 1e-30)


def bf_cells(case, t):
    Cb = read_field(os.path.join(case, "0", "C"))[1]
    _, tw = read_field(os.path.join(case, t, "wallShearStress"))
    xf = Cb["airfoil"][:, :2]
    tau = -tw["airfoil"][:, :2]                       # shear on the fluid side, along the flow
    P, _, sP = airfoil()
    from scipy.spatial import cKDTree
    _, k = cKDTree(P).query(xf)
    s = sP[k]
    o = np.argsort(s)
    ds = np.empty(len(s))
    so = s[o]
    mid = 0.5*(so[1:] + so[:-1])
    ds[o] = np.diff(np.concatenate([[so[0] - (mid[0] - so[0])], mid, [so[-1] + (so[-1] - mid[-1])]]))
    return xf, tau, ds*DZ


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref")
    ap.add_argument("--time")
    ap.add_argument("cases", nargs="+")
    a = ap.parse_args()
    hdr = "".join("%7s" % ("<%g" % b) for b in BINS[1:]) + "   total"
    rows = []
    if a.ref:
        t = latest_time(a.ref)
        e, U0 = drag_dir(a.ref)
        xf, tau, A = bf_cells(a.ref, t)
        x, _ = chord_x(xf)
        cd = (tau @ e)*A/(0.5*U0**2*DZ)*1e4
        ref = np.histogram(x, BINS, weights=cd)[0]
        rows.append(("ref " + os.path.basename(a.ref.rstrip("/")), ref))
    print("viscous drag by x bin, 1e-4 C_d (upper + lower)")
    print("%-34s" % "" + hdr)
    if a.ref:
        print("%-34s" % rows[0][0] + "".join("%7.2f" % v for v in ref) + "  %7.2f" % ref.sum())
    for case in a.cases:
        t = a.time or latest_time(case)
        e, U0 = drag_dir(case)
        c = ibm_cells(case, t)
        x, side = chord_x(c["C"][:, :2])
        cd = c["tau"]*c["Aw"]*(c["t"] @ e)/(0.5*U0**2*DZ)*1e4
        name = "/".join(os.path.abspath(case).split("/")[-2:])
        tot = np.histogram(x, BINS, weights=cd)[0]
        line = "%-34s" % name + "".join("%7.2f" % v for v in tot) + "  %7.2f" % tot.sum()
        if a.ref:
            line += "   - ref: " + " ".join("%+.2f" % v for v in tot - ref) + "  %+.2f" % (tot.sum() - ref.sum())
        print(line)
        for k in ("cut", "master", "slave", "full"):
            m = c["cls"] == k
            if not m.any():
                continue
            h = np.histogram(x[m], BINS, weights=cd[m])[0]
            ut = np.sqrt(c["tau"][m])
            A = c["Aw"][m]
            print("    %-8s n=%6d  w %.3f  d+ %5.2f" % (k, m.sum(), np.sum(A*c["w"][m])/A.sum(),
                  np.sum(A*ut*c["d"][m]/NU)/A.sum()) + "".join("%7.2f" % v for v in h) + "  %7.2f" % h.sum())


if __name__ == "__main__":
    main()
