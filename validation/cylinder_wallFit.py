#!/usr/bin/env python
"""Pressure and viscous drag of the periodic cylinder (tutorials/cylinderRe1, centre (5, 5),
radius 1, N x N cells on 10 x 10) from quadratic least-squares fits over FULL cells around each
wall point (post-processing prototype; doc/secondOrderStudy.md, "Chasing the cut-cell pressure").

    cylinder_wallFit.py nu radius case_N41 case_N81 ...

Each case needs 0/C, 0/ibmAlpha, 0/ibmSw (writeFields true, ibmSetGeometry,
postProcess -func writeCellCentres -time 0) and the solution in time directory 1 (p, U,
ibmCentroidOffset). The wall point of a wall cell is the closest point of the circle to its
fluid centroid; the fit (value and gradient at the wall point, inverse-distance weights) uses
the full cells (alpha = 1) within radius*h of it. Drag per unit depth: pressure sum p_w Sw_x,
viscous -nu sum ((G + G^T) Sw)_x, divided by the depth h.
"""
import sys, os, re
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from foamio import read_field

c0 = np.array([5.0, 5.0]); R = 1.0


def qfit(X, V, x0):
    d = X - x0
    M = np.array([np.ones(len(d)), d[:, 0], d[:, 1], d[:, 0]**2, d[:, 0]*d[:, 1], d[:, 1]**2]).T
    w = 1/np.maximum(np.hypot(d[:, 0], d[:, 1]), 1e-12)
    return np.linalg.lstsq(M*w[:, None], V*w[:, None], rcond=None)[0][:3]


nu, rad = float(sys.argv[1]), float(sys.argv[2])
rows = []
for c in sys.argv[3:]:
    N = int(re.search(r"_N(\d+)", c).group(1)); h = 10.0/N
    C, _ = read_field(os.path.join(c, "0", "C")); a, _ = read_field(os.path.join(c, "0", "ibmAlpha"))
    Sw, _ = read_field(os.path.join(c, "0", "ibmSw")); p, _ = read_field(os.path.join(c, "1", "p"))
    U, _ = read_field(os.path.join(c, "1", "U"))
    off = read_field(os.path.join(c, "1", "ibmCentroidOffset"))[0]
    X = C[:, :2] + off[:, :2]; full = a >= 1
    wc = np.where((a > 0) & (np.hypot(Sw[:, 0], Sw[:, 1]) > 0))[0]
    Xf = X[full]; Vf = np.column_stack([p[full], U[full, 0], U[full, 1]]); tree = cKDTree(Xf)
    Fp = Fv = 0.0
    for ci in wc:
        rr = X[ci] - c0; xw = c0 + R*rr/np.linalg.norm(rr); s = Sw[ci, :2]
        idx = tree.query_ball_point(xw, rad*h)
        f = qfit(Xf[idx], Vf[idx], xw)
        G = np.array([[f[1, 1], f[2, 1]], [f[1, 2], f[2, 2]]])      # G[i, j] = du_i/dx_j
        Fp += f[0, 0]*s[0]
        Fv += -nu*((G + G.T) @ s)[0]
    rows.append((N, Fp/h, Fv/h))
print("%6s %10s %10s %10s" % ("N", "pressure", "viscous", "total"))
for N, fp, fv in rows:
    print("%6d %10.5f %10.5f %10.5f" % (N, fp, fv, fp + fv))
for k, name in ((1, "pressure"), (2, "viscous")):
    d = np.diff([r[k] for r in rows])
    if len(d) > 1:
        print("%-8s orders %s" % (name, " ".join("%.2f" % x for x in np.log2(np.abs(d[:-1]/d[1:])))))
