#!/usr/bin/env python3
"""Prepare a NASA TMR NACA 0012 C-grid (2-D PLOT3D, .p2dfmt) for plot3dToFoam.

    tmrGrid.py <in.p2dfmt> <out.p2dfmt> [--firstCell h]

The grids (e.g. n0012_449-129.p2dfmt from https://tmbwg.github.io/turbmodels/
naca0012_grids.html) are wall-resolved for Re = 6e6. With --firstCell h the
grid lines next to the wall are removed until the first cell at mid-chord is
at least h high (a wall-function mesh with the same outer grid). The output is
written with the i direction reversed if needed so that plot3dToFoam -2D
produces cells of positive volume.
"""
import argparse
import numpy as np


def read(path):
    tok = open(path).read().split()
    nb = int(tok[0])
    if nb != 1:
        raise SystemExit("expected a single block")
    ni, nj = int(tok[1]), int(tok[2])
    v = np.array(tok[3:3 + 2 * ni * nj], float)
    return v[:ni * nj].reshape(nj, ni), v[ni * nj:].reshape(nj, ni)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp")
    ap.add_argument("out")
    ap.add_argument("--firstCell", type=float, default=0.0)
    a = ap.parse_args()
    x, y = read(a.inp)
    nj, ni = x.shape

    if a.firstCell > 0:
        wall = np.where((x[0] >= -1e-9) & (x[0] <= 1 + 1e-9) & (y[0] >= 0))[0]
        im = wall[np.argmin(np.abs(x[0, wall] - 0.5))]
        d = np.hypot(x[:, im] - x[0, im], y[:, im] - y[0, im])
        k = int(np.argmax(d >= a.firstCell))
        keep = [0] + list(range(k, nj))
        x, y = x[keep], y[keep]
        print("removed %d wall lines: first cell at mid-chord %.3g" % (k - 1, d[k]))

    # orientation: cell (0,0) must be anticlockwise in (i, j) for positive volumes
    def area(x, y):
        return 0.5 * ((x[0, 1] - x[0, 0]) * (y[1, 0] - y[0, 0]) - (y[0, 1] - y[0, 0]) * (x[1, 0] - x[0, 0]))
    if area(x, y) < 0:
        x, y = x[:, ::-1], y[:, ::-1]
    nj, ni = x.shape
    with open(a.out, "w") as f:
        f.write("1\n%d %d\n" % (ni, nj))
        for arr in (x, y):
            vals = arr.ravel()
            for k in range(0, len(vals), 4):
                f.write(" ".join("%.15e" % v for v in vals[k:k + 4]) + "\n")
    print("wrote %s: %d x %d" % (a.out, ni, nj))


if __name__ == "__main__":
    main()
