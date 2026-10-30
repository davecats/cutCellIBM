#!/usr/bin/env python
"""Compare Spalart-Allmaras channel runs (tutorials/turbulence/channelSA/makeCase.py).

    channelSA_compare.py [--plot out.png] case1 case2 ...

For each case: u_tau from the force balance of the mean-velocity force,
u_tau^2 = G delta (delta = 1), Re_tau, and the profiles u+(y+) and nut/nu at
the cell values (fluid centroids in cut cells), y the distance to the nearer
wall at y = -1 or 1. Prints u_tau and the relative difference of each case to
the first one, and u+ at a few y+ (interpolated).
"""
import argparse, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from foamio import read_field, latest_time

NU = 1e-4


def load(case):
    log = [f for f in os.listdir(case) if f.startswith("log.") and "Foam" in f][0]
    G = float(re.findall(r"pressure gradient = (\S+)", open(os.path.join(case, log)).read())[-1])
    t = latest_time(case)
    C, _ = read_field(os.path.join(case, "0", "C"))
    U, _ = read_field(os.path.join(case, t, "U"))
    nut, _ = read_field(os.path.join(case, t, "nut"))
    y = C[:, 1].copy()
    a = np.ones(len(y))
    if os.path.exists(os.path.join(case, "0", "ibmAlpha")):
        a, _ = read_field(os.path.join(case, "0", "ibmAlpha"))
    fl = a > 0
    utau = np.sqrt(G * 1.0)
    d = 1 - np.abs(y)
    yw = os.path.join(case, "0", "yWall")
    if os.path.exists(yw):          # IBM: distance of the fluid centroid to the wall
        d, _ = read_field(yw)
    return dict(name=os.path.basename(os.path.normpath(case)), utau=utau, G=G,
                yp=(d * utau / NU)[fl], up=(U[fl, 0] / utau), nut=(nut[fl] / NU),
                alpha=a[fl], bottom=(y[fl] < 0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plot")
    ap.add_argument("cases", nargs="+")
    a = ap.parse_args()
    runs = [load(c) for c in a.cases]
    ref = runs[0]
    for r in runs:
        line = "%-24s u_tau %.5f  Re_tau %.1f  d(u_tau) %+.3f%%" % (
            r["name"], r["utau"], r["utau"] / NU, 100 * (r["utau"] / ref["utau"] - 1))
        for side, m in (("bottom", r["bottom"]), ("top", ~r["bottom"])):
            o = np.argsort(r["yp"][m])
            yp, up = r["yp"][m][o], r["up"][m][o]
            first = "first cell y+ %.2f (alpha %.2f)" % (yp[0], r["alpha"][m][o][0])
            vals = " ".join("u+(%g)=%.3f" % (q, np.interp(q, yp, up)) for q in (5, 30, 100, 300)
                            if yp[0] <= q <= yp[-1])
            line += "\n    %-6s %s  %s  u+ max %.3f" % (side, first, vals, up.max())
        print(line)
    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
        for r in runs:
            for side, m, mk in (("bottom", r["bottom"], "o"), ("top", ~r["bottom"], "s")):
                o = np.argsort(r["yp"][m])
                lab = r["name"] + ("" if side == "bottom" else " (top)")
                ax[0].semilogx(r["yp"][m][o], r["up"][m][o], "-" + mk, ms=3, lw=1, label=lab)
                ax[1].semilogx(r["yp"][m][o], r["nut"][m][o], "-" + mk, ms=3, lw=1, label=lab)
        yp = np.logspace(-1, 3, 200)
        ax[0].semilogx(yp, yp, "k:", lw=0.7)
        ax[0].semilogx(yp, np.log(yp) / 0.41 + 5.0, "k--", lw=0.7)
        ax[0].set_ylim(0, 30); ax[0].set_xlabel("$y^+$"); ax[0].set_ylabel("$u^+$")
        ax[1].set_xlabel("$y^+$"); ax[1].set_ylabel(r"$\nu_t/\nu$")
        ax[0].legend(fontsize=6, frameon=False)
        fig.tight_layout(); fig.savefig(a.plot, dpi=130)
        print("wrote", a.plot)


if __name__ == "__main__":
    main()
