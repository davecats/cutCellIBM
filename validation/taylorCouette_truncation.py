"""Cut-cell pressure of Taylor-Couette against the local truncation error of the exact solution
(doc/secondOrderStudy.md, "The cut-cell pressure is a continuity error").

    taylorCouette_truncation.py classes  res_case ...
        truncation errors written by ibmTruncationTC, multiplied by alpha (per cell volume V:
        the force or volume flux a cell sees), rms by class: inner/outer wall cut cells with
        2, 3 or 4 open neighbour faces; continuity split into face-flux and wall-source parts
    taylorCouette_truncation.py fit  solver_case res_case
        solver pressure error of the inner cut cells (bulk mean removed, exact pressure at the
        fluid centroids Xvalue of res_case) fitted cell by cell with the continuity error
        alpha*Rcont and the viscous error alpha*h*Rvisc (tangential, radial)

res_case: a case of the same mesh and ibmProperties after ibmTruncationTC (time 0.02).
Uniform 2-D mesh of the 3.2 x 3.2 box, cells numbered row by row.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from foamio import read_field, latest_time

Ri, Ro, cx, cy = 0.5, 1.5, 1.6, 1.6
A = -Ri**2/(Ro**2 - Ri**2); B = Ri**2*Ro**2/(Ro**2 - Ri**2)


def pex(r):
    return A*A*r*r/2 + 2*A*B*np.log(r) - B*B/(2*r*r)


def res(case):
    rd = lambda f: read_field(os.path.join(case, "0.02", f))[0]
    d = {k: rd(k) for k in ("Xvalue", "ibmAlpha", "Rcont", "RcontFace", "RcontWall", "Rvisc")}
    n = int(round(np.sqrt(len(d["Xvalue"])))); h = 3.2/n
    a = d["ibmAlpha"]; I = np.arange(len(a)); i = I % n; j = I // n
    nopen = np.zeros(len(a), int)
    for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        nopen += a[np.clip(j+dj, 0, n-1)*n + np.clip(i+di, 0, n-1)] > 0
    x = d["Xvalue"][:, 0] - cx; y = d["Xvalue"][:, 1] - cy
    return d, n, h, nopen, x, y, np.hypot(x, y)


def classes(case):
    d, n, h, nopen, x, y, r = res(case)
    a = d["ibmAlpha"]
    s = "%-24s N%4d" % (case.rstrip("/").split("/")[-1], n)
    for wall, near in (("i", r < 1), ("o", r > 1)):
        for no in (2, 3, 4):
            m = (a > 0) & (a < 1) & near & (nopen == no)
            rms = lambda k: np.sqrt(np.mean((a*d[k])[m]**2))
            s += " | %s%d cont %7.1e face %7.1e wall %7.1e" % (wall, no, rms("Rcont"), rms("RcontFace"), rms("RcontWall"))
    print(s)


def fit(sc, rc):
    d, n, h, nopen, x, y, r = res(rc)
    a = d["ibmAlpha"]
    t = latest_time(sc)
    p, _ = read_field(os.path.join(sc, t, "p"))
    e = p - pex(r)
    bulk = (a >= 1) & (np.minimum(r - Ri, Ro - r) > 3*h); e -= e[bulk].mean()
    cut = (a > 0) & (a < 1) & (r < 1)
    q = (a*d["Rcont"])[cut]
    Rv = d["Rvisc"]
    vt = (a*h*(-Rv[:, 0]*y + Rv[:, 1]*x)/r)[cut]; vr = (a*h*(Rv[:, 0]*x + Rv[:, 1]*y)/r)[cut]
    for name, M in (("cont", np.c_[q]), ("cont+visc", np.c_[q, vt, vr])):
        c, *_ = np.linalg.lstsq(M, e[cut], rcond=None)
        rr = e[cut] - M @ c
        print("%-24s N%4d %-10s coef %s  rms %8.2e -> %8.2e (explained %4.1f %%)" % (
            sc.rstrip("/").split("/")[-1], n, name, " ".join("%+7.3f" % v for v in c),
            np.sqrt(np.mean(e[cut]**2)), np.sqrt(np.mean(rr**2)), 100*(1 - np.mean(rr**2)/np.mean(e[cut]**2))))


if __name__ == "__main__":
    if sys.argv[1] == "classes":
        for c in sys.argv[2:]:
            classes(c)
    else:
        fit(*sys.argv[2:])
