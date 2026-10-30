#!/usr/bin/env python
"""Surface pressure and skin friction of NACA 0012 runs (body-fitted and IBM).

    naca0012_surface.py [--plot out.png] [--window w] [--image D] case1 case2 ...

Body-fitted: C_p from the wall-adjacent cells (p is zeroGradient on the wall),
C_f from the wallShearStress function object, at the wall-face centres.
IBM: C_p from the pressure of the wall cells, C_f from Spalding's law with the
tangential slip of the wall cell at its wall distance (the wall shear the
solver imposes; 0/ibmWallDistance, which includes linking), at the closest
points of the airfoil to the fluid centroids. C_f carries the sign of the
streamwise wall shear (negative where the flow is reversed). Upper and lower
surfaces from the sign of y.

These cell values follow the staircase of the cut cells (wall distances from
almost 0 to 1.5 h, merged groups): a sawtooth in C_f and small jumps in C_p,
although the field itself is smooth.
  --window w   area-weighted moving average over an arc length w of the wall
               shear and pressure the solver applies (preserves the
               integrated force; w = 4 h removes the sawtooth)
  --image D    image-point estimate instead (both kinds of runs): U and p
               interpolated (linear, over cell centres or fluid centroids) at
               the distance D along the normal from the airfoil, C_f from
               Spalding's law at D, C_p = 2 p(D) (dp/dn = 0 across the
               boundary layer). Smooth, but a different quantity: it assumes
               the wall law up to D (about 5 % high at y+ 50 on a
               wall-resolved SA solution, more further out).
"""
import argparse, math, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from foamio import read_field, latest_time

NU = 1e-6


def spalding_uplus(Rey):
    k, E = 0.41, 9.8
    if Rey <= 0:
        return 0.0
    g = lambda u: u + (math.exp(min(k * u, 50)) - 1 - k * u - (k * u)**2 / 2 - (k * u)**3 / 6) / E - Rey / u
    lo, hi = 1e-12, math.sqrt(Rey)
    if g(hi) <= 0:
        return hi
    for _ in range(100):
        m = 0.5 * (lo + hi)
        if g(m) > 0: hi = m
        else: lo = m
    return 0.5 * (lo + hi)


def owner_of_patch(case, patch):
    d = os.path.join(case, "constant", "polyMesh")
    b = open(os.path.join(d, "boundary")).read()
    m = re.search(patch + r"\s*\{[^}]*nFaces\s+(\d+);\s*startFace\s+(\d+);", b)
    n, s = int(m.group(1)), int(m.group(2))
    txt = open(os.path.join(d, "owner")).read()
    txt = txt[txt.index("(", txt.index("}")) + 1:]
    own = np.array(txt[:txt.index(")")].split(), int)
    return own[s:s + n]


def yt(x):
    x = np.asarray(x, float)
    return 0.594689181 * (0.298222773 * np.sqrt(x) - 0.127125232 * x - 0.357907906 * x**2
                          + 0.291984971 * x**3 - 0.105174606 * x**4)


def airfoil(n=4001):
    """Dense airfoil points (cosine spacing), their outward unit normals and
    arc length, from the lower trailing edge round the nose to the upper one."""
    xs = 0.5 * (1 - np.cos(np.linspace(0, math.pi, n)))
    x = np.concatenate([xs[::-1], xs[1:]])
    y = np.concatenate([-yt(xs[::-1]), yt(xs[1:])])
    P = np.column_stack([x, y])
    t = np.gradient(P, axis=0)
    t /= np.linalg.norm(t, axis=1)[:, None]
    N = np.column_stack([t[:, 1], -t[:, 0]])          # outward for this orientation
    N *= np.where(np.einsum("ij,ij->i", N, P - [0.5, 0]) < 0, -1, 1)[:, None]
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])
    return P, N, s


def windowed(s, area, *fields, w):
    """Area-weighted moving average over the arc length w (fields on sorted s)."""
    out = []
    cs = np.concatenate([[0], np.cumsum(area)])
    lo = np.searchsorted(s, s - 0.5 * w, "left")
    hi = np.searchsorted(s, s + 0.5 * w, "right")
    for f in fields:
        cf = np.concatenate([[0], np.cumsum(area * f)])
        out.append((cf[hi] - cf[lo]) / (cs[hi] - cs[lo]))
    return out


def image(case, t, D, ibm):
    """Image-point C_f and C_p at the dense airfoil points."""
    from scipy.interpolate import LinearNDInterpolator
    C, _ = read_field(os.path.join(case, "0", "C"))
    if ibm:
        a, _ = read_field(os.path.join(case, "0", "ibmAlpha"))
        off = os.path.join(case, t, "ibmCentroidOffset")
        C = C + (read_field(off)[0] if os.path.exists(off) else 0)
    else:
        a = np.ones(len(C))
    p, _ = read_field(os.path.join(case, t, "p"))
    U, _ = read_field(os.path.join(case, t, "U"))
    box = (a > 0) & (C[:, 0] > -0.1) & (C[:, 0] < 1.1) & (abs(C[:, 1]) < 0.2)
    ip = LinearNDInterpolator(C[box, :2], np.column_stack([U[box, :2], p[box]]))
    P, N, s = airfoil()
    keep = (P[:, 0] > 1e-4) & (P[:, 0] < 1 - 1e-4)
    P, N = P[keep], N[keep]
    v = ip(P + D * N)
    tx = np.column_stack([N[:, 1], -N[:, 0]])
    tx *= np.sign(tx[:, 0:1] + 1e-30)                   # along +x
    ut = np.einsum("ij,ij->i", v[:, :2], tx)
    cf = np.array([2 * np.sign(u) * (abs(u) / spalding_uplus(abs(u) * D / NU))**2
                   if np.isfinite(u) and u != 0 else np.nan for u in ut])
    return P[:, 0], P[:, 1], 2 * v[:, 2], cf


def surface(case, window=0.0, D=0.0):
    t = latest_time(case)
    ibm = os.path.exists(os.path.join(case, "constant", "ibmProperties"))
    if D > 0:
        return image(case, t, D, ibm)
    if ibm:
        C, _ = read_field(os.path.join(case, "0", "C"))
        a, _ = read_field(os.path.join(case, "0", "ibmAlpha"))
        Sw, _ = read_field(os.path.join(case, "0", "ibmSw"))
        dw, _ = read_field(os.path.join(case, "0", "ibmWallDistance"))
        p, _ = read_field(os.path.join(case, t, "p"))
        U, _ = read_field(os.path.join(case, t, "U"))
        off = os.path.join(case, t, "ibmCentroidOffset")
        X = C + (read_field(off)[0] if os.path.exists(off) else 0)
        Aw = np.linalg.norm(Sw, axis=1)
        wall = np.where((a > 0) & (Aw > 1e-10 * 1e-4) & (dw < 1))[0]
        P, _, sP = airfoil()
        from scipy.spatial import cKDTree
        _, k = cKDTree(P).query(X[wall, :2])
        cp, cf = np.empty(len(wall)), np.empty(len(wall))
        for j, i in enumerate(wall):
            n = Sw[i] / Aw[i]
            ut = U[i] - (U[i] @ n) * n
            mag = np.linalg.norm(ut)
            up = spalding_uplus(mag * dw[i] / NU)
            utau2 = (mag / up)**2 if up > 0 else 0.0
            # streamwise sign: along +x on both surfaces for attached flow
            sgn = np.sign(ut[0]) if mag > 0 else 1.0
            cp[j], cf[j] = 2 * p[i], 2 * sgn * utau2
        o = np.argsort(sP[k])
        x, y, s, area, cp, cf = P[k[o], 0], P[k[o], 1], sP[k[o]], Aw[wall[o]], cp[o], cf[o]
    else:
        C, Cb = read_field(os.path.join(case, "0", "C"))
        p, _ = read_field(os.path.join(case, t, "p"))
        _, tw = read_field(os.path.join(case, t, "wallShearStress"))
        own = owner_of_patch(case, "airfoil")
        xf = Cb["airfoil"]
        # wallShearStress is the stress on the wall from the fluid, -nu_eff dU/dn:
        # the streamwise wall shear is -tau_x
        x, y, cp, cf = xf[:, 0], xf[:, 1], 2 * p[own], -2 * tw["airfoil"][:, 0]
        if window <= 0:
            return x, y, cp, cf
        P, _, sP = airfoil()
        from scipy.spatial import cKDTree
        _, k = cKDTree(P).query(xf[:, :2])
        o = np.argsort(sP[k])
        x, y, s, cp, cf = x[o], y[o], sP[k[o]], cp[o], cf[o]
        area = np.ones(len(x))
    if window > 0:
        cp, cf = windowed(s, area, cp, cf, w=window)
    return x, y, cp, cf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plot")
    ap.add_argument("--window", type=float, default=0.0, help="arc-length averaging window")
    ap.add_argument("--image", type=float, default=0.0, help="image-point distance")
    ap.add_argument("cases", nargs="+")
    a = ap.parse_args()
    data = []
    for c in a.cases:
        x, y, cp, cf = surface(c, a.window, a.image)
        ok = np.isfinite(cp) & np.isfinite(cf)
        data.append((os.path.basename(os.path.normpath(c)), (x[ok], y[ok], cp[ok], cf[ok])))
    ref = data[0][1]
    for name, (x, y, cp, cf) in data:
        line = "%-20s wall points %4d  min Cp %.3f" % (name, len(x), cp.min())
        # sawtooth: rms deviation of C_f from its mean over 0.05 c, 0.2 < x < 0.8
        m = (x > 0.2) & (x < 0.8) & (y >= 0)
        o = np.argsort(x[m]); xm, cm = x[m][o], cf[m][o]
        sm = np.array([cm[abs(xm - xi) < 0.025].mean() for xi in xm])
        line += "  Cf noise %.3f" % (np.sqrt(np.mean((cm - sm)**2)) / np.mean(cm))
        for side, m, mr in (("upper", y >= 0, ref[1] >= 0), ("lower", y < 0, ref[1] < 0)):
            o, orf = np.argsort(x[m]), np.argsort(ref[0][mr])
            xs = np.linspace(0.05, 0.95, 19)
            cfi = np.interp(xs, x[m][o], cf[m][o]); cfr = np.interp(xs, ref[0][mr][orf], ref[3][mr][orf])
            cpi = np.interp(xs, x[m][o], cp[m][o]); cpr = np.interp(xs, ref[0][mr][orf], ref[2][mr][orf])
            line += "  %s: Cf(0.5) %.5f, rms dCf/Cf %.3f, rms dCp %.4f" % (
                side, np.interp(0.5, x[m][o], cf[m][o]), np.sqrt(np.mean(((cfi - cfr) / cfr)**2)),
                np.sqrt(np.mean((cpi - cpr)**2)))
        print(line)
    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        for name, (x, y, cp, cf) in data:
            for m, ls in ((y >= 0, "-"), (y < 0, "--")):
                o = np.argsort(x[m])
                l, = ax[0].plot(x[m][o], -cp[m][o], ls, lw=1, label=name if ls == "-" else None)
                ax[1].plot(x[m][o], cf[m][o], ls, lw=1, color=l.get_color())
        ax[0].set_xlabel("x/c"); ax[0].set_ylabel("$-C_p$"); ax[0].legend(fontsize=7, frameon=False)
        ax[1].set_xlabel("x/c"); ax[1].set_ylabel("$C_f$ (upper solid, lower dashed)")
        ax[1].set_ylim(-0.002, 0.012); ax[1].set_xlim(0, 1)
        fig.tight_layout(); fig.savefig(a.plot, dpi=120)
        print("wrote", a.plot)


if __name__ == "__main__":
    main()
