#!/usr/bin/env python3
"""Numerical check of the two fixes proposed in the handout, on real 2-D cut-cell
geometry (cylinder of radius 1, Cartesian mesh, sub-sampled alpha, sliver absorption).

1. Mass source of a cylinder rotating about its axis, per wall cell:
   centre u_b(x_P).Sw, movingBody/scalar u_b(x_P + d n_w).Sw, scalar_secondOrder
   u_b(x_c + d n_w).Sw, proposed -sum_f theta_f u_b(x'_f).S_f.
2. Cut-cell pressure gradient of a linear pressure, relative error against the
   exact integral V_geo G over the fluid polygon: first order, current second order,
   proposed.
"""
import sys
import numpy as np

sys.path.insert(0, ".")
import ccgeom as cg

R = 1.0
body = cg.Circle(0.0, 0.0, R)


def polygon(g, i, j):
    c = [g.corner(i, j), g.corner(i + 1, j), g.corner(i + 1, j + 1), g.corner(i, j + 1)]
    pts = []
    for k in range(4):
        a, b = c[k], c[(k + 1) % 4]
        if body(*a) > 0:
            pts.append(a)
        if (body(*a) > 0) != (body(*b) > 0):
            pts.append(cg.bisect(body, a, b) if body(*a) > 0 else cg.bisect(body, b, a))
    return np.array(pts)


def area_centroid(P):
    x, y = P[:, 0], P[:, 1]
    xs, ys = np.roll(x, -1), np.roll(y, -1)
    cr = x * ys - xs * y
    A = cr.sum() / 2
    return A, np.array([((x + xs) * cr).sum(), ((y + ys) * cr).sum()]) / (6 * A)


def on_boundary(g, i, j, p, tol=1e-12):
    x0, y0 = g.corner(i, j)
    return (abs(p[0] - x0) < tol or abs(p[0] - x0 - g.h) < tol or
            abs(p[1] - y0) < tol or abs(p[1] - y0 - g.h) < tol)


def faces(g, i, j):
    """(theta-array, index, outward area vector, end points) of the 4 faces."""
    h = g.h
    return [("x", (i, j), np.array([-h, 0.0]), (g.corner(i, j), g.corner(i, j + 1)), (i - 1, j)),
            ("x", (i + 1, j), np.array([h, 0.0]), (g.corner(i + 1, j), g.corner(i + 1, j + 1)), (i + 1, j)),
            ("y", (i, j), np.array([0.0, -h]), (g.corner(i, j), g.corner(i + 1, j)), (i, j - 1)),
            ("y", (i, j + 1), np.array([0.0, h]), (g.corner(i, j + 1), g.corner(i + 1, j + 1)), (i, j + 1))]


def analyse(h, ox, oy, ns=32):
    n = int(np.ceil(3.0 / h)) + 2
    g = cg.Grid(-1.5 + ox, -1.5 + oy, h, n, n)
    thx0, thy0, a0, _, _ = cg.fractions(g, body, ns=ns, absorb=False)
    thx, thy, alpha, cen, absd = cg.fractions(g, body, ns=ns, absorb=True)
    th0 = {"x": thx0, "y": thy0}
    th = {"x": thx, "y": thy}
    G = np.array([0.7, -0.4])                       # linear pressure p = G.x
    p = lambda x: G @ x
    rows = []
    for i in range(1, n - 1):
        for j in range(1, n - 1):
            if alpha[i, j] <= 0:
                continue
            Sw = cg.wall_vector(g, thx, thy, i, j)
            Aw = np.hypot(*Sw)
            if Aw < 1e-10 * h:
                continue
            nw = Sw / Aw
            xP = g.centre(i, j)
            xc = cen[i, j]                            # sub-sampled centroid (as the code)
            dwall = alpha[i, j] * h * h / (2 * Aw)
            dc = max(body(*xc), 1e-3 * h)
            ub = lambda x: np.array([-x[1], x[0]])    # omega = 1 about the body axis
            # ---- 1. mass source
            S_centre = ub(xP) @ Sw
            S_mb = ub(xP + dwall * nw) @ Sw
            S_so = ub(xc + dwall * nw) @ Sw
            S_new = 0.0
            # ---- 2. gradient
            G1 = p(xc) * Sw                          # first order: p_P Sw + sum theta p_f Sf
            G2 = p(xc) * Sw + (G @ nw) * dc * Sw     # current second order, wall part
            G3 = np.zeros(2)                          # proposed, faces
            wall_pieces = []                          # proposed, wall: (vector, centroid)
            for kind, idx, S, (e0, e1), nb in faces(g, i, j):
                t, t0 = th[kind][idx], th0[kind][idx]
                if t > 0:
                    fr, seg = cg.wet_segment(body, e0, e1)
                    xw = 0.5 * (seg[0] + seg[1])     # wet-face centroid
                    xf = 0.5 * (e0 + e1)
                    cN = cen[nb]
                    pf = 0.5 * (p(xc) + p(cN))       # linear interpolate: value at centroid midpoint
                    G1 += t * pf * S
                    corr = G @ (xw - xf) if 0 < t < 1 else 0.0
                    G2 += t * (pf + corr) * S
                    G3 += t * (pf + G @ (xw - 0.5 * (xc + cN))) * S
                    S_new -= t * (ub(xw) @ S)
                else:                                 # geometric wet part closed by the post-processing
                    fr, seg = cg.wet_segment(body, e0, e1)
                    if fr > 0:
                        wall_pieces.append((fr * S, 0.5 * (seg[0] + seg[1])))
            P = polygon(g, i, j)
            chords = [(P[k], P[(k + 1) % len(P)]) for k in range(len(P))
                      if not on_boundary(g, i, j, 0.5 * (P[k] + P[(k + 1) % len(P)]))]
            for a, b in chords:                       # chord: outward normal into the body
                t = b - a
                wall_pieces.append((np.array([t[1], -t[0]]), 0.5 * (a + b)))
            Ssum = sum(v for v, _ in wall_pieces) if wall_pieces else np.zeros(2)
            if np.hypot(*(Ssum - Sw)) > 1e-12 * h:    # orientation check of the chord normal
                wall_pieces = [(-v if k >= len(wall_pieces) - len(chords) else v, x)
                               for k, (v, x) in enumerate(wall_pieces)]
                Ssum = sum(v for v, _ in wall_pieces)
            assert np.hypot(*(Ssum - Sw)) < 1e-10 * h, (i, j, Ssum, Sw)
            for v, x in wall_pieces:
                G3 += (p(xc) + G @ (x - xc)) * v
            Vgeo, _ = area_centroid(P)
            ref = Vgeo * G
            err = [np.hypot(*(Gx - ref)) / np.hypot(*ref) for Gx in (G1, G2, G3)]
            phi = np.degrees(np.arctan2(*(xP[::-1])))
            rows.append([phi, alpha[i, j], S_centre, S_mb, S_so, S_new, Aw, *err])
    return np.array(rows)


if __name__ == "__main__":
    r = analyse(0.1, -0.013, -0.037)
    print("cells", len(r))
    norm = 0.1 * r[:, 6]                              # |omega| h A_wall
    for k, lab in zip(range(2, 6), ["centre", "x_P + d n", "x_c + d n", "proposed"]):
        print("%-10s max|S|/(w h Aw) = %.3e   sum S = %.3e" % (lab, np.abs(r[:, k] / norm).max(), r[:, k].sum()))
    print("S centre vs x_P+d n: max diff %.2e" % np.abs(r[:, 2] - r[:, 3]).max())
    for k, lab in zip(range(7, 10), ["first order", "current 2nd", "proposed"]):
        print("%-12s grad rel. error: max %.3e median %.3e" % (lab, r[:, k].max(), np.median(r[:, k])))
    np.save("analysis_h01.npy", r)
