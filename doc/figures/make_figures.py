#!/usr/bin/env python3
"""Figures of doc/cutCellIBM.tex.

    python make_figures.py [name ...]      (all figures when no name is given)

Schematics are drawn from real cut-cell data computed by ccgeom.py with the
algorithm of the library (bisection, sub-sampling, sliver absorption, closure).
"""
import sys
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, FancyArrowPatch, Circle as MCircle

import ccgeom as cg

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "font.size": 9,
    "axes.titlesize": 9,
    "axes.labelsize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})

C_SOLID = "0.80"
C_CUT = "#cfe3f3"
C_WET = "#1b9e3e"
C_WALL = "#d62728"
C_BODY = "#303030"
C_ARROW = "#1f4e9c"
C_GREY = "0.45"

FIGS = {}


def figure(fn):
    FIGS[fn.__name__[4:]] = fn
    return fn


def arrow(ax, p, q, color=C_ARROW, lw=1.2, ms=8, **kw):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=ms,
                                 color=color, lw=lw, shrinkA=0, shrinkB=0, **kw))


def cell_polygon(g, phi, i, j):
    """Fluid polygon of cell (i, j): fluid corners and crossings, in order."""
    c = [g.corner(i, j), g.corner(i + 1, j), g.corner(i + 1, j + 1), g.corner(i, j + 1)]
    pts = []
    for k in range(4):
        a, b = c[k], c[(k + 1) % 4]
        if phi(*a) > 0:
            pts.append(a)
        if (phi(*a) > 0) != (phi(*b) > 0):
            pts.append(cg.bisect(phi, a, b) if phi(*a) > 0 else cg.bisect(phi, b, a))
    return np.array(pts)


def draw_grid(ax, g, color="0.6", lw=0.5):
    for i in range(g.nx + 1):
        x = g.x0 + i * g.h
        ax.plot([x, x], [g.y0, g.y0 + g.ny * g.h], color=color, lw=lw, zorder=1)
    for j in range(g.ny + 1):
        y = g.y0 + j * g.h
        ax.plot([g.x0, g.x0 + g.nx * g.h], [y, y], color=color, lw=lw, zorder=1)


def draw_body(ax, body, **kw):
    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(body.xc + body.R * np.cos(t), body.yc + body.R * np.sin(t),
            **({"color": C_BODY, "lw": 1.3, "zorder": 5} | kw))


# ---------------------------------------------------------------------------
@figure
def fig_geometry():
    """(a) cell classes around a cylinder, (b) one cut cell: crossings, wet
    faces, sub-samples, (c) the closure polygon."""
    body = cg.Circle(0.0, 0.0, 1.0)
    g = cg.Grid(-0.17, -0.13, 0.2, 8, 8)
    thx, thy, alpha, cen, absd = cg.fractions(g, body, ns=32)

    fig, axs = plt.subplots(1, 3, figsize=(6.6, 2.35),
                            gridspec_kw={"width_ratios": [1, 1, 0.9]})
    ax = axs[0]
    for i in range(g.nx):
        for j in range(g.ny):
            c = g.corner(i, j)
            a = alpha[i, j]
            col = "white" if a == 1 else (C_SOLID if a == 0 else C_CUT)
            if absd[i, j]:
                col = "#f4c7a1"
            ax.add_patch(Rectangle(c, g.h, g.h, fc=col, ec="none", zorder=0))
    draw_grid(ax, g)
    draw_body(ax, body)
    ax.text(0.12, 0.12, "body\n$\\varphi<0$", ha="center", va="center", fontsize=8)
    ax.text(1.2, 1.25, "fluid\n$\\varphi>0$", ha="center", va="center", fontsize=8)
    handles = [Rectangle((0, 0), 1, 1, fc="white", ec="0.5", lw=0.5),
               Rectangle((0, 0), 1, 1, fc=C_CUT, ec="none"),
               Rectangle((0, 0), 1, 1, fc=C_SOLID, ec="none")]
    ax.legend(handles, ["fluid $\\alpha=1$", "cut $0<\\alpha<1$", "solid $\\alpha=0$"],
              loc="upper left", bbox_to_anchor=(-0.02, -0.02), ncol=2, frameon=False,
              handlelength=1.0, columnspacing=0.8, fontsize=7)
    ax.set_title("(a)", loc="left")
    ax.set_xlim(g.x0, g.x0 + g.nx * g.h)
    ax.set_ylim(g.y0, g.y0 + g.ny * g.h)

    # (b) one cut cell
    ax = axs[1]
    i, j = 4, 4
    gz = cg.Grid(*g.corner(i - 1, j - 1), g.h, 3, 3)
    draw_grid(ax, gz, lw=0.4)
    c0 = g.corner(i, j)
    h = g.h
    ax.add_patch(Rectangle(c0, h, h, fc="none", ec="k", lw=0.9, zorder=3))
    poly = cell_polygon(g, body, i, j)
    ax.add_patch(Polygon(poly, fc=C_CUT, ec="none", zorder=0))
    t = np.linspace(0, np.pi / 2, 200)
    ax.plot(np.cos(t), np.sin(t), color=C_BODY, lw=1.2, zorder=5)
    ax.add_patch(Polygon(np.r_[np.c_[np.cos(t), np.sin(t)], [[0, 0]]], fc=C_SOLID,
                         ec="none", zorder=-1))
    ns = 8
    s = (np.arange(ns) + 0.5) / ns
    for a in s:
        for b in s:
            x, y = c0[0] + a * h, c0[1] + b * h
            wet = body(x, y) > 0
            ax.plot(x, y, "o", ms=1.6, color=C_ARROW if wet else "0.5", zorder=4)
    corners = [g.corner(i, j), g.corner(i + 1, j), g.corner(i + 1, j + 1), g.corner(i, j + 1)]
    for p in corners:
        wet = body(*p) > 0
        ax.plot(*p, "s", ms=3.5, mfc="white" if wet else "k", mec="k", zorder=6)
    for k in range(4):
        a, b = corners[k], corners[(k + 1) % 4]
        fr, seg = cg.wet_segment(body, a, b)
        if fr > 0:
            seg = np.array(seg)
            ax.plot(seg[:, 0], seg[:, 1], color=C_WET, lw=3, zorder=2, solid_capstyle="butt")
        if 0 < fr < 1:
            xa, xb = (a, b) if body(*a) > 0 else (b, a)
            lo, hi = np.array(xa, float), np.array(xb, float)
            for it in range(5):
                m = 0.5 * (lo + hi)
                ax.plot(*m, "|" if a[1] == b[1] else "_", ms=5 - 0.6 * it, color=C_WALL,
                        mew=0.8, zorder=6)
                if body(*m) > 0:
                    lo = m
                else:
                    hi = m
            x = cg.bisect(body, xa, xb)
            ax.plot(*x, "o", ms=4, mfc=C_WALL, mec="k", mew=0.5, zorder=7)
    crosses = [p for p in poly if body(*p) ** 2 < 1e-20]
    crosses = np.array(crosses)
    ax.plot(crosses[:, 0], crosses[:, 1], "--", color=C_WALL, lw=0.9, zorder=6)
    ax.text(c0[0] + 0.78 * h, c0[1] - 0.03, "$\\theta_S|\\mathbf{S}_S|$", ha="center",
            va="top", color=C_WET, fontsize=8)
    ax.text(c0[0] + h + 0.01, c0[1] + 0.55 * h, "$\\theta_E|\\mathbf{S}_E|$", ha="left",
            va="center", color=C_WET, fontsize=8)
    ax.set_xlim(gz.x0, gz.x0 + 3 * h)
    ax.set_ylim(gz.y0, gz.y0 + 3 * h)
    ax.set_title("(b)", loc="left")

    # (c) closure polygon, cell (i, j)
    ax = axs[2]
    vec = [("$\\theta_E\\mathbf{S}_E$", np.array([thx[i + 1, j] * h, 0.0])),
           ("$\\theta_N\\mathbf{S}_N$", np.array([0.0, thy[i, j + 1] * h])),
           ("$\\theta_W\\mathbf{S}_W$", np.array([-thx[i, j] * h, 0.0])),
           ("$\\theta_S\\mathbf{S}_S$", np.array([0.0, -thy[i, j] * h]))]
    p = np.zeros(2)
    for lab, v in vec:
        if np.hypot(*v) < 1e-12:
            continue
        arrow(ax, p, p + v, color=C_WET, lw=1.6)
        m = p + 0.5 * v
        off = np.array([v[1], -v[0]]) / np.hypot(*v) * 0.022
        ax.text(*(m + off), lab, color=C_WET, ha="center", va="center", fontsize=8)
        p = p + v
    Sw = cg.wall_vector(g, thx, thy, i, j)
    arrow(ax, p, p + Sw, color=C_WALL, lw=1.6)
    m = p + 0.5 * Sw
    ax.text(m[0] - 0.02, m[1] + 0.01, "$\\mathbf{S}_w=\\mathbf{n}_wA_{\\rm wall}$",
            color=C_WALL, ha="right", va="bottom", fontsize=8)
    ax.plot(0, 0, "ko", ms=2.5)
    ax.set_aspect("equal")
    ax.set_xlim(-0.07, 0.27)
    ax.set_ylim(-0.09, 0.25)
    ax.set_title("(c)", loc="left")
    ax.text(0.10, -0.07, "$\\sum_f\\theta_f\\mathbf{S}_f+\\mathbf{S}_w=\\mathbf{0}$",
            ha="center", fontsize=8)
    for a in axs:
        a.set_aspect("equal")
        a.set_xticks([])
        a.set_yticks([])
        for sp in a.spines.values():
            sp.set_visible(a is not axs[2])
    fig.savefig("fig_geometry.pdf")


# ---------------------------------------------------------------------------
def draw_cells(ax, g, alpha, absd=None, ij=None):
    for i in range(g.nx):
        for j in range(g.ny):
            if ij is not None and (i, j) not in ij:
                continue
            a = alpha[i, j]
            col = "white" if a == 1 else (C_SOLID if a == 0 else C_CUT)
            if absd is not None and absd[i, j]:
                col = "#f4c7a1"
            ax.add_patch(Rectangle(g.corner(i, j), g.h, g.h, fc=col, ec="none", zorder=0))


def draw_faces(ax, g, thx, thy, lw=2.4):
    """Open part of every face in green, closed faces (theta = 0) in black."""
    h = g.h
    for i in range(g.nx + 1):
        for j in range(g.ny):
            c = g.corner(i, j)
            t = thx[i, j]
            if t == 0:
                ax.plot([c[0], c[0]], [c[1], c[1] + h], color="k", lw=0.9, zorder=2)
    for i in range(g.nx):
        for j in range(g.ny + 1):
            c = g.corner(i, j)
            if thy[i, j] == 0:
                ax.plot([c[0], c[0] + h], [c[1], c[1]], color="k", lw=0.9, zorder=2)


@figure
def fig_sliver():
    """Sliver absorption: before and after, the wall vectors of the
    neighbours, and the order of absorption and coupled-face sync."""
    body = cg.Circle(0.0, 0.0, 1.0)
    g = cg.Grid(-0.19, -0.17, 0.2, 8, 8)
    b0 = cg.fractions(g, body, ns=32, absorb=False)
    b1 = cg.fractions(g, body, ns=32, absorb=True)
    isl = tuple(np.argwhere(b1[4])[0])
    fig, axs = plt.subplots(1, 3, figsize=(6.6, 2.4),
                            gridspec_kw={"width_ratios": [1, 1, 1.05]})
    t = np.linspace(0, np.pi / 2, 300)
    zoom = [(isl[0] + a, isl[1] + b) for a in (-1, 0, 1) for b in (-1, 0, 1)]
    for k, (ax, (thx, thy, alpha, cen, absd)) in enumerate(zip(axs[:2], (b0, b1))):
        draw_cells(ax, g, alpha, absd if k == 1 else None)
        draw_grid(ax, g, lw=0.3)
        draw_faces(ax, g, thx, thy)
        ax.plot(np.cos(t), np.sin(t), color=C_BODY, lw=1.2, zorder=5)
        for (i, j) in zoom:
            if not (0 <= i < g.nx and 0 <= j < g.ny):
                continue
            if 0 < alpha[i, j] < 1 or (alpha[i, j] == 1 and np.hypot(*cg.wall_vector(g, thx, thy, i, j)) > 1e-12):
                Sw = cg.wall_vector(g, thx, thy, i, j)
                p = cen[i, j] if alpha[i, j] < 1 else g.centre(i, j)
                arrow(ax, p, p + 0.8 * Sw, color=C_WALL, lw=1.0, ms=6, zorder=8)
        lo = g.corner(isl[0] - 1, isl[1] - 1)
        ax.set_xlim(lo[0], lo[0] + 3 * g.h)
        ax.set_ylim(lo[1], lo[1] + 3 * g.h)
        a = b0[2][isl]
        ax.set_title("(a) before, $\\alpha_P=%.3f$" % a if k == 0 else
                     "(b) after absorption", loc="left")
    axs[0].annotate("sliver", xy=g.corner(*isl) + np.array([0.19, 0.19]),
                    xytext=g.corner(*isl) + np.array([0.02, 0.05]), fontsize=8,
                    arrowprops=dict(arrowstyle="-", lw=0.6))
    axs[1].text(*(g.corner(*isl) + 0.5 * g.h), "absorbed\n$\\alpha:=0$\n$\\theta_f:=0$",
                ha="center", va="center", fontsize=7)

    # (c) periodic seam: two orders of the operations
    ax = axs[2]
    h = 1.0
    for row, (y0, good) in enumerate(((2.2, False), (0.0, True))):
        for x0 in (0.0, 2.4):
            ax.add_patch(Rectangle((x0, y0), h, h, fc="white", ec="0.5", lw=0.6))
        ax.add_patch(Rectangle((0.0, y0), 0.8, h, fc=C_SOLID, ec="none"))
        ax.add_patch(Rectangle((0.8, y0), 0.2, h, fc="#f4c7a1", ec="none"))
        ax.plot([1, 1], [y0, y0 + h], color="k", lw=2.2)
        ax.plot([2.4, 2.4], [y0, y0 + h], color="k" if good else C_WET, lw=2.2)
        ax.annotate("", xy=(2.35, y0 + 0.5), xytext=(1.05, y0 + 0.5),
                    arrowprops=dict(arrowstyle="<->", lw=0.6, color=C_GREY, ls="--"))
        ax.text(1.7, y0 + 0.58, "coupled", ha="center", fontsize=6.5, color=C_GREY)
        ax.text(0.5, y0 - 0.06, "$x=L$", ha="center", va="top", fontsize=7)
        ax.text(2.9, y0 - 0.06, "$x=0$", ha="center", va="top", fontsize=7)
        ax.text(3.45, y0 + 0.5, "ok" if good else "wrong", fontsize=7,
                color=C_WET if good else C_WALL, va="center")
    ax.text(0.0, 3.3, "sync, then absorb: open at $x=0$", fontsize=7, va="bottom")
    ax.text(0.0, 1.1, "absorb, then $\\theta_f:=\\min(\\theta_f,\\theta_{f'})$", fontsize=7,
            va="bottom")
    ax.set_xlim(-0.1, 3.9)
    ax.set_ylim(-0.4, 3.7)
    ax.set_title("(c) periodic or processor seam", loc="left")
    for a in axs:
        a.set_aspect("equal")
        a.set_xticks([])
        a.set_yticks([])
    for sp in axs[2].spines.values():
        sp.set_visible(False)
    fig.savefig("fig_sliver.pdf")


# ---------------------------------------------------------------------------
@figure
def fig_walldist():
    """(a) the two wall-distance estimates in a cut cell; (b) their ratio over
    all cut cells of a cylinder for many grid offsets."""
    body = cg.Circle(0.0, 0.0, 1.0)
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.6), gridspec_kw={"width_ratios": [1, 1.35]})
    ax = axs[0]
    h = 1.0
    ax.add_patch(Rectangle((0, 0), h, h, fc="none", ec="k", lw=0.9))
    wall = np.array([[0.0, 0.55], [1.0, 0.05]])
    poly = np.array([[0, 0.55], [1.0, 0.05], [1, 1], [0, 1]])
    ax.add_patch(Polygon(poly, fc=C_CUT, ec="none", zorder=0))
    ax.add_patch(Polygon([[0, 0], [1, 0], [1, 0.05], [0, 0.55]], fc=C_SOLID, ec="none"))
    ax.plot(wall[:, 0], wall[:, 1], color=C_WALL, lw=2)
    # centroid of polygon
    x, y = poly[:, 0], poly[:, 1]
    xs, ys = np.roll(x, -1), np.roll(y, -1)
    cr = x * ys - xs * y
    A = cr.sum() / 2
    cx, cy = ((x + xs) * cr).sum() / (6 * A), ((y + ys) * cr).sum() / (6 * A)
    t = wall[1] - wall[0]
    Aw = np.hypot(*t)
    n = np.array([t[1], -t[0]]) / Aw          # into the body (downwards)
    dc = np.dot(np.array([cx, cy]) - wall[0], -n)
    foot = np.array([cx, cy]) + dc * n
    ax.plot(cx, cy, "o", color=C_ARROW, ms=4, zorder=6)
    ax.text(cx + 0.03, cy + 0.04, "$\\mathbf{x}_c$", color=C_ARROW, fontsize=9)
    ax.plot(0.5, 0.5, "+", color="k", ms=6)
    ax.text(0.53, 0.46, "$\\mathbf{x}_P$", fontsize=8)
    ax.annotate("", xy=foot, xytext=(cx, cy), arrowprops=dict(arrowstyle="<->", lw=0.9, color=C_ARROW))
    ax.text(*(0.5 * (foot + [cx, cy]) + [0.03, -0.02]), "$d_c$", color=C_ARROW, fontsize=9)
    # slab of equal volume
    slab = A / Aw
    p0 = wall[0] + 0.80 * t
    q0 = p0 - n * slab
    ax.annotate("", xy=q0, xytext=p0, arrowprops=dict(arrowstyle="<->", lw=0.8, color=C_GREY))
    ax.text(*(0.5 * (p0 + q0) + [0.03, 0.0]), "$\\alpha V/A_{\\rm wall}$", color=C_GREY,
            fontsize=8, ha="left", va="center")
    arrow(ax, 0.5 * (wall[0] + wall[1]), 0.5 * (wall[0] + wall[1]) + 0.22 * n, color=C_WALL)
    ax.text(*(0.5 * (wall[0] + wall[1]) + 0.24 * n + [0.03, 0]), "$\\mathbf{n}_w$", color=C_WALL,
            fontsize=9)
    ax.set_xlim(-0.05, 1.45)
    ax.set_ylim(-0.05, 1.05)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("(a)", loc="left")
    ax.text(0.5, -0.02, "$d_{\\rm wall}=\\alpha V/(2A_{\\rm wall})$ (first order)\n"
            "$d_c=$ distance of $\\mathbf{x}_c$ to the wall (second order)",
            ha="center", va="top", fontsize=7.5)

    ax = axs[1]
    rng = np.random.default_rng(1)
    A_, R_ = [], []
    for k in range(40):
        ox, oy = rng.uniform(-0.1, 0.0, 2)
        g = cg.Grid(ox - 1.2, oy - 1.2, 0.1, 25, 25)
        thx, thy, alpha, cen, absd = cg.fractions(g, body, ns=32, absorb=True)
        for i, j in zip(*np.nonzero((alpha > 0) & (alpha < 1))):
            Sw = cg.wall_vector(g, thx, thy, i, j)
            Aw = np.hypot(*Sw)
            if Aw < 1e-12:
                continue
            dwall = alpha[i, j] * g.h ** 2 / (2 * Aw)
            dc = body(*cen[i, j])
            A_.append(alpha[i, j])
            R_.append(dwall / dc)
    A_, R_ = np.array(A_), np.array(R_)
    ax.scatter(A_, R_, s=2, color=C_ARROW, alpha=0.5, lw=0)
    ax.axhline(1, color="k", lw=0.6)
    ax.set_xlabel("cell fluid fraction $\\alpha_P$")
    ax.set_ylabel("$d_{\\rm wall}/d_c$")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 3)
    ax.set_title("(b) cylinder $R=1$, $h=0.1$, 40 grid offsets", loc="left")
    fig.tight_layout()
    fig.savefig("fig_walldist.pdf")


# ---------------------------------------------------------------------------
@figure
def fig_swept():
    """A cut cell over one time step: fixed faces, moving wall, swept band."""
    fig, ax = plt.subplots(figsize=(3.4, 2.7))
    h = 1.0
    ub = np.array([0.18, 0.36])
    w0 = np.array([[0.0, 0.40], [1.0, 0.05]])
    w1 = w0 + ub

    def clip_poly(w):
        # fluid above the line through w, inside the unit square
        pts = []
        sq = [[0, 0], [1, 0], [1, 1], [0, 1]]
        t = w[1] - w[0]
        side = lambda p: t[0] * (p[1] - w[0][1]) - t[1] * (p[0] - w[0][0])
        for k in range(4):
            a, b = np.array(sq[k], float), np.array(sq[(k + 1) % 4], float)
            sa, sb = side(a), side(b)
            if sa > 0:
                pts.append(a)
            if (sa > 0) != (sb > 0):
                pts.append(a + (b - a) * sa / (sa - sb))
        return np.array(pts)

    P0, P1 = clip_poly(w0), clip_poly(w1)
    ax.add_patch(Polygon(P1, fc=C_CUT, ec="none", zorder=0))
    band = np.array([w0[0] if w0[0][0] >= 0 else w0[0], ])
    # swept band = P0 minus P1
    from matplotlib.path import Path
    ax.add_patch(Polygon(P0, fc="#f6d5a8", ec="none", zorder=-1))
    ax.add_patch(Polygon(np.array([[0, 0], [1, 0], [1, 1], [0, 1]]), fc="none", ec="k", lw=1.0, zorder=3))
    for w, col, lab in ((w0, C_WALL, "$t^n$"), (w1, C_ARROW, "$t^{n+1}$")):
        t = w[1] - w[0]
        s0 = max(0, (0 - w[0][0]) / t[0])
        ax.plot(*np.c_[w[0] - 0.15 * t, w[1] + 0.15 * t], color=col, lw=1.6, zorder=4)
        ax.text(*(w[1] + 0.18 * t + [0.02, -0.02]), "wall at " + lab, color=col, fontsize=8,
                va="center")
    m = 0.5 * (w0[0] + w0[1]) + np.array([-0.1, 0.05])
    arrow(ax, m, m + ub * 0.95, color="k", lw=1.2)
    ax.text(*(m + 0.5 * ub + [-0.12, 0.02]), "$\\mathbf{u}_b\\Delta t$", fontsize=9)
    t = w1[1] - w1[0]
    n = np.array([t[1], -t[0]]) / np.hypot(*t)
    c = w1[0] + 0.18 * (w1[1] - w1[0])
    arrow(ax, c, c + 0.25 * n, color=C_ARROW, lw=1.2)
    ax.text(*(c + 0.27 * n + [-0.1, 0.0]), "$\\mathbf{n}_w$", color=C_ARROW, fontsize=9)
    ax.text(0.62, 0.78, "$\\alpha_P^{n+1}V_P$", fontsize=9, ha="center")
    ax.text(0.84, 0.21, "swept: $\\Delta t\\,\\mathbf{u}_b\\cdot\\mathbf{S}_w$", fontsize=8,
            ha="center", rotation=-19)
    ax.text(0.5, -0.07, "faces fixed: only $\\theta_f$ changes", ha="center", va="top", fontsize=8)
    ax.set_xlim(-0.25, 1.55)
    ax.set_ylim(-0.2, 1.1)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.savefig("fig_swept.pdf")


# ---------------------------------------------------------------------------
@figure
def fig_gcl():
    """Geometry-first vs flux-first mass source for a translating cylinder."""
    body0 = cg.Circle(0.0, 0.0, 1.0)
    h = 0.1
    g = cg.Grid(-1.5, -1.5, h, 30, 30)
    ub = np.array([0.5, 0.0])
    dt = 0.02
    V = h * h
    nsteps = 25
    prev = cg.fractions(g, body0, ns=32)
    sg, sf, err = [], [], []
    for n in range(1, nsteps + 1):
        body = body0.moved(*(ub * n * dt))
        cur = cg.fractions(g, body, ns=32)
        thx, thy, alpha = cur[0], cur[1], cur[2]
        Sgeo = (alpha - prev[2]) * V / dt
        Sflx = np.zeros_like(alpha)
        for i in range(g.nx):
            for j in range(g.ny):
                if alpha[i, j] > 0:
                    Sflx[i, j] = np.dot(ub, cg.wall_vector(g, thx, thy, i, j))
        # dV/dt = -S  (a body advancing into the fluid removes volume)
        sg.append(abs(Sgeo.sum()))
        sf.append(abs(Sflx.sum()))
        m = (alpha > 0) | (prev[2] > 0)
        cutm = m & ((alpha < 1) | (prev[2] < 1))
        err.append(np.abs(Sgeo[cutm] + Sflx[cutm]) / (np.abs(Sflx[cutm]).max()))
        prev = cur
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.4))
    ax = axs[0]
    k = np.arange(1, nsteps + 1)
    ax.semilogy(k, np.maximum(sg, 1e-18), "o-", ms=3, color=C_WALL,
                label="geometry first, $(\\alpha^{n+1}-\\alpha^n)V/\\Delta t$")
    ax.semilogy(k, np.maximum(sf, 1e-18), "s-", ms=3, color=C_ARROW,
                label="flux first, $\\mathbf{u}_b\\cdot\\mathbf{S}_w$")
    ax.set_xlabel("time step $n$")
    ax.set_ylabel("$|\\sum_P S_P|$")
    ax.set_ylim(1e-19, 1e-1)
    ax.legend(loc="center right", frameon=False, fontsize=7)
    ax.set_title("(a) compatibility of the source", loc="left")
    ax = axs[1]
    e = np.concatenate(err)
    ax.hist(np.log10(np.maximum(e, 1e-17)), bins=40, color=C_WALL, alpha=0.8)
    ax.set_xlabel("$\\log_{10}\\,|S^{\\rm geo}_P-S^{\\rm flux}_P|\\,/\\max|S^{\\rm flux}|$")
    ax.set_ylabel("cut cells")
    ax.set_title("(b) per-cell discrepancy", loc="left")
    fig.tight_layout()
    fig.savefig("fig_gcl.pdf")
    print("  geometry-first sum: median %.2e, flux-first max %.2e" % (np.median(sg), max(sf)))


# ---------------------------------------------------------------------------
@figure
def fig_centroids():
    """Second-order ingredients: fluid centroids, centroid vector d, its
    normal and tangential parts, wet-face centroid, centroid wall distance."""
    body = cg.Circle(0.0, 0.0, 1.0)
    g = cg.Grid(-0.17, -0.13, 0.2, 8, 8)
    thx, thy, alpha, cen, absd = cg.fractions(g, body, ns=64)
    P, N = (5, 2), (5, 3)
    fig, ax = plt.subplots(figsize=(3.4, 3.8))
    lo = g.corner(4, 1)
    gz = cg.Grid(*lo, g.h, 3, 4)
    draw_grid(ax, gz, lw=0.4)
    for c in (P, N):
        ax.add_patch(Polygon(cell_polygon(g, body, *c), fc=C_CUT, ec="none", zorder=0))
    t = np.linspace(0, np.pi / 2, 300)
    ax.add_patch(Polygon(np.r_[np.c_[np.cos(t), np.sin(t)], [[0, 0]]], fc=C_SOLID, ec="none", zorder=-1))
    ax.plot(np.cos(t), np.sin(t), color=C_BODY, lw=1.2, zorder=5)
    f0, f1 = g.corner(5, 3), g.corner(6, 3)      # y-face between P (below) and N (above)
    fr, seg = cg.wet_segment(body, f0, f1)
    seg = np.array(seg)
    ax.plot(seg[:, 0], seg[:, 1], color=C_WET, lw=3, zorder=3)
    xf = seg.mean(axis=0)
    ax.plot(*xf, "D", ms=3.5, color=C_WET, zorder=7)
    ax.text(xf[0] + 0.012, xf[1] + 0.004, "$\\mathbf{x}_f'$", color=C_WET, fontsize=9)
    ax.text(seg[-1][0] + 0.004, f0[1] - 0.004, "wet face $f$", color=C_WET, fontsize=7, va="top")
    for c, lab in ((P, "P"), (N, "N")):
        xp = g.centre(*c)
        ax.plot(*xp, "+", color="k", ms=7, zorder=6)
        ax.text(xp[0] + 0.006, xp[1] + 0.004, "$\\mathbf{x}_%s$" % lab, fontsize=9, ha="left",
                va="bottom")
        ax.plot(*cen[c], "o", color=C_ARROW, ms=4, zorder=7)
        ax.text(cen[c][0] + 0.01, cen[c][1] - 0.004, "$\\mathbf{c}_%s$" % lab, color=C_ARROW,
                fontsize=9, ha="left", va="top")
    cP, cN = cen[P], cen[N]
    d = cN - cP
    arrow(ax, cP, cN, color=C_ARROW, lw=1.3)
    ax.text(*(cP + 0.3 * d + [0.01, 0.0]), "$\\mathbf{d}$", color=C_ARROW, fontsize=10)
    nvec = np.array([0.0, 1.0])
    corner = cP + np.dot(d, nvec) * nvec
    ax.plot(*np.c_[cP, corner], ":", color=C_GREY, lw=1)
    ax.plot(*np.c_[corner, cN], ":", color=C_GREY, lw=1)
    ax.text(corner[0] - 0.004, cP[1] + 0.3 * d[1], "$(\\mathbf{n}\\cdot\\mathbf{d})\\mathbf{n}$",
            color=C_GREY, fontsize=8, ha="right", va="center",
            bbox=dict(fc="white", ec="none", pad=0.5, alpha=0.8))
    ax.text(0.5 * (corner[0] + cN[0]), corner[1] + 0.006, "$\\mathbf{k}$", color=C_GREY,
            fontsize=9, ha="center", va="bottom")
    arrow(ax, f0 + [0.03, 0], f0 + [0.03, 0.06], color="k", lw=1.0, ms=6)
    ax.text(f0[0] + 0.037, f0[1] + 0.05, "$\\mathbf{n}$", fontsize=9)
    r = np.hypot(*cP)
    foot = cP / r
    ax.annotate("", xy=foot, xytext=cP, arrowprops=dict(arrowstyle="<->", lw=0.8, color=C_WALL))
    ax.text(*(0.5 * (foot + cP) + [-0.004, -0.02]), "$d_c$", color=C_WALL, fontsize=9)
    c52 = g.corner(5, 2)
    ax.set_xlim(c52[0] - 0.5 * g.h, c52[0] + 1.6 * g.h)
    ax.set_ylim(c52[1] - 0.1 * g.h, c52[1] + 2.1 * g.h)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.savefig("fig_centroids.pdf")


# ---------------------------------------------------------------------------
ORDER_N = np.array([20, 40, 80, 160, 320])
ORDER = {  # L2 errors, doc/secondOrderStudy.md (Python study; the OpenFOAM operators
           # reproduce them to the printed digits for N <= 160)
    "LD": {"V0/F0": [1.79e-2, 1.14e-2, 4.92e-3, 2.40e-3, 1.23e-3],
           "V1/F0": [1.65e-2, 9.45e-3, 4.09e-3, 2.14e-3, np.nan],
           "V1/F1": [1.55e-3, 4.09e-4, 9.81e-5, 2.31e-5, 5.52e-6],
           "V3/F1": [2.32e-3, 6.00e-4, 1.52e-4, 3.48e-5, np.nan]},
    "CD": {"V0/F0": [1.13e-2, 7.71e-3, 3.76e-3, 1.62e-3, 9.45e-4],
           "V1/F1": [3.46e-3, 1.06e-3, 2.78e-4, 6.97e-5, 1.77e-5],
           "V3/F1": [5.45e-3, 1.71e-3, 4.55e-4, 1.12e-4, np.nan]},
    "LN": {"V0/F0": [4.09e-3, 8.45e-4, 1.56e-4, 2.01e-5, 5.57e-6],
           "V1/F1": [3.14e-3, 7.41e-4, 1.48e-4, 2.66e-5, 6.44e-6]},
    "CN": {"V0/F0": [3.28e-3, 8.53e-4, 2.29e-4, 6.45e-5, 1.75e-5],
           "V1/F1": [2.95e-3, 8.16e-4, 2.22e-4, 6.36e-5, 1.73e-5]},
}
TITLES = {"LD": "Laplacian, Dirichlet wall", "CD": "convection--diffusion, Dirichlet",
          "LN": "Laplacian, Neumann wall", "CN": "convection--diffusion, Neumann"}
STYLE = {"V0/F0": ("o-", C_WALL, "first order (V0/F0)"),
         "V1/F0": ("v--", "#e7a33e", "centroid $d_c$ only (V1/F0)"),
         "V1/F1": ("s-", C_ARROW, "second order (V1/F1)"),
         "V3/F1": ("^:", C_GREY, "V1/F1 + normal probe (V3/F1)")}


def slope(ax, x0, y0, p, L=2.0, label=True):
    x = np.array([x0, x0 * L])
    ax.plot(x, y0 * (x / x0) ** (-p), color="k", lw=0.6)
    if label:
        ax.text(x0 * L ** 0.5, y0 * L ** (-p / 2) * 0.72, "%d" % p, fontsize=7, ha="center",
                va="top")


OF_N = np.array([20, 40, 80, 160])
OF = {  # OpenFOAM ibmOrderTest, L2 errors (tutorials/order/laplaceDisk), switch off / on
    "LD": ([1.7914e-2, 1.1436e-2, 4.9168e-3, 2.3975e-3], [1.5474e-3, 4.0919e-4, 9.8092e-5, 2.3094e-5]),
    "LN": ([4.0891e-3, 8.4503e-4, 1.5630e-4, 2.0065e-5], [3.1415e-3, 7.4060e-4, 1.4751e-4, 2.6643e-5]),
    "CD": ([1.1300e-2, 7.7073e-3, 3.7591e-3, 1.6219e-3], [4.0577e-3, 1.1950e-3, 3.1110e-4, 7.8439e-5]),
    "CN": ([3.2762e-3, 8.5301e-4, 2.2855e-4, 6.4525e-5], [3.3315e-3, 8.8509e-4, 2.3652e-4, 6.6669e-5]),
}
OF_LINF = {
    "LD": ([2.2075e-2, 1.3052e-2, 6.6080e-3, 3.3128e-3], [1.7212e-3, 4.5010e-4, 1.1392e-4, 2.6239e-5]),
    "LN": ([1.6125e-2, 4.2247e-3, 1.1307e-3, 2.5174e-4], [1.2879e-2, 4.7041e-3, 1.2893e-3, 3.2518e-4]),
    "CD": ([2.4890e-2, 2.1532e-2, 1.2231e-2, 4.9777e-3], [8.4355e-3, 2.8732e-3, 7.8843e-4, 2.1025e-4]),
    "CN": ([8.5684e-3, 2.4747e-3, 7.0345e-4, 2.0234e-4], [1.1785e-2, 3.5292e-3, 1.0289e-3, 2.8711e-4]),
}


@figure
def fig_order_operators():
    fig, axs = plt.subplots(1, 4, figsize=(6.8, 2.25))
    for ax, key in zip(axs, ("LD", "CD", "LN", "CN")):
        off, on = OF[key]
        loff, lon = OF_LINF[key]
        ax.loglog(OF_N, off, "o-", color=C_WALL, ms=3, lw=1, label="first order, $L_2$")
        ax.loglog(OF_N, on, "s-", color=C_ARROW, ms=3, lw=1, label="second order, $L_2$")
        ax.loglog(OF_N, loff, "o:", color=C_WALL, ms=2.5, lw=0.8, mfc="white",
                  label="first order, $L_\\infty$")
        ax.loglog(OF_N, lon, "s:", color=C_ARROW, ms=2.5, lw=0.8, mfc="white",
                  label="second order, $L_\\infty$")
        ax.set_title(TITLES[key], fontsize=7.5)
        ax.set_xlabel("$N$")
        slope(ax, 45, on[1] * 0.3, 2)
        slope(ax, 45, max(off[1], loff[1]) * 2.0, 1)
        ax.set_xticks(list(OF_N))
        ax.set_xticklabels([str(n) for n in OF_N], fontsize=6.5)
        ax.minorticks_off()
    axs[0].set_ylabel("error")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.08),
               fontsize=7)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig("fig_order_operators.pdf")


@figure
def fig_order_ablation():
    fig, axs = plt.subplots(1, 2, figsize=(6.0, 2.5))
    for ax, key in zip(axs, ("LD", "CD")):
        for var, e in ORDER[key].items():
            st, col, lab = STYLE[var]
            if key == "CD" and var == "V1/F1":
                lab = "second order + C1 (V1/C1/F1)"
            ax.loglog(ORDER_N, e, st, color=col, ms=3, lw=1, label=lab)
        if key == "LD":
            ax.loglog(ORDER_N, [2.49e-3, 6.58e-4, 1.63e-4, 5.18e-5, np.nan], "D-.", color="#7b3294",
                      ms=3, lw=1, label="V1/F1 + row scaling (V2/F1)")
        ax.set_title(TITLES[key], fontsize=8)
        ax.set_xlabel("$N$")
        slope(ax, 60, ORDER[key]["V1/F1"][1] * 0.35, 2)
        slope(ax, 60, ORDER[key]["V0/F0"][1] * 1.7, 1)
        ax.set_xticks(list(ORDER_N))
        ax.set_xticklabels([str(n) for n in ORDER_N], fontsize=7)
        ax.minorticks_off()
        ax.legend(frameon=False, fontsize=6.3, loc="lower left")
    axs[0].set_ylabel("$L_2$ error")
    fig.tight_layout()
    fig.savefig("fig_order_ablation.pdf")


@figure
def fig_order_coupled():
    """Taylor--Couette errors and cylinder drag with its components."""
    fig, axs = plt.subplots(1, 3, figsize=(6.8, 2.4))
    ax = axs[0]
    N = np.array([40, 80, 160])
    ax.loglog(N, [6.479e-3, 4.478e-3, 2.104e-3], "o-", color=C_WALL, ms=3, label="$L_2(u)$, first")
    ax.loglog(N, [1.159e-3, 3.353e-4, 1.274e-4], "s-", color=C_ARROW, ms=3, label="$L_2(u)$, second")
    ax.loglog(N, [3.633e-2, 3.337e-2, 2.295e-2], "o--", color=C_WALL, ms=3, mfc="white",
              label="$L_2(p)$, first")
    ax.loglog(N, [2.169e-2, 1.697e-2, 1.296e-2], "s--", color=C_ARROW, ms=3, mfc="white",
              label="$L_2(p)$, second")
    slope(ax, 60, 1.6e-4, 2, L=1.6)
    slope(ax, 60, 8.5e-3, 1, L=1.6)
    ax.set_xlabel("$N$")
    ax.set_ylabel("error")
    ax.set_xticks([40, 80, 160])
    ax.set_xticklabels(["40", "80", "160"])
    ax.minorticks_off()
    ax.legend(frameon=False, fontsize=6, loc="lower left", ncol=1)
    ax.set_title("(a) Taylor--Couette", loc="left")
    ax = axs[1]
    N = np.array([41, 81, 161, 321])
    comp = {"first": ([2.7965, 2.8689, 2.8928, 2.9074], [2.7465, 2.7915, 2.8649, 2.8878]),
            "second": ([2.8017, 2.8599, 2.8903, 2.9065], [2.9817, 2.9606, 2.9418, 2.9288])}
    for tr, (fp, fv) in comp.items():
        col = C_WALL if tr == "first" else C_ARROW
        ax.plot(1 / N, fp, "-", color=col, marker="^", ms=3, label="pressure, %s" % tr)
        ax.plot(1 / N, fv, "--", color=col, marker="v", ms=3, label="viscous, %s" % tr)
    ax.set_xlabel("$1/N$")
    ax.set_ylabel("drag component per unit depth")
    ax.set_xlim(0, 0.026)
    ax.legend(frameon=False, fontsize=5.8, loc="lower left")
    ax.set_title("(b) no slip, components", loc="left")
    ax = axs[2]
    tot = {("no slip", "first"): np.add(*comp["first"]),
           ("no slip", "second"): np.add(*comp["second"]),
           ("free slip", "first"): np.array([2.837673, 2.942117, 3.000386, 3.035406]),
           ("free slip", "second"): np.array([2.916107, 2.984759, 3.023749, 3.051058])}
    for (w, tr), F in tot.items():
        col = C_WALL if tr == "first" else C_ARROW
        ref = tot[(w, "second")][-1]
        ax.plot(1 / N, F / ref, ("o" if tr == "first" else "s") + ("-" if w == "no slip" else "--"),
                color=col, ms=3, mfc=col if w == "no slip" else "white", label="%s, %s" % (w, tr))
    ax.set_xlabel("$1/N$")
    ax.set_ylabel("$F/F_{321}^{\\rm second}$")
    ax.set_xlim(0, 0.026)
    ax.legend(frameon=False, fontsize=5.8, loc="lower left")
    ax.set_title("(c) total drag", loc="left")
    fig.tight_layout()
    fig.savefig("fig_order_coupled.pdf")


def load(name):
    return np.loadtxt("data/" + name)


@figure
def fig_scalar_budget():
    """Budget terms of the scalar in three moving/static runs (tracked solver logs)."""
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.5))
    ax = axs[0]
    b = load("budget_oscillating_fixedValue.dat")
    t = 0.025 * np.arange(1, len(b) + 1)
    ax.plot(t, b[:, 0], color=C_ARROW, label="content")
    ax.plot(t, b[:, 1], color="k", ls="--", label="injected")
    ax.plot(t, b[:, 2], color=C_WALL, label="through walls")
    ax.set_xlabel("$t$")
    ax.legend(frameon=False, fontsize=7)
    ax.set_title("(a) oscillating cylinder, fixed value", loc="left")
    ax = axs[1]
    cases = [("budget_oscillating_zeroFlux.dat", 0.025, "oscillating, zero flux", C_ARROW),
             ("budget_oscillating_fixedValue.dat", 0.025, "oscillating, fixed value", C_WALL),
             ("budget_galilean_zeroFlux.dat", 0.025, "translating, zero flux", "#1b9e3e")]
    for fn, dt, lab, col in cases:
        b = load(fn)
        t = dt * np.arange(1, len(b) + 1)
        rel = np.maximum(np.abs(b[:, 0]), 1e-30)
        ax.semilogy(t, np.abs(b[:, 3]) / rel, color=col, lw=1, label=lab + ": swept defect")
        ax.semilogy(t, np.maximum(np.abs(b[:, 4]), 1e-18) / rel, color=col, lw=0.6, ls=":",
                    label="budget error")
    ax.set_xlabel("$t$")
    ax.set_ylabel("relative to content")
    ax.set_ylim(1e-15, 1e-1)
    ax.legend(frameon=False, fontsize=5.8, ncol=1, loc="upper right")
    ax.set_title("(b) defect and budget error", loc="left")
    fig.tight_layout()
    fig.savefig("fig_scalar_budget.pdf")


@figure
def fig_moving_diag():
    """Per-step diagnostics of the oscillating cylinder."""
    d = load("moving_oscillating.dat")
    t = d[:, 0]
    fig, axs = plt.subplots(1, 2, figsize=(6.4, 2.3))
    ax = axs[0]
    ax.semilogy(t, np.maximum(np.abs(d[:, 1]), 1e-22), ".", ms=2, color=C_ARROW,
                label="$|\\sum_PS_P|$")
    ax.semilogy(t, np.maximum(d[:, 4], 1e-22), ".", ms=2, color=C_WALL,
                label="continuity error (sum local)")
    ax.set_ylim(1e-21, 1e-12)
    ax.set_xlabel("$t$")
    ax.legend(frameon=False, fontsize=7, loc="upper right", markerscale=3)
    ax.set_title("(a) compatibility and continuity", loc="left")
    ax = axs[1]
    ax.plot(t, d[:, 2], color=C_ARROW, lw=0.8)
    ax.axhline(0, color="k", lw=0.5)
    neg = d[:, 3] > 0
    ax.plot(t[neg], d[neg, 2], "o", ms=2.5, color=C_WALL, label="steps with $\\alpha^n_P<0$")
    ax.set_xlabel("$t$")
    ax.set_ylabel("$\\min_P\\alpha^n_P$")
    ax.legend(frameon=False, fontsize=7)
    ax.set_title("(b) smallest swept fraction", loc="left")
    fig.tight_layout()
    fig.savefig("fig_moving_diag.pdf")


# ---------------------------------------------------------------------------
def main(names):
    for n in names or FIGS:
        FIGS[n]()
        print("wrote", n)


if __name__ == "__main__":
    main(sys.argv[1:])
