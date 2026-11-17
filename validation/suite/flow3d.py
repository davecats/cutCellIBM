"""Flow-field pictures of the 3-D cases for the report: the cells cut by an axis-aligned plane of the
(reconstructed, Cartesian) mesh drawn as their bounding rectangles, coloured by a cell field; solid
cells (ibmAlpha = 0) grey, cut cells (0 < ibmAlpha < 1) outlined, the exact section of the immersed
wall as a line, streamlines of the in-plane velocity.

    flow3d.figures(out, R)   (from latexreport.py; R: result dictionaries by id)"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

import common
import flowfield as ff
from flowfield import _faces, _tokens

sys.path.insert(0, str(common.VALIDATION))
from foamio import read_field, latest_time   # noqa: E402

CUT = "#202020"


class Slice:
    """cells of case whose bounding box contains the plane x[axis] = c (within the window
    (u0, u1, v0, v1) of the two in-plane axes, in increasing order)"""

    def __init__(self, case, axis, c, window, time="latest"):
        self.case = Path(case)
        self.time = latest_time(str(self.case)) if time == "latest" else str(time)
        m = self.case / "constant" / "polyMesh"
        pts = _tokens(m / "points", float)[1:].reshape(-1, 3)
        off, verts = _faces(m / "faces")
        own = _tokens(m / "owner", np.int64)[1:]
        nei = _tokens(m / "neighbour", np.int64)[1:]
        n = max(own.max(), nei.max()) + 1
        lo, hi = np.full((n, 3), np.inf), np.full((n, 3), -np.inf)
        size = np.diff(off)
        for s in np.unique(size):
            f = np.flatnonzero(size == s)
            P = pts[verts[off[f][:, None] + np.arange(s)]]
            fl, fh = P.min(axis=1), P.max(axis=1)
            np.minimum.at(lo, own[f], fl)
            np.maximum.at(hi, own[f], fh)
            k = f < len(nei)
            np.minimum.at(lo, nei[f[k]], fl[k])
            np.maximum.at(hi, nei[f[k]], fh[k])
        self.ij = [i for i in range(3) if i != axis]
        i, j = self.ij
        u0, u1, v0, v1 = window
        sel = ((lo[:, axis] <= c) & (hi[:, axis] > c) & (hi[:, i] >= u0) & (lo[:, i] <= u1) &
               (hi[:, j] >= v0) & (lo[:, j] <= v1))
        self.cells = np.flatnonzero(sel)
        L, H = lo[self.cells][:, [i, j]], hi[self.cells][:, [i, j]]
        self.rect = np.stack([L, np.column_stack([H[:, 0], L[:, 1]]), H, np.column_stack([L[:, 0], H[:, 1]])], axis=1)
        self.C = 0.5*(L + H)
        self.window, self.axis, self.c = window, axis, c
        self.alpha = self.field("ibmAlpha")

    def field(self, name):
        for t in (self.time, "0"):
            f = self.case / t / name
            if f.exists():
                v, _ = read_field(str(f))
                return np.asarray(v)[self.cells]
        return None

    def draw(self, ax, v, cmap=ff.SEQ, norm=None, edges=False):
        a = self.alpha
        fl = a > 0
        style = dict(edgecolors="#ffffff", linewidths=0.25) if edges else dict(linewidths=0, antialiaseds=False)
        pc = PolyCollection(self.rect[fl], array=v[fl], cmap=cmap, norm=norm, rasterized=True, **style)
        ax.add_collection(pc)
        if (~fl).any():
            ax.add_collection(PolyCollection(self.rect[~fl], facecolors=ff.SOLID, rasterized=True, **style))
        cut = fl & (a < 1 - 1e-9)
        if edges and cut.any():
            ax.add_collection(PolyCollection(self.rect[cut], facecolors="none", edgecolors=CUT, linewidths=0.5))
        u0, u1, v0, v1 = self.window
        ax.set_xlim(u0, u1)
        ax.set_ylim(v0, v1)
        ax.set_aspect("equal")
        ax.grid(False)
        return pc

    def streamlines(self, ax, U, n=200, density=1.0):
        from scipy.interpolate import griddata
        i, j = self.ij
        u0, u1, v0, v1 = self.window
        f = self.alpha > 0
        ny = max(int(n*(v1 - v0)/(u1 - u0)), 20)
        X, Y = np.meshgrid(np.linspace(u0, u1, n), np.linspace(v0, v1, ny))
        u = griddata(self.C[f], U[f, i], (X, Y))
        v = griddata(self.C[f], U[f, j], (X, Y))
        solid = griddata(self.C, self.alpha, (X, Y), method="nearest") < 0.5
        u[solid] = np.nan
        v[solid] = np.nan
        ax.streamplot(X, Y, np.ma.masked_invalid(u), np.ma.masked_invalid(v), color=ff.WALL, linewidth=0.45,
                      density=density, arrowsize=0.5, broken_streamlines=True)


def sphere_section(ax, centre, R, axis, c, ij, **kw):
    d = c - centre[axis]
    if abs(d) >= R:
        return
    r = np.sqrt(R*R - d*d)
    t = np.linspace(0, 2*np.pi, 300)
    ax.plot(centre[ij[0]] + r*np.cos(t), centre[ij[1]] + r*np.sin(t), **{"color": ff.WALL, "lw": 0.8, **kw})


def mag(U):
    return np.linalg.norm(U, axis=1)


# ---------------------------------------------------------------------------------------- figures

def v11(out, r):
    """symmetry plane z = 0: |U| and streamlines at Re 100, 300 (D/h 48), the cut cells at D/h 24"""
    root = Path(r["cases"])
    fig, ax = plt.subplots(1, 3, figsize=(10.0, 3.0), layout="constrained",
                           gridspec_kw={"width_ratios": [1.5, 1.5, 1.0]})
    w = (-1.2, 2.8, 0.0, 1.5)
    for a_, Re in zip(ax[:2], (100, 300)):
        s = Slice(root / ("re%d_dh48" % Re), 2, 1e-9, w)
        U = s.field("U")
        pc = s.draw(a_, mag(U), norm=matplotlib.colors.Normalize(0, 1.15))
        s.streamlines(a_, U, density=0.9)
        sphere_section(a_, (0, 0, 0), 0.5, 2, 0.0, s.ij)
        a_.set_title("Re %d, D/h 48" % Re)
        a_.set_xlabel("x/D")
    ax[0].set_ylabel("y/D")
    wz = (-0.62, -0.22, 0.12, 0.52)
    s = Slice(root / "re100_dh24", 2, 1e-9, wz)
    U = s.field("U")
    s.draw(ax[2], mag(U), norm=matplotlib.colors.Normalize(0, 1.15), edges=True)
    sphere_section(ax[2], (0, 0, 0), 0.5, 2, 0.0, s.ij, lw=1.0)
    ax[2].set_title("cut cells, D/h 24 (Re 100)")
    ax[2].set_xlabel("x/D")
    fig.colorbar(pc, ax=ax, fraction=0.02, pad=0.01).set_label("|U|/U")
    out.save(fig, "V11_flow")


def v14(out, r):
    """3D-1Z (Re 20) on D/h 23.4: mid-plane z = H/2 (|U|, streamlines, cut cells) and the plane
    y = 0.2 through the cylinder axis, where the immersed cylinder meets the channel walls"""
    root = Path(r["cases"])
    case = root / "re20_n24"
    H = 0.41
    fig = plt.figure(figsize=(10.0, 4.6), layout="constrained")
    g = fig.add_gridspec(2, 3, width_ratios=[2.2, 1.0, 1.0], height_ratios=[1, 1])
    a0 = fig.add_subplot(g[0, :])
    s = Slice(case, 2, 0.5*H, (0.25, 1.6, 0.0, H))
    U = s.field("U")
    pc = s.draw(a0, mag(U), norm=matplotlib.colors.Normalize(0, 0.5))
    s.streamlines(a0, U, n=300, density=1.2)
    sphere_section(a0, (0.5, 0.2, 0.0), 0.05, 2, 0.0, s.ij)   # cylinder: the section is the circle
    a0.set_title("3D-1Z, Re 20, D/h 23.4: mid-plane z = H/2")
    a0.set_xlabel("x")
    a0.set_ylabel("y")
    fig.colorbar(pc, ax=a0, fraction=0.02, pad=0.01).set_label("|U|")
    a1 = fig.add_subplot(g[1, 0])
    s = Slice(case, 1, 0.2, (0.3, 0.8, 0.0, H))
    U = s.field("U")
    s.draw(a1, mag(U), norm=matplotlib.colors.Normalize(0, 0.5))
    for xx in (0.45, 0.55):
        a1.axvline(xx, color=ff.WALL, lw=0.8)
    a1.set_title("plane y = 0.2 through the cylinder axis")
    a1.set_xlabel("x")
    a1.set_ylabel("z")
    a2 = fig.add_subplot(g[1, 1])
    s = Slice(case, 2, 0.5*H, (0.43, 0.53, 0.18, 0.28))
    s.draw(a2, mag(s.field("U")), norm=matplotlib.colors.Normalize(0, 0.5), edges=True)
    sphere_section(a2, (0.5, 0.2, 0.0), 0.05, 2, 0.0, s.ij, lw=1.0)
    a2.set_title("cut cells, z = H/2")
    a2.set_xlabel("x")
    a3 = fig.add_subplot(g[1, 2])
    s = Slice(case, 1, 0.2, (0.43, 0.57, 0.0, 0.07))
    s.draw(a3, mag(s.field("U")), norm=matplotlib.colors.Normalize(0, 0.5), edges=True)
    for xx in (0.45, 0.55):
        a3.axvline(xx, color=ff.WALL, lw=1.0)
    a3.set_title("junction with the wall z = 0")
    a3.set_xlabel("x")
    a3.set_ylabel("z")
    out.save(fig, "V14_flow")


FIGURES = {"V11": v11, "V14": v14}


def figures(out, R, only=None):
    for k, fn in FIGURES.items():
        if k not in R or (only and k not in only):
            continue
        try:
            fn(out, R[k])
        except Exception as e:   # noqa: BLE001
            print("  flow figure %s skipped: %s: %s" % (k, type(e).__name__, e))
