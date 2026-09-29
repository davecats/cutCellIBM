"""Cut-cell geometry on a 2-D Cartesian grid, mirroring src/cutCellIBM/cutCellGeometry.

Used only to draw the figures of doc/cutCellIBM.tex from real cut-cell data:
edge crossings by bisection, face fractions theta, cell fractions alpha by
sub-sampling, sliver absorption, closure wall vector S_w, wall distances.
"""
import numpy as np


class Circle:
    """Level set of a cylinder: negative inside the body, positive in the fluid."""

    def __init__(self, xc, yc, R):
        self.xc, self.yc, self.R = xc, yc, R

    def __call__(self, x, y):
        return np.hypot(np.asarray(x) - self.xc, np.asarray(y) - self.yc) - self.R

    def wall_distance(self, x, y):
        return self(x, y)

    def moved(self, dx, dy=0.0):
        return Circle(self.xc + dx, self.yc + dy, self.R)


def bisect(phi, a, b, n=40):
    """Crossing of phi on the segment a-b (phi(a), phi(b) of opposite sign)."""
    a = np.array(a, float)
    b = np.array(b, float)
    fa = phi(*a)
    for _ in range(n):
        m = 0.5 * (a + b)
        fm = phi(*m)
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    return 0.5 * (a + b)


def wet_segment(phi, p0, p1):
    """Fluid part of the edge p0-p1 as (fraction, list of points)."""
    f0, f1 = phi(*p0) > 0, phi(*p1) > 0
    if f0 and f1:
        return 1.0, [np.array(p0, float), np.array(p1, float)]
    if not f0 and not f1:
        return 0.0, []
    x = bisect(phi, p0, p1) if f0 else bisect(phi, p1, p0)
    L = np.hypot(*(np.subtract(p1, p0)))
    if f0:
        return np.hypot(*(x - p0)) / L, [np.array(p0, float), x]
    return np.hypot(*(p1 - x)) / L, [x, np.array(p1, float)]


class Grid:
    """Uniform 2-D grid, cells (i, j) with lower-left corner x0 + (i, j) h."""

    def __init__(self, x0, y0, h, nx, ny):
        self.x0, self.y0, self.h, self.nx, self.ny = x0, y0, h, nx, ny

    def corner(self, i, j):
        return np.array([self.x0 + i * self.h, self.y0 + j * self.h])

    def centre(self, i, j):
        return self.corner(i, j) + 0.5 * self.h


def fractions(g, phi, ns=32, alpha_min=0.02, absorb=True):
    """theta on x-faces (nx+1, ny) and y-faces (nx, ny+1), alpha (nx, ny),
    fluid centroids (nx, ny, 2), absorbed mask."""
    h = g.h
    thx = np.zeros((g.nx + 1, g.ny))
    thy = np.zeros((g.nx, g.ny + 1))
    for i in range(g.nx + 1):
        for j in range(g.ny):
            thx[i, j] = wet_segment(phi, g.corner(i, j), g.corner(i, j + 1))[0]
    for i in range(g.nx):
        for j in range(g.ny + 1):
            thy[i, j] = wet_segment(phi, g.corner(i, j), g.corner(i + 1, j))[0]
    s = (np.arange(ns) + 0.5) / ns
    SX, SY = np.meshgrid(s, s, indexing="ij")
    alpha = np.zeros((g.nx, g.ny))
    cen = np.zeros((g.nx, g.ny, 2))
    for i in range(g.nx):
        for j in range(g.ny):
            c = g.corner(i, j)
            X, Y = c[0] + SX * h, c[1] + SY * h
            wet = phi(X, Y) > 0
            alpha[i, j] = wet.mean()
            if wet.any():
                cen[i, j] = [X[wet].mean(), Y[wet].mean()]
            else:
                cen[i, j] = g.centre(i, j)
    absorbed = np.zeros_like(alpha, bool)
    if absorb:
        absorbed = (alpha > 0) & (alpha < alpha_min)
        for i, j in zip(*np.nonzero(absorbed)):
            alpha[i, j] = 0.0
            thx[i, j] = thx[i + 1, j] = 0.0
            thy[i, j] = thy[i, j + 1] = 0.0
    for i, j in zip(*np.nonzero(alpha == 0)):   # a cell without fluid closes its faces
        thx[i, j] = thx[i + 1, j] = 0.0
        thy[i, j] = thy[i, j + 1] = 0.0
    return thx, thy, alpha, cen, absorbed


def wall_vector(g, thx, thy, i, j):
    """S_w = -sum_f theta_f S_f (outward from the fluid into the body)."""
    h = g.h
    s = np.array([thx[i + 1, j] * h - thx[i, j] * h, thy[i, j + 1] * h - thy[i, j] * h])
    return -s
