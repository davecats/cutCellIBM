#!/usr/bin/env python
"""
Order-of-accuracy study of the cut-cell immersed boundary, on the notebook's
discretisation (2-D Cartesian, periodic box, body = level set). The OpenFOAM
code is not involved; this explores which wall treatments reach second order.

Fluid = interior of a disk of radius R centred in the box (the paper's Poisson
test), so that manufactured solutions are simple polynomials.

Tests (D = 1):
  LD  Laplacian, Dirichlet wall:  lap(w) = 1,         w(R) = 0
  LN  Laplacian, Neumann wall:    lap(w) = 16r^2-8R^2, dw/dn(R) = 0
  CD  convection-diffusion, Dirichlet, u = Omega(-y, x) (tangential to the wall)
      T = (R^2-r^2) x y
  CN  convection-diffusion, Neumann,   T = (R^2-r^2)^2 x y

Dirichlet wall variants:
  V0  d = alpha V/(2 A_wall)                              (current method)
  V1  d = signed distance of the fluid centroid to the wall
  V2  V1 + curvature correction from the cell's own Laplacian (row scaling)
  V3  V1 + deferred three-point quadratic along the wall normal (one probe)
Convective face value variants:
  C0  face value = mean of the two cell values           (current method)
  C1  value reconstructed at the wet-face centroid with cell gradients (deferred)
Diffusive face flux variants:
  F0  theta A (T_N - T_P)/h                              (current method)
  F1  F0 + deferred correction for the true centroid positions (least-squares
      gradients on the fluid centroids, minimum-correction non-orthogonal form)

Usage: order_study.py [--N 20 40 80 160] [--omega 10] [--nSub 32]
"""
import argparse, sys, time
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import splu

ap = argparse.ArgumentParser()
ap.add_argument('--N', type=int, nargs='+', default=[20, 40, 80, 160])
ap.add_argument('--omega', type=float, default=10.0)
ap.add_argument('--nSub', type=int, default=32)
ap.add_argument('--alphaMin', type=float, default=0.02)
ap.add_argument('--tests', nargs='+', default=['LD', 'LN', 'CD', 'CN'])
ap.add_argument('--noProject', action='store_true', help='do not make the Neumann source compatible with the discrete fluid volume')
args = ap.parse_args()

R, L = 1.0, 3.0
ALPHA_MIN = args.alphaMin


def levelSet(x, y):
    """Signed distance, positive in the fluid (inside the disk)."""
    return R - np.hypot(x - L/2, y - L/2)


def gradLevelSet(x, y):
    """Unit gradient of the level set (points into the fluid)."""
    dx, dy = x - L/2, y - L/2
    r = np.maximum(np.hypot(dx, dy), 1e-300)
    return -dx/r, -dy/r


# ------------------------------------------------------------------ geometry
class Geometry:
    def __init__(self, N, nSub):
        self.N = N
        h = L/N
        self.h = h
        self.V = h*h
        xp = (np.arange(N) + 0.5)*h
        self.XP, self.YP = np.meshgrid(xp, xp)
        xf = np.arange(N + 1)*h
        # face fractions and crossing positions by bisection
        def frac(x0, y0, x1, y1):
            f0, f1 = levelSet(x0, y0), levelSet(x1, y1)
            th = (f0 > 0).astype(float)
            t = np.full(f0.shape, np.nan)            # crossing parameter from (x0,y0)
            cut = (f0 > 0) != (f1 > 0)
            if cut.any():
                lo = np.zeros(cut.sum()); hi = np.ones(cut.sum())
                xa, ya, sx, sy, fa = x0[cut], y0[cut], (x1 - x0)[cut], (y1 - y0)[cut], f0[cut]
                for _ in range(50):
                    mid = 0.5*(lo + hi)
                    same = (levelSet(xa + mid*sx, ya + mid*sy) > 0) == (fa > 0)
                    lo = np.where(same, mid, lo); hi = np.where(same, hi, mid)
                tc = 0.5*(lo + hi)
                t[cut] = tc
                th[cut] = np.where(fa > 0, tc, 1 - tc)
            return th, t, cut, f0 > 0
        XF1 = np.tile(xf, (N, 1)); Y0 = np.tile(xf[:-1][:, None], (1, N + 1))
        self.thX, self.tX, self.cutX, self.fl0X = frac(XF1, Y0, XF1, Y0 + h)   # (N, N+1)
        YF1 = np.tile(xf[:, None], (1, N)); X0 = np.tile(xf[:-1], (N + 1, 1))
        self.thY, self.tY, self.cutY, self.fl0Y = frac(X0, YF1, X0 + h, YF1)   # (N+1, N)
        # wet-face centroids (segment midpoint of the wet part)
        # x-face at (xf[i], y from yf[j] to yf[j]+h): wet part is [0,t] if fluid at 0 else [t,1]
        sX = np.where(self.cutX, np.where(self.fl0X, 0.5*self.tX, 0.5*(1 + self.tX)), 0.5)
        self.cwX = (XF1, Y0 + sX*h)
        sY = np.where(self.cutY, np.where(self.fl0Y, 0.5*self.tY, 0.5*(1 + self.tY)), 0.5)
        self.cwY = (X0 + sY*h, YF1)
        # cell fraction and fluid centroid by sub-sampling
        sd = levelSet(self.XP, self.YP)
        cand = np.abs(sd) < 0.5*np.hypot(h, h)
        alpha = (sd > 0).astype(float)
        cx, cy = self.XP.copy(), self.YP.copy()
        if cand.any():
            x0 = self.XP[cand] - 0.5*h; y0 = self.YP[cand] - 0.5*h
            s = (np.arange(nSub) + 0.5)/nSub
            acc = np.zeros(x0.shape); sx = np.zeros(x0.shape); sy = np.zeros(x0.shape)
            for a in s:
                for b in s:
                    inside = levelSet(x0 + a*h, y0 + b*h) > 0
                    acc += inside; sx += inside*(x0 + a*h); sy += inside*(y0 + b*h)
            alpha[cand] = acc/(nSub*nSub)
            ok = acc > 0
            cx[cand] = np.where(ok, sx/np.maximum(acc, 1), self.XP[cand])
            cy[cand] = np.where(ok, sy/np.maximum(acc, 1), self.YP[cand])
        tiny = alpha < ALPHA_MIN
        alpha[tiny] = 0
        self.thX[:, :-1][tiny] = 0; self.thX[:, 1:][tiny] = 0
        self.thY[:-1, :][tiny] = 0; self.thY[1:, :][tiny] = 0
        self.alpha = alpha
        self.alphaV = alpha*self.V
        self.fluid = alpha > 0
        self.cx, self.cy = cx, cy
        # wall segment from closure (points into the body), and distances
        self.Swx = -(self.thX[:, 1:] - self.thX[:, :-1])*h
        self.Swy = -(self.thY[1:, :] - self.thY[:-1, :])*h
        self.Awall = np.hypot(self.Swx, self.Swy)
        self.wall = self.fluid & (self.Awall > 1e-12*h)
        self.d0 = np.where(self.wall, self.alphaV/(2*np.maximum(self.Awall, 1e-300)), np.inf)
        self.dc = np.where(self.wall, np.maximum(levelSet(cx, cy), 1e-3*h), np.inf)
        self.nCut = int((self.fluid & (alpha < 1)).sum())


# ---------------------------------------------------------------- operators
def assemble(G, wallType, variant, omega, T_old=None, convVariant='C0', wallValue=0.0, faceVariant='F0'):
    """
    Sparse matrix and right-hand side of  -lap(T) [+ div(u T)] on the fluid cells,
    volume-integrated. Returns (A, b_explicit) where b_explicit holds the deferred
    corrections evaluated from T_old (zero if T_old is None).
    """
    N, h, V = G.N, G.h, G.V
    ID = np.arange(N*N).reshape(N, N)
    E, W, No, S = np.roll(ID, -1, 1), np.roll(ID, 1, 1), np.roll(ID, -1, 0), np.roll(ID, 1, 0)
    rows, cols, vals = [], [], []
    diag = np.zeros((N, N))
    b = np.zeros((N, N))

    def add(r, c, v):
        rows.append(r.ravel()); cols.append(c.ravel()); vals.append(v.ravel())

    # diffusion: -sum_f theta_f (T_N - T_P)/h * h  (D = 1); with F1 the implicit
    # coefficient uses the centroid distance along the face normal (over-relaxed form)
    if faceVariant == 'F1':
        ndX, ndY = centroidDistances(G)           # (N, N+1), (N+1, N): n.d per face
        cX = G.thX*h/ndX; cY = G.thY*h/ndY
        gE, gW = cX[:, 1:], cX[:, :-1]
        gN, gS = cY[1:, :], cY[:-1, :]
    else:
        gE, gW = G.thX[:, 1:], G.thX[:, :-1]
        gN, gS = G.thY[1:, :], G.thY[:-1, :]
    diag += gE + gW + gN + gS
    add(ID, E, -gE); add(ID, W, -gW); add(ID, No, -gN); add(ID, S, -gS)

    # convection: u = omega(-(y-yc), x-xc), flux through the wet part of each face
    # evaluated at the wet-face centroid (exact for a linear field)
    if omega != 0:
        ux = lambda x, y: -omega*(y - L/2)
        uy = lambda x, y: omega*(x - L/2)
        phiX = G.thX*h*ux(*G.cwX)           # (N, N+1), positive in +x
        phiY = G.thY*h*uy(*G.cwY)
        fE, fW, fN, fS = phiX[:, 1:], phiX[:, :-1], phiY[1:, :], phiY[:-1, :]
        # central: T_f = 0.5 (T_P + T_N)
        diag += 0.5*(fE - fW + fN - fS)
        add(ID, E, 0.5*fE); add(ID, W, -0.5*fW); add(ID, No, 0.5*fN); add(ID, S, -0.5*fS)
        if convVariant == 'C1' and T_old is not None:
            # deferred: value at the wet-face centroid from a linear reconstruction
            # about the fluid centroids of the two cells, gradient by cut-cell Gauss
            gx, gy = cutCellGradient(G, T_old)
            def faceCorr(TP, TN, gxP, gyP, gxN, gyN, cxP, cyP, cxN, cyN, xw, yw):
                recP = TP + gxP*(xw - cxP) + gyP*(yw - cyP)
                recN = TN + gxN*(xw - cxN) + gyN*(yw - cyN)
                return 0.5*(recP + recN) - 0.5*(TP + TN)
            # x faces between W-cell (i-1) and cell i, including periodic wrap
            TW = np.roll(T_old, 1, 1); gxW = np.roll(gx, 1, 1); gyW = np.roll(gy, 1, 1)
            cxW = np.roll(G.cx, 1, 1) - np.where(np.arange(N) == 0, L, 0)[None, :]
            cyW = np.roll(G.cy, 1, 1)
            xw, yw = G.cwX[0][:, :-1], G.cwX[1][:, :-1]
            corrX = faceCorr(T_old, TW, gx, gy, gxW, gyW, G.cx, G.cy, cxW, cyW, xw, yw)  # face W of cell
            corrX = np.concatenate([corrX, corrX[:, :1]], 1)   # (N, N+1)
            TS = np.roll(T_old, 1, 0); gxS = np.roll(gx, 1, 0); gyS = np.roll(gy, 1, 0)
            cxS = np.roll(G.cx, 1, 0)
            cyS = np.roll(G.cy, 1, 0) - np.where(np.arange(N) == 0, L, 0)[:, None]
            xw, yw = G.cwY[0][:-1, :], G.cwY[1][:-1, :]
            corrY = faceCorr(T_old, TS, gx, gy, gxS, gyS, G.cx, G.cy, cxS, cyS, xw, yw)
            corrY = np.concatenate([corrY, corrY[:1, :]], 0)
            # explicit flux correction moved to the right-hand side: -(sum_f phi_f corr_f)
            b -= (phiX[:, 1:]*corrX[:, 1:] - phiX[:, :-1]*corrX[:, :-1]
                  + phiY[1:, :]*corrY[1:, :] - phiY[:-1, :]*corrY[:-1, :])

    # Dirichlet wall: -A_wall (T_w - T_P)/d  added to -lap
    if wallType == 'D':
        d = G.d0 if variant == 'V0' else G.dc
        c = np.where(G.wall, G.Awall/d, 0.0)
        diag += c
        b += c*wallValue
        if variant == 'V2':
            # two-point gradient sits at d/2 from the wall: (T_P - T_w)/d = T'_w + (d/2) T''.
            # Take T'' from the cell's own Laplacian: -lap = (F_faces + F_wall)/(alpha V)
            # => F_wall (1 + k) = A (T_P-T_w)/d - k F_faces, k = A d/(2 alpha V), i.e. the
            # whole diffusive row is divided by (1 + k).
            k = np.zeros_like(d); k[G.wall] = (G.Awall*d)[G.wall]/(2*G.alphaV[G.wall])
            scale = 1.0/(1.0 + k)
            # rescale the diffusive part only: rebuild diffusive entries with the scale
            diag -= (gE + gW + gN + gS) + c
            diag += ((gE + gW + gN + gS) + c)*scale
            b -= c*wallValue; b += c*wallValue*scale
            add(ID, E, -gE*(scale - 1)); add(ID, W, -gW*(scale - 1))
            add(ID, No, -gN*(scale - 1)); add(ID, S, -gS*(scale - 1))
        if variant == 'V3' and T_old is not None:
            b += threePointCorrection(G, T_old, wallValue)

    if faceVariant == 'F1' and T_old is not None:
        b += faceFluxCorrection(G, T_old)

    # solid cells: identity rows
    diag[~G.fluid] = 1.0
    b[~G.fluid] = 0.0
    add(ID, ID, diag)
    A = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                   shape=(N*N, N*N)).tocsr()
    return A, b


def cutCellGradient(G, T):
    """Volume-averaged gradient of the fluid part: sum_f theta_f (T_f - T_P) S_f / (alpha V)."""
    h = G.h
    dE = 0.5*(np.roll(T, -1, 1) - T)*G.thX[:, 1:]
    dW = 0.5*(np.roll(T, 1, 1) - T)*G.thX[:, :-1]
    dN = 0.5*(np.roll(T, -1, 0) - T)*G.thY[1:, :]
    dS = 0.5*(np.roll(T, 1, 0) - T)*G.thY[:-1, :]
    aV = np.maximum(G.alphaV, 1e-300)
    return h*(dE - dW)/aV, h*(dN - dS)/aV


def lsqGradient(G, T):
    """Weighted least-squares gradient of the cell values about the fluid centroids,
    using the four face neighbours that are fluid (exact for a linear field)."""
    N = G.N
    a11 = np.zeros_like(T); a12 = np.zeros_like(T); a22 = np.zeros_like(T)
    r1 = np.zeros_like(T); r2 = np.zeros_like(T)
    for axis, shift, sign in ((1, -1, 1), (1, 1, -1), (0, -1, 1), (0, 1, -1)):
        TN = np.roll(T, shift, axis)
        cxN = np.roll(G.cx, shift, axis); cyN = np.roll(G.cy, shift, axis)
        flN = np.roll(G.fluid, shift, axis)
        wrap = np.zeros(N); wrap[-1 if shift == -1 else 0] = sign*L
        if axis == 1:
            cxN = cxN + wrap[None, :]
        else:
            cyN = cyN + wrap[:, None]
        dx = cxN - G.cx; dy = cyN - G.cy
        w = np.where(flN, 1.0/np.maximum(dx*dx + dy*dy, 1e-300), 0.0)
        dT = TN - T
        a11 += w*dx*dx; a12 += w*dx*dy; a22 += w*dy*dy
        r1 += w*dx*dT; r2 += w*dy*dT
    det = a11*a22 - a12*a12
    ok = det > 1e-12*np.maximum(a11*a22, 1e-300)
    gx = np.where(ok, (a22*r1 - a12*r2)/np.where(ok, det, 1), 0.0)
    gy = np.where(ok, (a11*r2 - a12*r1)/np.where(ok, det, 1), 0.0)
    return gx, gy


def centroidDistances(G):
    """n.d of every face: distance between the fluid centroids of the two cells
    along the face normal (h for uncut cells). Faces of solid cells keep h."""
    N, h = G.N, G.h
    dX = np.roll(G.cx, -1, 1) - G.cx; dX[:, -1] += L      # cell i -> i+1, periodic
    dY = np.roll(G.cy, -1, 0) - G.cy; dY[-1, :] += L
    both = G.fluid & np.roll(G.fluid, -1, 1)
    dX = np.where(both, np.maximum(dX, 0.1*h), h)
    bothY = G.fluid & np.roll(G.fluid, -1, 0)
    dY = np.where(bothY, np.maximum(dY, 0.1*h), h)
    ndX = np.concatenate([dX[:, -1:], dX], 1)             # face i between cells i-1 and i
    ndY = np.concatenate([dY[-1:, :], dY], 0)
    return ndX, ndY


def faceFluxCorrection(G, T):
    """
    Deferred correction of the diffusive flux through the wet part of a face
    (over-relaxed non-orthogonal form on the fluid centroids): with d the
    centroid-to-centroid vector and g_f the mean least-squares gradient,
    n.grad T = (T_N - T_P)/(n.d) + g_f.n - (g_f.d)/(n.d). The first term is
    implicit (see assemble), the rest is returned here as a per-cell source of
    -lap (volume-integrated). Exact for a linear field.
    """
    N, h = G.N, G.h
    gx, gy = lsqGradient(G, T)
    src = np.zeros_like(T)
    for axis, th, nxf, nyf in ((1, G.thX, 1.0, 0.0), (0, G.thY, 0.0, 1.0)):
        TN = np.roll(T, -1, axis); gxN = np.roll(gx, -1, axis); gyN = np.roll(gy, -1, axis)
        cxN = np.roll(G.cx, -1, axis); cyN = np.roll(G.cy, -1, axis)
        wrap = np.zeros(N); wrap[-1] = L
        if axis == 1:
            cxN = cxN + wrap[None, :]
        else:
            cyN = cyN + wrap[:, None]
        dx = cxN - G.cx; dy = cyN - G.cy
        d2 = np.maximum(dx*dx + dy*dy, 1e-300)
        thf = th[:, 1:] if axis == 1 else th[1:, :]
        both = G.fluid & np.roll(G.fluid, -1, axis) & (thf > 0)
        gfx = 0.5*(gx + gxN); gfy = 0.5*(gy + gyN)
        nd = nxf*dx + nyf*dy
        # over-relaxed: implicit (T_N - T_P)/(n.d) already in the matrix; explicit remainder
        ndp = np.where(both, np.maximum(nd, 0.1*h), h)
        expl = (nxf*gfx + nyf*gfy) - (gfx*dx + gfy*dy)/ndp
        corr = np.where(both, thf*h*expl, 0.0)                  # extra flux out of P through f
        src += corr                        # -lap row of P: -flux_out; explicit part to rhs
        src -= np.roll(corr, 1, axis)      # neighbour receives the flux
    return src


def bilinear(G, T, x, y):
    """Bilinear interpolation of cell values (at cell centres) at points (x, y);
    returns nan where any of the four cells is solid."""
    N, h = G.N, G.h
    fx = x/h - 0.5; fy = y/h - 0.5
    i0 = np.floor(fx).astype(int); j0 = np.floor(fy).astype(int)
    tx = fx - i0; ty = fy - j0
    def val(i, j):
        i = i % N; j = j % N
        v = T[j, i].astype(float).copy()
        v[~G.fluid[j, i]] = np.nan
        return v
    return ((1 - tx)*(1 - ty)*val(i0, j0) + tx*(1 - ty)*val(i0 + 1, j0)
            + (1 - tx)*ty*val(i0, j0 + 1) + tx*ty*val(i0 + 1, j0 + 1))


def threePointCorrection(G, T, wallValue):
    """
    Explicit correction of the wall flux: quadratic through the wall value, the
    cell value at distance d_c and a probe at distance d_c + h along the wall
    normal, minus the implicit two-point estimate. Cells whose probe touches a
    solid cell keep the two-point flux.
    """
    h = G.h
    w = G.wall
    corr = np.zeros_like(T)
    if not w.any():
        return corr
    # unit normal pointing into the fluid, from the level set at the centroid
    nx, ny = gradLevelSet(G.cx[w], G.cy[w])
    # wall point along the normal through the centroid
    dc = G.dc[w]
    xw = G.cx[w] - dc*nx; yw = G.cy[w] - dc*ny
    a = dc; bq = dc + h
    xp = xw + bq*nx; yp = yw + bq*ny
    T2 = bilinear(G, T, xp, yp)
    T1 = T[w]
    grad3 = (-wallValue*(a + bq)/(a*bq) + T1*bq/(a*(bq - a)) - T2*a/(bq*(bq - a)))
    grad2 = (T1 - wallValue)/a
    dgrad = np.where(np.isfinite(grad3), grad3 - grad2, 0.0)
    # wall flux of -lap is -A (dT/dn)_wall with n into the fluid: -A grad; the implicit
    # part already carries -A grad2 on the row, so the deferred part is -A (grad3 - grad2)
    # moved to the right-hand side: b += A (grad3 - grad2)
    corr[w] = G.Awall[w]*dgrad
    return corr


# --------------------------------------------------------------------- tests
def exactAndSource(test, x, y, omega):
    X, Y = x - L/2, y - L/2
    r2 = X*X + Y*Y
    if test == 'LD':
        return -(R*R - r2)/4, np.ones_like(x)                       # lap w = 1
    if test == 'LN':
        return r2*r2 - 2*R*R*r2, 16*r2 - 8*R*R
    if test == 'CD':
        T = (R*R - r2)*X*Y
        lap = -12*X*Y
        conv = omega*(X*X - Y*Y)*(R*R - r2)
        return T, conv - lap                                       # u.grad T - lap T = q
    if test == 'CN':
        T = (R*R - r2)**2*X*Y
        lap = X*Y*(32*r2 - 24*R*R)
        conv = omega*(R*R - r2)**2*(X*X - Y*Y)
        return T, conv - lap
    raise ValueError(test)


def solve(G, test, variant, convVariant, omega, faceVariant='F0'):
    wallType = 'D' if test in ('LD', 'CD') else 'N'
    om = omega if test in ('CD', 'CN') else 0.0
    Tex, q = exactAndSource(test, G.cx, G.cy, om)     # source at the fluid centroid
    # -lap T + div(uT) = -q for LD/LN (lap w = q), or = q for CD/CN (u.grad T - lap T = q)
    rhs0 = (-q if test in ('LD', 'LN') else q)*G.alphaV
    if wallType == 'N' and not args.noProject:
        # pure Neumann problem: the discrete fluid volume differs from the disk by
        # the absorbed slivers, so the source must be projected onto the compatible
        # space (as the pressure equation's right-hand side is, by construction);
        # otherwise the reference pinning deposits the mismatch as a point source
        rhs0 = rhs0 - G.alphaV*np.sum(rhs0[G.fluid])/np.sum(G.alphaV[G.fluid])
    deferred = (variant == 'V3' and wallType == 'D') or (convVariant == 'C1' and om != 0) or faceVariant == 'F1'
    T = np.zeros_like(Tex)
    A, b = assemble(G, wallType, variant, om, None, convVariant, faceVariant=faceVariant)
    if wallType == 'N':
        # pure Neumann: pin one fluid cell (reference doubling, as fvMatrix::setReference)
        iref = np.flatnonzero(G.fluid.ravel())[0]
        A = A.tolil(); A[iref, iref] *= 2; A = A.tocsr()
    lu = splu(A.tocsc())
    nIter = 1
    for it in range(60):
        A_, b = assemble(G, wallType, variant, om, T if deferred else None, convVariant, faceVariant=faceVariant)
        r = rhs0 + b
        r[~G.fluid] = 0
        Tn = lu.solve(r.ravel()).reshape(G.N, G.N)
        change = np.max(np.abs(Tn - T)[G.fluid])
        T = Tn
        nIter = it + 1
        if not deferred or change < 1e-11*max(1.0, np.max(np.abs(T))):
            break
    err = T - Tex
    if wallType == 'N':
        m = np.sum(err*G.alphaV)/np.sum(G.alphaV)
        err = err - m
    fl = G.fluid
    L2 = np.sqrt(np.sum((err**2*G.alphaV)[fl])/np.sum(G.alphaV[fl]))
    Linf = np.max(np.abs(err[fl]))
    return L2, Linf, nIter


def order(errs, Ns):
    return [np.nan] + [np.log(errs[i-1]/errs[i])/np.log(Ns[i]/Ns[i-1]) for i in range(1, len(errs))]


if __name__ == '__main__':
    geos = {}
    for N in args.N:
        t0 = time.time()
        geos[N] = Geometry(N, args.nSub)
        print(f'geometry N={N}: {geos[N].nCut} cut cells, {time.time()-t0:.1f}s', flush=True)
    configs = {
        'LD': [('V0', 'C0', 'F0'), ('V1', 'C0', 'F0'), ('V1', 'C0', 'F1'), ('V2', 'C0', 'F1'), ('V3', 'C0', 'F1')],
        'LN': [('N0', 'C0', 'F0'), ('N0', 'C0', 'F1')],
        'CD': [('V0', 'C0', 'F0'), ('V1', 'C1', 'F1'), ('V2', 'C1', 'F1'), ('V3', 'C0', 'F1'), ('V3', 'C1', 'F1')],
        'CN': [('N0', 'C0', 'F0'), ('N0', 'C0', 'F1'), ('N0', 'C1', 'F1')],
    }
    for test in args.tests:
        print(f'\n=== {test}  (omega = {args.omega if test in ("CD", "CN") else 0})')
        print(f'{"variant":10s} ' + ' '.join(f'{"N=%d" % N:>22s}' for N in args.N) + '   iterations')
        for v, cv, fv in configs[test]:
            L2s, Lis, its = [], [], []
            for N in args.N:
                L2, Li, it = solve(geos[N], test, v, cv, args.omega, fv)
                L2s.append(L2); Lis.append(Li); its.append(it)
            oL2 = order(L2s, args.N); oLi = order(Lis, args.N)
            cells = ' '.join(f'{L2s[i]:9.2e} ({oL2[i]:4.2f}) {Lis[i]:9.2e}({oLi[i]:4.2f})'.replace('( nan)', '(  - )')
                             for i in range(len(args.N)))
            print(f'{v+"/"+cv+"/"+fv:10s} {cells}   {its}', flush=True)
        print('    (per mesh: L2 error (order)  Linf error (order); errors of the fluid cells, alpha-weighted)')
