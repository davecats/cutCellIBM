#!/usr/bin/env python
"""
Reference solution: the cut-cell IBM of .asset/fvm/SIMPLEC_faceArea_noSlip.ipynb run as a
script. The numerics are copied verbatim from the notebook (cells 1-6); only the
plotting is removed and the geometry, solution and forcing are saved to an .npz
for comparison with the OpenFOAM implementation.

Usage: reference_notebook.py [--nx 81] [--nIter 400] [--bc noSlip] [--out ref.npz]
"""
import os, sys, argparse, time
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu

ap = argparse.ArgumentParser()
ap.add_argument('--nx', type=int, default=81)
ap.add_argument('--nIter', type=int, default=400)
ap.add_argument('--bc', default='noSlip')
ap.add_argument('--alphaMin', type=float, default=0.02)
ap.add_argument('--dtMax', type=float, default=10.0)
ap.add_argument('--out', default='ref.npz')
args = ap.parse_args()

# Inputs (notebook cell 1)
R  = 1
Lx = 10
Ly = 10
nx = args.nx
ny = args.nx
dt = 0.1
nu = 1.0
uBulkTarget = 0.5
bcType = args.bc

# Grid (cell 2)
IU, IV, IP = 0, 1, 2
dx, dy = Lx/nx, Ly/ny
V      = dx*dy
nCells = nx*ny
xp = (np.arange(nx)+0.5)*dx
yp = (np.arange(ny)+0.5)*dy
xf = np.arange(nx+1)*dx
yf = np.arange(ny+1)*dy
XP, YP = np.meshgrid(xp, yp)
XF, YF = np.meshgrid(xf, yf)
VAR = np.zeros((3, ny, nx))
ID  = np.arange(nCells).reshape(ny, nx)
ROW = np.stack([ID]*5, -1).ravel()
COL = np.stack([ID, np.roll(ID, -1, 1), np.roll(ID, 1, 1),
                    np.roll(ID, -1, 0), np.roll(ID, 1, 0)], -1).ravel()
_pattern    = csr_matrix((np.arange(1.0, ROW.size+1), (ROW, COL)), shape=(nCells, nCells))
CSR_PERM    = (_pattern.data-1).astype(np.int64)
CSR_INDICES = _pattern.indices
CSR_INDPTR  = _pattern.indptr

def zeroStencil():
    return np.zeros((5, ny, nx))
def toCSR(a):
    return csr_matrix((a.reshape(5, -1).T.ravel()[CSR_PERM], CSR_INDICES, CSR_INDPTR), shape=(nCells, nCells))
def offDiagDot(a, phi):
    return (a[1]*np.roll(phi, -1, 1) + a[2]*np.roll(phi, 1, 1)
          + a[3]*np.roll(phi, -1, 0) + a[4]*np.roll(phi, 1, 0))
def wrapX(f):
    return np.concatenate([f, f[:, :1]], axis=1)
def wrapY(f):
    return np.concatenate([f, f[:1, :]], axis=0)

# Operators (cell 4)
def fvm_ddt_LHS(a, dt):
    a[0] += alphaV/dt
def fvm_ddt_RHS(rhs, phi, dt):
    rhs += alphaV*phi/dt
def fvm_lap(a, GammaX, GammaY=None):
    if GammaY is None:
        GammaY = GammaX
    gE = 0.5*(GammaX+np.roll(GammaX, -1, 1))*(dy/dx)*thX[:, 1:]
    gW = 0.5*(GammaX+np.roll(GammaX,  1, 1))*(dy/dx)*thX[:, :-1]
    gN = 0.5*(GammaY+np.roll(GammaY, -1, 0))*(dx/dy)*thY[1:, :]
    gS = 0.5*(GammaY+np.roll(GammaY,  1, 0))*(dx/dy)*thY[:-1, :]
    a[0] -= gE+gW+gN+gS
    a[1] += gE
    a[2] += gW
    a[3] += gN
    a[4] += gS
def fvm_divConvFlux(a, phiX, phiY):
    fE, fW = phiX[:, 1:], phiX[:, :-1]
    fN, fS = phiY[1:, :], phiY[:-1, :]
    a[0] += 0.5*(fE-fW+fN-fS)
    a[1] += 0.5*fE
    a[2] -= 0.5*fW
    a[3] += 0.5*fN
    a[4] -= 0.5*fS
def fvc_grd(phi):
    dE = 0.5*(np.roll(phi, -1, 1)-phi)*thX[:, 1:]
    dW = 0.5*(np.roll(phi,  1, 1)-phi)*thX[:, :-1]
    dN = 0.5*(np.roll(phi, -1, 0)-phi)*thY[1:, :]
    dS = 0.5*(np.roll(phi,  1, 0)-phi)*thY[:-1, :]
    return np.stack([dy*(dE-dW), dx*(dN-dS)])
def fvc_grdFlux(phi, GammaX, GammaY=None):
    if GammaY is None:
        GammaY = GammaX
    fx = (phi-np.roll(phi, 1, 1))/dx*0.5*(GammaX+np.roll(GammaX, 1, 1))*dy
    fy = (phi-np.roll(phi, 1, 0))/dy*0.5*(GammaY+np.roll(GammaY, 1, 0))*dx
    return wrapX(fx)*thX, wrapY(fy)*thY
def fvc_divFlux(phiX, phiY):
    return (phiX[:, 1:]-phiX[:, :-1]) + (phiY[1:, :]-phiY[:-1, :])
def interpolatePhi(Ux, Uy):
    return (0.5*(wrapX(np.roll(Ux, 1, 1))+wrapX(Ux))*dy*thX,
            0.5*(wrapY(np.roll(Uy, 1, 0))+wrapY(Uy))*dx*thY)

# Cut cells (cell 5)
def levelSet(x, y):
    return (x-Lx/2)**2 + (y-Ly/2)**2 - R**2
ALPHA_MIN = args.alphaMin

def segmentFluidFraction(pA, pB, nBisect=40):
    fa, fb = pA(0.0), pB(1.0)
    lo = np.zeros_like(fa)
    hi = np.ones_like(fa)
    for _ in range(nBisect):
        mid = 0.5*(lo+hi)
        fm = pA(mid)
        sameAsA = (fm > 0.0) == (fa > 0.0)
        lo = np.where(sameAsA, mid, lo)
        hi = np.where(sameAsA, hi, mid)
    tCross = 0.5*(lo+hi)
    frac = np.where(fa > 0.0, tCross, 1.0-tCross)
    both = (fa > 0.0) == (fb > 0.0)
    return np.where(both, np.where(fa > 0.0, 1.0, 0.0), np.clip(frac, 0.0, 1.0))

def faceFractions():
    XF1 = np.tile(xf, (ny, 1))
    Y0 = np.tile(yf[:-1][:, None], (1, nx+1))
    Y1 = np.tile(yf[1:][:, None], (1, nx+1))
    thX = segmentFluidFraction(lambda t: levelSet(XF1, Y0+t*(Y1-Y0)),
                               lambda t: levelSet(XF1, Y0+t*(Y1-Y0)))
    YF1 = np.tile(yf[:, None], (1, nx))
    X0 = np.tile(xp*0 + xf[:-1], (ny+1, 1))
    X1 = np.tile(xp*0 + xf[1:], (ny+1, 1))
    thY = segmentFluidFraction(lambda t: levelSet(X0+t*(X1-X0), YF1),
                               lambda t: levelSet(X0+t*(X1-X0), YF1))
    thX[:, -1] = thX[:, 0]
    thY[-1, :] = thY[0, :]
    return thX, thY

def cellFraction(nSub=32):
    s = (np.arange(nSub)+0.5)/nSub
    acc = np.zeros((ny, nx))
    for a in s:
        for b in s:
            acc += (levelSet(XF[:-1, :-1]+a*dx, YF[:-1, :-1]+b*dy) > 0.0)
    return acc/(nSub*nSub)

thX, thY = faceFractions()
alpha = cellFraction()
tiny = alpha < ALPHA_MIN
alpha = np.where(tiny, 0.0, alpha)
thX[:, :-1] = np.where(tiny, 0.0, thX[:, :-1]); thX[:, 1:] = np.where(tiny, 0.0, thX[:, 1:])
thY[:-1, :] = np.where(tiny, 0.0, thY[:-1, :]); thY[1:, :] = np.where(tiny, 0.0, thY[1:, :])
fluid  = alpha > 0.0
solid  = ~fluid
alphaV = alpha*V
nearBody = fluid & (np.roll(solid, 1, 1) | np.roll(solid, -1, 1) |
                    np.roll(solid, 1, 0) | np.roll(solid, -1, 0))
IBM_SOLID_WEIGHT = 1e12
Swx   = (thX[:, 1:]-thX[:, :-1])*dy
Swy   = (thY[1:, :]-thY[:-1, :])*dx
Awall = np.hypot(Swx, Swy)
dWall = np.where(Awall > 0.0, alphaV/(2.0*np.maximum(Awall, 1e-300)), np.inf)

def ibm_noSlipWall(a, Gamma):
    a[0] -= np.where(Awall > 0.0, Gamma*Awall/dWall, 0.0)
def ibm_blankSolid(a):
    a[0] += np.where(solid, IBM_SOLID_WEIGHT, 0.0)

# Run controls (cell 6)
consistent  = True
nIterations = args.nIter
cflMax      = 100.0
foMax       = 65.0
dtMax       = args.dtMax
dtRelax     = 1.0
residOld    = np.inf
dtCflPrev   = 1.0
VAR[...] = 0.0
phiHbyAx, phiHbyAy = interpolatePhi(VAR[IU], VAR[IV])

lapl = zeroStencil()
fvm_lap(lapl, nu*np.ones((ny, nx)))
if bcType == 'noSlip':
    ibm_noSlipWall(lapl, nu)
elif bcType != 'freeSlip':
    raise ValueError("bcType must be 'noSlip' or 'freeSlip'")

jRef, iRef = divmod(int(np.flatnonzero(fluid.ravel())[0]), nx)
Afluid     = alphaV.sum()
unitSource = np.stack([alphaV, np.zeros((ny, nx))])
t          = 0.0
def bulkVelocity(u):
    return np.sum(u*alphaV)/Afluid

t0 = time.time()
gx = 0.0
history = []
for it in range(nIterations):
    VARold  = VAR.copy()
    phiOldX, phiOldY = phiHbyAx.copy(), phiHbyAy.copy()
    aMom = []
    for iC in range(2):
        a = zeroStencil()
        fvm_ddt_LHS(a, dt)
        a -= lapl
        aMom.append(a)
    for a in aMom:
        fvm_divConvFlux(a, phiHbyAx, phiHbyAy)
    for a in aMom:
        ibm_blankSolid(a)
    aP     = np.stack([a[0] for a in aMom])
    rAU    = 1.0/aP
    sumOff = np.stack([a[1:].sum(axis=0) for a in aMom])
    aP     = np.maximum(aP, alphaV/dt)
    rAU    = 1.0/aP
    rAtU   = 1.0/np.maximum(aP+sumOff, alphaV/dt) if consistent else rAU
    rAtUV  = rAtU*V
    pStencil = zeroStencil()
    fvm_lap(pStencil, rAtUV[IU], rAtUV[IV])
    pStencil[0][solid] = -1.0
    pStencil[0, jRef, iRef] *= 2.0
    Ap = toCSR(pStencil)
    luMom = [splu(toCSR(a).tocsc()) for a in aMom]
    luP   = splu(Ap.tocsc())

    def advance(sourceField, pOld):
        RHS   = sourceField.reshape(2, -1)
        gradP = fvc_grd(pOld)
        U = np.stack([luMom[iC].solve(RHS[iC]-gradP[iC].ravel()).reshape(ny, nx)
                      for iC in range(2)])
        HbyA = np.stack([(sourceField[iC] - offDiagDot(aMom[iC], U[iC]))*rAU[iC]
                         for iC in range(2)])
        pX, pY = interpolatePhi(HbyA[IU], HbyA[IV])
        if consistent:
            fX, fY = fvc_grdFlux(pOld, (rAtU[IU]-rAU[IU])*V, (rAtU[IV]-rAU[IV])*V)
            pX, pY = pX+fX, pY+fY
            HbyA = HbyA - (rAU-rAtU)*gradP
        div0 = fvc_divFlux(pX, pY)
        rhsP = div0.copy()
        rhsP[solid] = 0.0
        pNew = luP.solve(rhsP.ravel())
        pNew = pNew.reshape(ny, nx)
        fX, fY = fvc_grdFlux(pNew, rAtUV[IU], rAtUV[IV])
        return (HbyA - rAtU*fvc_grd(pNew), pNew, pX-fX, pY-fY, div0)

    rhs = np.zeros((2, ny, nx))
    for iC in range(2):
        fvm_ddt_RHS(rhs[iC], VAR[iC], dt)
    base = advance(rhs,        VAR[IP])
    resp = advance(unitSource, np.zeros((ny, nx)))
    gx   = (uBulkTarget-bulkVelocity(base[0][IU]))/bulkVelocity(resp[0][IU])
    VAR[0:2] = base[0] + gx*resp[0]
    VAR[IP]  = base[1] + gx*resp[1]
    phiHbyAx = base[2] + gx*resp[2]
    phiHbyAy = base[3] + gx*resp[3]
    meanDiv0 = np.mean(np.abs(base[4] + gx*resp[4]))

    resid = np.mean(np.abs(VARold[0:2]-VAR[0:2]))
    if not np.isfinite(resid) or resid > 10.0*residOld:
        VAR[...] = VARold
        phiHbyAx, phiHbyAy = phiOldX, phiOldY
        dtRelax = max(0.5*dtRelax, 1e-4)
        dt = dtRelax*min(dtCflPrev, dtMax)
        continue
    t += dt
    divU    = fvc_divFlux(phiHbyAx, phiHbyAy)
    invTime = max(np.max(np.abs(VAR[IU]))/dx, np.max(np.abs(VAR[IV]))/dy)
    history.append((it, t, dt, gx, resid))
    if it % 20 == 0 or it == nIterations-1:
        print(f'it {it:4d}  t={t:9.3f} dt={dt:7.3f}  gx={gx:.9e}  drag={gx*Afluid:.9f}  '
              f'Ub={bulkVelocity(VAR[IU]):.9f}  change={resid:.3e}  '
              f'|div|={np.mean(np.abs(divU)):.2e}  slip={np.max(np.abs(VAR[0:2][:, nearBody])):.5f}', flush=True)
    cellPe = invTime*min(dx, dy)**2/nu
    dtPe   = foMax*min(dx, dy)**2/nu/max(1.0, 0.5*cellPe)
    dtCfl  = dtMax if invTime == 0.0 else cflMax/invTime
    dtRelax = max(0.5*dtRelax, 1e-2) if resid > 1.2*residOld else min(1.5*dtRelax, 1.0)
    residOld = resid
    dtCflPrev = min(dtCfl, dtPe)
    dt = dtRelax*min(dtCfl, dtPe, dtMax)

print(f'done in {time.time()-t0:.1f}s; final gx={gx:.12e} drag={gx*Afluid:.12f} Afluid={Afluid:.12f} dt={dt}')
np.savez(args.out, nx=nx, ny=ny, Lx=Lx, Ly=Ly, R=R, nu=nu, dt=dt, bc=bcType, alphaMin=ALPHA_MIN,
         xp=xp, yp=yp, xf=xf, yf=yf, u=VAR[IU], v=VAR[IV], p=VAR[IP],
         phiX=phiHbyAx, phiY=phiHbyAy, thX=thX, thY=thY, alpha=alpha,
         Awall=Awall, dWall=dWall, Swx=Swx, Swy=Swy, gx=gx, Afluid=Afluid,
         drag=gx*Afluid, nearBody=nearBody, history=np.array(history))
