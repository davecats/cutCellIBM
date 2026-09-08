#!/usr/bin/env python
"""
Reference solution for moving bodies: the cut-cell IBM of
.asset/fvm/SIMPLEC_faceArea_moving.ipynb run as a script. The numerics are
the notebook's cells 1-6 verbatim (plus the free-stream test of cell 7);
only the plotting is removed and the results are saved for comparison with
ibmPimpleFoam.

Usage:
  reference_moving.py --motion oscillate --nSteps 400 --out ref_osc.npz
  reference_moving.py --motion translate --uTrans 0.5 --uBulk 0 --nSteps 1200 --out ref_gal.npz
  reference_moving.py --freeStream          (prints the free-stream table)
"""
import os, sys, argparse, time
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu

ap = argparse.ArgumentParser()
ap.add_argument('--nx', type=int, default=81)
ap.add_argument('--dt', type=float, default=0.025)
ap.add_argument('--nOuter', type=int, default=3)
ap.add_argument('--nSteps', type=int, default=400)
ap.add_argument('--bc', default='noSlip')
ap.add_argument('--motion', default='oscillate')
ap.add_argument('--uTrans', type=float, default=0.5)
ap.add_argument('--aOsc', type=float, default=1.0)
ap.add_argument('--fOsc', type=float, default=0.05)
ap.add_argument('--uBulk', type=float, default=0.0)
ap.add_argument('--alphaMin', type=float, default=0.02)
ap.add_argument('--freeStream', action='store_true')
ap.add_argument('--seamFix', action='store_true', help='close both copies of a periodic face when an absorbed cell touches it (see README)')
ap.add_argument('--out', default='ref_moving.npz')
args = ap.parse_args()

# Inputs (cell 1)
R  = 1
Lx = 10
Ly = 10
nx = args.nx
ny = args.nx
dt     = args.dt
nOuter = args.nOuter
nu = 1.0
uBulkTarget = args.uBulk
bcType = args.bc
motion   = args.motion
uTrans   = args.uTrans
aOsc     = args.aOsc
fOsc     = args.fOsc
xc0, yc0 = Lx/2, Ly/2

def bodyPosition(t):
    if motion == 'still':      return xc0, yc0
    if motion == 'translate':  return xc0+uTrans*t, yc0
    if motion == 'oscillate':  return xc0+aOsc*np.sin(2*np.pi*fOsc*t), yc0
    raise ValueError(f'unknown motion {motion!r}')

def bodyVelocity(t):
    if motion == 'still':      return np.array([0.0, 0.0])
    if motion == 'translate':  return np.array([uTrans, 0.0])
    if motion == 'oscillate':
        return np.array([2*np.pi*fOsc*aOsc*np.cos(2*np.pi*fOsc*t), 0.0])
    raise ValueError(f'unknown motion {motion!r}')

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
def fvm_ddt_LHS(a, alphaVNew, dt):
    a[0] += alphaVNew/dt
def fvm_ddt_RHS(rhs, phi, alphaVOld, dt):
    rhs += alphaVOld*phi/dt
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

# Cut cells on a moving body (cell 5)
def levelSet(x, y, xc, yc):
    dxb = (x-xc+Lx/2) % Lx - Lx/2
    dyb = (y-yc+Ly/2) % Ly - Ly/2
    return np.hypot(dxb, dyb) - R

ALPHA_MIN = args.alphaMin
IBM_SOLID_WEIGHT = 1e12

def segmentFluidFraction(x0, y0, x1, y1, xc, yc, nBisect=40):
    f0 = levelSet(x0, y0, xc, yc)
    f1 = levelSet(x1, y1, xc, yc)
    frac = (f0 > 0.0).astype(float)
    cut  = (f0 > 0.0) != (f1 > 0.0)
    if cut.any():
        xa, ya = x0[cut], y0[cut]
        sx, sy = (x1-x0)[cut], (y1-y0)[cut]
        fa = f0[cut]
        lo = np.zeros_like(fa)
        hi = np.ones_like(fa)
        for _ in range(nBisect):
            mid = 0.5*(lo+hi)
            sameAsA = (levelSet(xa+mid*sx, ya+mid*sy, xc, yc) > 0.0) == (fa > 0.0)
            lo = np.where(sameAsA, mid, lo)
            hi = np.where(sameAsA, hi, mid)
        tCross = 0.5*(lo+hi)
        frac[cut] = np.clip(np.where(fa > 0.0, tCross, 1.0-tCross), 0.0, 1.0)
    return frac

def faceFractions(xc, yc):
    XF1 = np.tile(xf, (ny, 1))
    Y0  = np.tile(yf[:-1][:, None], (1, nx+1))
    thX = segmentFluidFraction(XF1, Y0, XF1, Y0+dy, xc, yc)
    YF1 = np.tile(yf[:, None], (1, nx))
    X0  = np.tile(xf[:-1], (ny+1, 1))
    thY = segmentFluidFraction(X0, YF1, X0+dx, YF1, xc, yc)
    thX[:, -1] = thX[:, 0]
    thY[-1, :] = thY[0, :]
    return thX, thY

def cellFraction(xc, yc, nSub=32):
    sd  = levelSet(XP, YP, xc, yc)
    cut = np.abs(sd) < 0.5*np.hypot(dx, dy)
    a   = (sd > 0.0).astype(float)
    if cut.any():
        x0, y0 = XF[:-1, :-1][cut], YF[:-1, :-1][cut]
        s   = (np.arange(nSub)+0.5)/nSub
        acc = np.zeros(x0.shape)
        for sx in s:
            for sy in s:
                acc += levelSet(x0+sx*dx, y0+sy*dy, xc, yc) > 0.0
        a[cut] = acc/(nSub*nSub)
    return a

def updateGeometry(xc, yc):
    global thX, thY, alpha, alphaV, fluid, solid, nearBody, sdCell
    global Swx, Swy, Awall, dWall
    thX, thY = faceFractions(xc, yc)
    alpha    = cellFraction(xc, yc)
    sdCell   = levelSet(XP, YP, xc, yc)
    tiny  = alpha < ALPHA_MIN
    alpha = np.where(tiny, 0.0, alpha)
    thX[:, :-1] = np.where(tiny, 0.0, thX[:, :-1]); thX[:, 1:] = np.where(tiny, 0.0, thX[:, 1:])
    thY[:-1, :] = np.where(tiny, 0.0, thY[:-1, :]); thY[1:, :] = np.where(tiny, 0.0, thY[1:, :])
    if args.seamFix:
        # x = 0 and x = Lx are the same face: an absorbed cell on one side must
        # close it on both sides (the notebook copies theta across the seam
        # before absorbing slivers, so the two copies can disagree there)
        thX[:, 0] = thX[:, -1] = np.minimum(thX[:, 0], thX[:, -1])
        thY[0, :] = thY[-1, :] = np.minimum(thY[0, :], thY[-1, :])
    fluid    = alpha > 0.0
    solid    = ~fluid
    alphaV   = alpha*V
    nearBody = fluid & (np.roll(solid, 1, 1) | np.roll(solid, -1, 1) |
                        np.roll(solid, 1, 0) | np.roll(solid, -1, 0))
    Swx   = (thX[:, 1:]-thX[:, :-1])*dy
    Swy   = (thY[1:, :]-thY[:-1, :])*dx
    Awall = np.hypot(Swx, Swy)
    dWall = np.where(Awall > 0.0, alphaV/(2.0*np.maximum(Awall, 1e-300)), np.inf)

def ibm_noSlipWall(a, Gamma):
    a[0] -= np.where(Awall > 0.0, Gamma*Awall/dWall, 0.0)
def ibm_wallSource(Gamma, uWallComponent):
    return np.where(Awall > 0.0, Gamma*Awall/dWall, 0.0)*uWallComponent
def ibm_blankSolid(a):
    a[0] += np.where(solid, IBM_SOLID_WEIGHT, 0.0)
def ibm_blankSource(uWallComponent):
    return np.where(solid, IBM_SOLID_WEIGHT, 0.0)*uWallComponent
def ibm_movingWall(uWall, dt):
    S = -(uWall[0]*Swx + uWall[1]*Swy)
    return S, alphaV - dt*S

# Run controls (cell 6)
consistent = True

def bulkVelocity(u):
    return np.sum(u*alphaV)/np.sum(alphaV)

def referenceCell():
    return divmod(int(np.argmax(np.where(fluid, sdCell, -np.inf))), nx)

def bodyForce(p, U, uWall):
    fPres = -np.array([np.sum(p*Swx), np.sum(p*Swy)])
    if bcType == 'noSlip':
        c = np.where(Awall > 0.0, nu*Awall/dWall, 0.0)
        fVisc = np.array([np.sum(c*(U[iC]-uWall[iC])) for iC in range(2)])
    else:
        fVisc = np.zeros(2)
    return fPres+fVisc

def timeStep(t, dt, VAR, phiHbyAx, phiHbyAy):
    VARold = VAR.copy()
    updateGeometry(*bodyPosition(t+dt))
    uWall = bodyVelocity(t+dt)
    massSrc, alphaVOld = ibm_movingWall(uWall, dt)
    jRef, iRef = referenceCell()
    phiHbyAx, phiHbyAy = interpolatePhi(VARold[IU], VARold[IV])
    rhs = np.zeros((2, ny, nx))
    for iC in range(2):
        fvm_ddt_RHS(rhs[iC], VARold[iC], alphaVOld, dt)
        rhs[iC] += ibm_blankSource(uWall[iC])
        if bcType == 'noSlip':
            rhs[iC] += ibm_wallSource(nu, uWall[iC])
    unitSource = np.stack([alphaV, np.zeros((ny, nx))])
    noField    = np.zeros((ny, nx))
    lapl = zeroStencil()
    fvm_lap(lapl, nu*np.ones((ny, nx)))
    if bcType == 'noSlip':
        ibm_noSlipWall(lapl, nu)
    elif bcType != 'freeSlip':
        raise ValueError("bcType must be 'noSlip' or 'freeSlip'")
    for outer in range(nOuter):
        aMom = []
        for iC in range(2):
            a = zeroStencil()
            fvm_ddt_LHS(a, alphaV, dt)
            a -= lapl
            aMom.append(a)
        for a in aMom:
            fvm_divConvFlux(a, phiHbyAx, phiHbyAy)
        for a in aMom:
            ibm_blankSolid(a)
        aP     = np.stack([a[0] for a in aMom])
        sumOff = np.stack([a[1:].sum(axis=0) for a in aMom])
        aP    = np.maximum(aP, alphaV/dt)
        rAU   = 1.0/aP
        rAtU  = 1.0/np.maximum(aP+sumOff, alphaV/dt) if consistent else rAU
        rAtUV = rAtU*V
        pStencil = zeroStencil()
        fvm_lap(pStencil, rAtUV[IU], rAtUV[IV])
        pStencil[0][solid] = -1.0
        pStencil[0, jRef, iRef] *= 2.0
        Ap = toCSR(pStencil)
        luMom = [splu(toCSR(a).tocsc()) for a in aMom]
        luP   = splu(Ap.tocsc())

        def advance(sourceField, pOld, massSource):
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
            rhsP = div0 + massSource
            rhsP[solid] = 0.0
            pNew = luP.solve(rhsP.ravel()).reshape(ny, nx)
            fX, fY = fvc_grdFlux(pNew, rAtUV[IU], rAtUV[IV])
            return (HbyA - rAtU*fvc_grd(pNew), pNew, pX-fX, pY-fY, div0)

        base = advance(rhs,        VAR[IP], massSrc)
        resp = advance(unitSource, noField, noField)
        gx   = (uBulkTarget-bulkVelocity(base[0][IU]))/bulkVelocity(resp[0][IU])
        VARprev  = VAR.copy()
        VAR[0:2] = base[0] + gx*resp[0]
        VAR[IP]  = base[1] + gx*resp[1]
        phiHbyAx = base[2] + gx*resp[2]
        phiHbyAy = base[3] + gx*resp[3]
        outerResid = np.mean(np.abs(VARprev[0:2]-VAR[0:2]))
        meanDiv0   = np.mean(np.abs(base[4] + gx*resp[4]))
    contErr = fvc_divFlux(phiHbyAx, phiHbyAy) + massSrc
    slip    = (np.max(np.abs(VAR[0:2][:, nearBody] - uWall[:, None]))
               if nearBody.any() else 0.0)
    diag = dict(t=t+dt, body=bodyPosition(t+dt), uWall=uWall, gx=gx,
                uBulk=bulkVelocity(VAR[IU]), force=bodyForce(VAR[IP], VAR[0:2], uWall),
                contMax=np.max(np.abs(contErr)), contSum=contErr.sum(),
                meanDiv0=meanDiv0, outerResid=outerResid,
                alphaOldMin=alphaVOld[fluid].min()/V,
                nSwept=int(((alphaVOld < 0.0) & fluid).sum()), slip=slip,
                cfl=dt*max(np.max(np.abs(VAR[IU]))/dx, np.max(np.abs(VAR[IV]))/dy))
    return phiHbyAx, phiHbyAy, diag


if args.freeStream:
    def freeStreamTest(uBody, nStep=5, nOuterTest=None):
        global motion, uTrans, uBulkTarget, nOuter
        keep = (motion, uTrans, uBulkTarget, nOuter, VAR.copy())
        motion, uTrans, uBulkTarget = 'translate', uBody, uBody
        if nOuterTest is not None:
            nOuter = nOuterTest
        tT = 0.0
        updateGeometry(*bodyPosition(tT))
        VAR[IU], VAR[IV], VAR[IP] = uBody, 0.0, 0.0
        pX, pY = interpolatePhi(VAR[IU], VAR[IV])
        err = 0.0
        for _ in range(nStep):
            pX, pY, dg = timeStep(tT, dt, VAR, pX, pY)
            tT  = dg['t']
            err = max(err, np.max(np.abs(VAR[IU][fluid]-uBody)),
                           np.max(np.abs(VAR[IV][fluid])))
        motion, uTrans, uBulkTarget, nOuter = keep[0], keep[1], keep[2], keep[3]
        VAR[...] = keep[4]
        return err
    print(f'T1  free-stream preservation, bcType = {bcType}: max |u - u_body| after 5 steps, dt = {dt}')
    print( '    nOuter     u_body = 0.37     u_body = 2.00')
    for k in (1, 2, 3, 5):
        print(f'    {k:^6d}     {freeStreamTest(0.37, nOuterTest=k):.3e}         '
              f'{freeStreamTest(2.00, nOuterTest=k):.3e}', flush=True)
    sys.exit(0)

# Time loop
t = 0.0
updateGeometry(*bodyPosition(t))
VAR[...] = 0.0
phiHbyAx, phiHbyAy = interpolatePhi(VAR[IU], VAR[IV])
history = []
t0 = time.time()
for it in range(args.nSteps):
    phiHbyAx, phiHbyAy, d = timeStep(t, dt, VAR, phiHbyAx, phiHbyAy)
    t = d['t']
    history.append([d['t'], d['body'][0], d['uWall'][0], d['gx'], d['uBulk'],
                    d['force'][0], d['force'][1], d['contMax'], d['contSum'],
                    d['alphaOldMin'], d['nSwept'], d['slip'], d['outerResid']])
    if it % 20 == 0 or it == args.nSteps-1:
        print(f'step {it:4d} t={t:8.4f} xb={d["body"][0]:.4f} ub={d["uWall"][0]:+.4f} '
              f'gx={d["gx"]:+.6e} Ub={d["uBulk"]:+.2e} F=({d["force"][0]:+.5f}, {d["force"][1]:+.5f}) '
              f'cont={d["contMax"]:.1e}/{d["contSum"]:.1e} aOldMin={d["alphaOldMin"]:+.4f}/{d["nSwept"]} '
              f'slip={d["slip"]:.4f} outer={d["outerResid"]:.2e}', flush=True)
print(f'done in {time.time()-t0:.1f}s')
np.savez(args.out, nx=nx, Lx=Lx, Ly=Ly, R=R, nu=nu, dt=dt, nOuter=nOuter, bc=bcType,
         motion=motion, xp=xp, yp=yp, xf=xf, yf=yf, t=t, body=bodyPosition(t),
         uWall=bodyVelocity(t), u=VAR[IU], v=VAR[IV], p=VAR[IP],
         phiX=phiHbyAx, phiY=phiHbyAy, thX=thX, thY=thY, alpha=alpha,
         Awall=Awall, dWall=dWall, Swx=Swx, Swy=Swy, nearBody=nearBody,
         history=np.array(history),
         historyColumns='t xb ub gx uBulk Fx Fy contMax contSum alphaOldMin nSwept slip outerResid')
