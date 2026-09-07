#!/usr/bin/env python
"""
Compare an OpenFOAM cut-cell IBM case (cylinderRe1 or cylinderRe1_pimple) with
the notebook reference (ref_nx81.npz produced by reference_notebook.py).

Usage: compare.py <case> [--ref ref_nx81.npz] [--time latest] [--fig out.png]
"""
import sys, os, argparse
import numpy as np
from foamio import read_field, latest_time, cell_ij, to_grid

ap = argparse.ArgumentParser()
ap.add_argument('case')
ap.add_argument('--ref', default=os.path.join(os.path.dirname(__file__), 'ref_nx81.npz'))
ap.add_argument('--time', default=None)
ap.add_argument('--fig', default=None)
args = ap.parse_args()

ref = np.load(args.ref, allow_pickle=True)
nx, ny = int(ref['nx']), int(ref['ny'])
Lx, Ly = float(ref['Lx']), float(ref['Ly'])
dx, dy = Lx/nx, Ly/ny
case = args.case
t = args.time or latest_time(case)
print(f'case {case}, time {t}, reference {args.ref}')

C, _ = read_field(os.path.join(case, '0', 'C'))
ix, iy = cell_ij(C, dx, dy)
assert C.shape[0] == nx*ny
G = lambda v: to_grid(v, ix, iy, nx, ny)

# mesh thickness: OpenFOAM is 3D with one cell in z
Cf, Cfp = read_field(os.path.join(case, '0', 'Cf'))
dz = 2*Cf[:, 2].max()   # face centres of x/y faces sit at z = dz/2
V3 = dx*dy*dz

def report(name, of, nb, mask=None, scale=None):
    d = of - nb
    if mask is not None:
        d = d[mask]; nbm = nb[mask]
    else:
        nbm = nb
    scale = scale or max(np.max(np.abs(nbm)), 1e-300)
    print(f'  {name:28s} max|diff| = {np.max(np.abs(d)):.3e}   '
          f'rms = {np.sqrt(np.mean(d**2)):.3e}   (ref max |.| = {scale:.3e}, '
          f'rel max = {np.max(np.abs(d))/scale:.2e})')
    return d

# ---------------------------------------------------------------- geometry
print('\nGeometry (from ibmSetGeometry, time 0):')
alpha_of = G(read_field(os.path.join(case, '0', 'ibmAlpha'))[0])
report('alpha (cell fluid fraction)', alpha_of, ref['alpha'])

Sw_of, _ = read_field(os.path.join(case, '0', 'ibmSw'))
# the notebook stores +sum theta n A (pointing out of the body), OpenFOAM
# Sw = -sum theta Sf (pointing into the body)
report('Sw_x (wall area vector)', -G(Sw_of[:, 0])/dz, ref['Swx'])
report('Sw_y', -G(Sw_of[:, 1])/dz, ref['Swy'])

dW_of = G(read_field(os.path.join(case, '0', 'ibmWallDistance'))[0])
wall = ref['Awall'] > 0
report('wall distance (wall cells)', dW_of, ref['dWall'], mask=wall)

# theta face by face: x-faces have Cf.x on the face grid and Cf.y at cell centres
th_of, th_p = read_field(os.path.join(case, '0', 'ibmTheta'))
thX = np.full((ny, nx+1), np.nan); thY = np.full((ny+1, nx), np.nan)
def place(cf, th):
    if np.isscalar(th):
        th = np.full(len(cf), th)
    for c, v in zip(cf, th):
        fx = c[0]/dx; fy = c[1]/dy
        if abs(fx - round(fx)) < 1e-6 and abs(fy - round(fy) ) > 0.25:   # x-face
            thX[int(round(fy - 0.5)), int(round(fx))] = v
        elif abs(fy - round(fy)) < 1e-6 and abs(fx - round(fx)) > 0.25:  # y-face
            thY[int(round(fy)), int(round(fx - 0.5))] = v
place(Cf, th_of)
for pname in th_p:
    if pname in Cfp and len(Cfp[pname]):
        place(Cfp[pname], th_p[pname])
print(f'  x-faces matched: {np.isfinite(thX).sum()}/{thX.size}, y-faces: {np.isfinite(thY).sum()}/{thY.size}')
report('theta on x-faces', np.nan_to_num(thX, nan=1.0), ref['thX'])
report('theta on y-faces', np.nan_to_num(thY, nan=1.0), ref['thY'])
print(f'  cut faces (0<theta<1): OpenFOAM {((thX>0)&(thX<1)).sum()+((thY>0)&(thY<1)).sum()}, '
      f'notebook {((ref["thX"]>0)&(ref["thX"]<1)).sum()+((ref["thY"]>0)&(ref["thY"]<1)).sum()}')
Af_of = np.nansum(alpha_of)*dx*dy
print(f'  fluid area: OpenFOAM {Af_of:.9f}, notebook {float(ref["Afluid"]):.9f}')

# ---------------------------------------------------------------- solution
print(f'\nSolution at time {t}:')
U_of, _ = read_field(os.path.join(case, t, 'U'))
p_of, _ = read_field(os.path.join(case, t, 'p'))
u = G(U_of[:, 0]); v = G(U_of[:, 1]); p = G(p_of)
fluid = ref['alpha'] > 0
alphaV = ref['alpha']*dx*dy
Afl = alphaV.sum()
# pressure level: compare after removing the fluid-volume mean
def demean(f):
    return f - np.sum(f*alphaV)/Afl
report('u', u, ref['u'])
report('v', v, ref['v'])
report('p (fluid, mean removed)', demean(p), demean(ref['p']), mask=fluid)
ub_of = np.sum(u*alphaV)/Afl
print(f'  bulk velocity: OpenFOAM {ub_of:.9f}, notebook {np.sum(ref["u"]*alphaV)/Afl:.9f}')
near = ref['nearBody']
print(f'  max |U| near body: OpenFOAM {np.max(np.hypot(u, v)[near]):.6f}, notebook {np.max(np.hypot(ref["u"], ref["v"])[near]):.6f}')
print(f'  max |U|: OpenFOAM {np.max(np.hypot(u, v)):.6f}, notebook {np.max(np.hypot(ref["u"], ref["v"])):.6f}')

# forcing and drag
gx_nb = float(ref['gx']); drag_nb = float(ref['drag'])
props = None
for root, dirs, files in os.walk(os.path.join(case, t)):
    for f in files:
        if f.endswith('Properties') and root.endswith('uniform'):
            props = os.path.join(root, f)
gx_of = None
import re, glob
if props:
    m = re.search(r'gradient\s+([-+0-9.eE]+)', open(props).read())
    gx_of = float(m.group(1))
else:
    logs = glob.glob(os.path.join(case, 'log.ibm*Foam'))
    if logs:
        vals = re.findall(r'pressure gradient = ([-+0-9.eE]+)', open(logs[0]).read())
        if vals:
            gx_of = float(vals[-1])
if gx_of is not None:
    print(f'  driving force gx: OpenFOAM {gx_of:.9e}, notebook {gx_nb:.9e}, rel diff {(gx_of-gx_nb)/gx_nb:.2e}')
    print(f'  drag gx*Afluid : OpenFOAM {gx_of*Af_of:.9f}, notebook {drag_nb:.9f}')
fdat = os.path.join(case, 'postProcessing', 'ibmForces')
if os.path.isdir(fdat):
    files = [os.path.join(r, f) for r, d, fs in os.walk(fdat) for f in fs]
    last = np.loadtxt(sorted(files)[-1])[-1]
    Fp, Fv = last[1:4]/dz, last[4:7]/dz
    print(f'  drag from wall terms (per unit depth): pressure {Fp[0]:.6f} + viscous {Fv[0]:.6f} = {Fp[0]+Fv[0]:.6f} '
          f'(notebook drag {drag_nb:.6f}, rel diff {(Fp[0]+Fv[0]-drag_nb)/drag_nb:.2e}); lift {Fp[1]+Fv[1]:.2e}')

if args.fig:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    xp, yp = ref['xp'], ref['yp']
    fig, ax = plt.subplots(2, 3, figsize=(15, 9))
    th = np.linspace(0, 2*np.pi, 200)
    xb, yb = Lx/2 + np.cos(th), Ly/2 + np.sin(th)
    for a, f, ttl in [(ax[0,0], u, 'u OpenFOAM'), (ax[0,1], ref['u'], 'u notebook'),
                      (ax[0,2], u-ref['u'], 'u difference')]:
        im = a.pcolormesh(xp, yp, f, shading='nearest'); plt.colorbar(im, ax=a); a.plot(xb, yb, 'w', lw=0.8)
        a.set_title(ttl); a.set_aspect('equal')
    for a, f, ttl in [(ax[1,0], demean(p), 'p OpenFOAM'), (ax[1,1], demean(ref['p']), 'p notebook')]:
        im = a.pcolormesh(xp, yp, np.where(fluid, f, np.nan), shading='nearest'); plt.colorbar(im, ax=a); a.plot(xb, yb, 'w', lw=0.8)
        a.set_title(ttl); a.set_aspect('equal')
    j = ny//2
    a = ax[1,2]
    a.plot(yp, u[:, nx//2], 'o', ms=3, label='u(x=5) OpenFOAM')
    a.plot(yp, ref['u'][:, nx//2], '-', label='u(x=5) notebook')
    a.plot(xp, u[j, :], 's', ms=3, label='u(y=5) OpenFOAM')
    a.plot(xp, ref['u'][j, :], '-', label='u(y=5) notebook')
    a.legend(); a.set_xlabel('x or y'); a.set_title('profiles through the cylinder centre')
    fig.suptitle(f'{os.path.basename(os.path.abspath(case))} vs notebook, time {t}')
    fig.tight_layout()
    fig.savefig(args.fig, dpi=110)
    print('figure written to', args.fig)
