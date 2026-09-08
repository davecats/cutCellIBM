#!/usr/bin/env python
"""
Compare an ibmPimpleFoam moving-body case with the notebook reference produced
by reference_moving.py: force history, driving force, and the final fields.

Usage: compare_moving.py <case> --ref ref_osc.npz [--fig out.png] [--tail 400]
"""
import os, re, glob, argparse
import numpy as np
from foamio import read_field, latest_time, cell_ij, to_grid

ap = argparse.ArgumentParser()
ap.add_argument('case')
ap.add_argument('--ref', required=True)
ap.add_argument('--fig', default=None)
ap.add_argument('--tail', type=int, default=0, help='average the force over the last N steps')
ap.add_argument('--time', default=None, help='time directory to compare the fields at (default: latest)')
args = ap.parse_args()

ref = np.load(args.ref, allow_pickle=True)
nx = int(ref['nx']); ny = nx
Lx, Ly = float(ref['Lx']), float(ref['Ly'])
dx, dy = Lx/nx, Ly/ny
case = args.case
h = ref['history']          # t xb ub gx uBulk Fx Fy contMax contSum alphaOldMin nSwept slip outerResid

# ---- force history from postProcessing/ibmForces/*/forces.dat
files = sorted(glob.glob(os.path.join(case, 'postProcessing', 'ibmForces', '*', 'forces.dat')))
dat = np.loadtxt(files[-1])
Cf, _ = read_field(os.path.join(case, '0', 'Cf'))
dz = 2*Cf[:, 2].max()
t_of = dat[:, 0]; Fx_of = (dat[:, 1] + dat[:, 4])/dz; Fy_of = (dat[:, 2] + dat[:, 5])/dz
Fp_of = dat[:, 1]/dz; Fv_of = dat[:, 4]/dz
t_nb = h[:, 0]; Fx_nb = h[:, 5]
n = min(len(t_of), len(t_nb))
assert np.allclose(t_of[:n], t_nb[:n], atol=1e-6), 'time levels differ'
d = Fx_of[:n] - Fx_nb[:n]
print(f'case {case}: {n} common steps, t = {t_of[0]:.4f} .. {t_of[n-1]:.4f}')
print(f'  F_x history: max|diff| = {np.max(np.abs(d)):.3e}, rms = {np.sqrt(np.mean(d**2)):.3e}, '
      f'max|F_x| = {np.max(np.abs(Fx_nb[:n])):.3f}')
print(f'  F_x history, excluding the impulsive first 10 steps: max|diff| = {np.max(np.abs(d[10:])):.3e}, '
      f'rel to max|F| = {np.max(np.abs(d[10:]))/np.max(np.abs(Fx_nb[10:n])):.2e}')
for i in [0, 1, 2, 10, 50, 100, 200, 300, n-1]:
    if i < n:
        print(f'    t={t_of[i]:8.4f}  F_x OpenFOAM {Fx_of[i]:+10.5f} (p {Fp_of[i]:+9.5f}, v {Fv_of[i]:+9.5f})   notebook {Fx_nb[i]:+10.5f}')
if args.tail:
    m = args.tail
    print(f'  mean F_x over the last {m} steps: OpenFOAM {Fx_of[n-m:n].mean():+.5f} +- {Fx_of[n-m:n].std():.5f}, '
          f'notebook {Fx_nb[n-m:n].mean():+.5f} +- {Fx_nb[n-m:n].std():.5f}')

# ---- driving force from the log
logs = glob.glob(os.path.join(case, 'log.ibm*Foam'))
if logs:
    gx = [float(x) for x in re.findall(r'pressure gradient = ([-+0-9.eE]+)', open(logs[0]).read())]
    print(f'  driving force gx at the end: OpenFOAM {gx[-1]:+.6e}, notebook {h[n-1, 3]:+.6e}')
    ub = [float(x) for x in re.findall(r'uncorrected Ubar = ([-+0-9.eE]+)', open(logs[0]).read())]
    print(f'  residual bulk velocity before the last correction: OpenFOAM {ub[-1]:+.2e} (notebook pins it exactly)')

# ---- final fields
t = args.time or latest_time(case)
if abs(float(t) - float(ref['t'])) < 1e-6:
    C, _ = read_field(os.path.join(case, '0', 'C'))
    ix, iy = cell_ij(C, dx, dy)
    G = lambda v: to_grid(v, ix, iy, nx, ny)
    U_of, _ = read_field(os.path.join(case, t, 'U'))
    p_of, _ = read_field(os.path.join(case, t, 'p'))
    a_of = G(read_field(os.path.join(case, t, 'ibmAlpha'))[0])
    u = G(U_of[:, 0]); v = G(U_of[:, 1]); p = G(p_of)
    fluid = ref['alpha'] > 0
    alphaV = ref['alpha']*dx*dy
    dm = lambda f: f - np.sum(f*alphaV)/alphaV.sum()
    print(f'\nFinal fields at t = {t} (body at x = {float(ref["body"][0]):.4f}):')
    print(f'  alpha        max|diff| = {np.max(np.abs(a_of - ref["alpha"])):.3e}')
    for name, of, nb in [('u', u, ref['u']), ('v', v, ref['v'])]:
        e = np.abs(of - nb)[fluid]
        E = np.where(fluid, np.abs(of - nb), 0)
        j, i = np.unravel_index(np.argmax(E), E.shape)
        print(f'  {name:12s} max|diff| = {e.max():.3e}  rms = {np.sqrt((e**2).mean()):.3e}  (max|{name}| = {np.max(np.abs(nb)):.3e})'
              f'  at cell ({i},{j}) x={(i+0.5)*dx:.3f} y={(j+0.5)*dy:.3f}: OF {of[j,i]:+.5f} nb {nb[j,i]:+.5f}, alpha nb {ref["alpha"][j,i]:.4f} OF {a_of[j,i]:.4f}')
    e = np.abs(dm(p) - dm(ref['p']))[fluid]
    print(f'  p (mean rm.) max|diff| = {e.max():.3e}  rms = {np.sqrt((e**2).mean()):.3e}  (max|p| = {np.max(np.abs(dm(ref["p"])[fluid])):.3e})')
    near = ref['nearBody']
    uW = ref['uWall']
    print(f'  max |U - u_b| near the body (component-wise): OpenFOAM {max(np.max(np.abs(u - uW[0])[near]), np.max(np.abs(v - uW[1])[near])):.5f}, notebook {float(h[n-1, 11]):.5f}')
else:
    print(f'\n(final time {t} differs from the notebook {float(ref["t"])}: fields not compared)')

if args.fig:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    ax[0].plot(t_nb[:n], Fx_nb[:n], '-', label='notebook')
    ax[0].plot(t_of[:n], Fx_of[:n], '--', label='ibmPimpleFoam')
    ax[0].plot(t_nb[:n], h[:n, 2]*10, ':', color='0.5', label='10 x body velocity')
    ax[0].set_ylabel('F_x on the body'); ax[0].legend(); ax[0].set_ylim(-12, 12)
    ax[1].semilogy(t_of[:n], np.abs(d) + 1e-18, label='|F_x difference|')
    ax[1].set_xlabel('t'); ax[1].legend()
    fig.suptitle(os.path.basename(os.path.abspath(case)) + ' vs notebook')
    fig.tight_layout(); fig.savefig(args.fig, dpi=110); print('figure written to', args.fig)
