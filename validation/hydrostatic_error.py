#!/usr/bin/env python
"""Spurious velocity and pressure error of a tutorials/order/hydrostatic run: fluid
at rest under g = (0, -1, 0), exact U = 0, p = -y + const (evaluated at the fluid
centroids), max over all fluid cells and over the wall cells."""
import sys, os, numpy as np
from foamio import read_field, latest_time
case = sys.argv[1]
t = latest_time(case)
C, _ = read_field(os.path.join(case, '0', 'C')); U, _ = read_field(os.path.join(case, t, 'U'))
p, _ = read_field(os.path.join(case, t, 'p'))
geo = t if os.path.exists(os.path.join(case, t, 'ibmAlpha')) else '0'   # static body: written at 0
a, _ = read_field(os.path.join(case, geo, 'ibmAlpha'))
offf = os.path.join(case, t, 'ibmCentroidOffset')
off = read_field(offf)[0] if os.path.exists(offf) else np.zeros_like(C)
y = C[:, 1] + off[:, 1]
fl = a > 0; cut = fl & (a < 1)
ep = p + y; ep -= np.mean(ep[fl & (a == 1)])
mU = np.linalg.norm(U, axis=1)
n = int(round(np.sqrt(len(C))))
i = np.argmax(np.where(fl, mU, 0))
print(f'{os.path.basename(os.path.normpath(case))}: N={n} t={t} max|U| {mU[fl].max():.3e} (cut cells {mU[cut].max():.3e}), '
      f'rms|U| {np.sqrt(np.mean(mU[fl]**2)):.3e}; p error max {np.abs(ep[fl]).max():.3e} (cut cells {np.abs(ep[cut]).max():.3e}); max|U| at alpha={a[i]:.3f} ({C[i,0]:.3f}, {C[i,1]:.3f})')
