#!/usr/bin/env python
"""Error of a taylorCouette run against the exact solution u_theta = A r + B/r,
p = A^2 r^2/2 + 2AB ln r - B^2/(2 r^2), evaluated at the fluid centroids."""
import sys, os, numpy as np
from foamio import read_field, latest_time
case = sys.argv[1]; Ri, Ro, Om, cx, cy = 0.5, 1.5, 1.0, 1.6, 1.6
A = -Om*Ri**2/(Ro**2 - Ri**2); B = Om*Ri**2*Ro**2/(Ro**2 - Ri**2)
t = latest_time(case)
C, _ = read_field(os.path.join(case, '0', 'C')); U, _ = read_field(os.path.join(case, t, 'U')); p, _ = read_field(os.path.join(case, t, 'p'))
a, _ = read_field(os.path.join(case, t, 'ibmAlpha'))
offf = os.path.join(case, t, 'ibmCentroidOffset')
off = read_field(offf)[0] if os.path.exists(offf) else np.zeros_like(C)
x = C[:, 0] + off[:, 0] - cx; y = C[:, 1] + off[:, 1] - cy; r = np.hypot(x, y)
fl = a > 0
ut = A*r + B/r
Uex = np.stack([-ut*y/r, ut*x/r], 1)
pex = A*A*r*r/2 + 2*A*B*np.log(np.maximum(r, 1e-300)) - B*B/(2*r*r)
w = a*C[:, 2]*0 + a   # alpha weights (uniform cells)
eU = np.hypot(U[:, 0] - Uex[:, 0], U[:, 1] - Uex[:, 1])
ep = (p - pex); ep = ep - np.sum(ep[fl]*w[fl])/np.sum(w[fl])
def norms(m):
    return (np.sqrt(np.sum(eU[m]**2*w[m])/np.sum(w[m])), eU[m].max(),
            np.sqrt(np.sum(ep[m]**2*w[m])/np.sum(w[m])), np.abs(ep[m]).max())
n = int(round(np.sqrt(len(C))))
a1 = norms(fl); a2 = norms(fl & (a >= 0.5))
print(f'{os.path.basename(os.path.normpath(case))}: N={n} t={t} all fluid cells: L2(U)={a1[0]:.3e} Linf(U)={a1[1]:.3e} L2(p)={a1[2]:.3e} Linf(p)={a1[3]:.3e} | '
      f'cells with alpha>=0.5: L2(U)={a2[0]:.3e} Linf(U)={a2[1]:.3e} L2(p)={a2[2]:.3e} Linf(p)={a2[3]:.3e}   (p range {pex[fl].max()-pex[fl].min():.3f})')
