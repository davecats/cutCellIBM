#!/usr/bin/env python
"""Plot the scalar field of a tutorials/scalar case at its latest time."""
import sys, os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from foamio import read_field, latest_time, cell_ij, to_grid
cases = sys.argv[1:-1]; out = sys.argv[-1]
fig, axs = plt.subplots(1, len(cases), figsize=(5.2*len(cases), 4.6))
for ax, c in zip(np.atleast_1d(axs), cases):
    t = latest_time(c); nx = 81; dx = 10/81
    C, _ = read_field(os.path.join(c, '0', 'C')) if os.path.exists(os.path.join(c, '0', 'C')) else (None, None)
    T, _ = read_field(os.path.join(c, t, 'T')); a, _ = read_field(os.path.join(c, t, 'ibmAlpha'))
    ix = np.arange(nx*nx) % nx; iy = np.arange(nx*nx)//nx           # blockMesh ordering
    G = lambda v: to_grid(v, ix, iy, nx, nx)
    Tg = np.where(G(a) > 0, G(T), np.nan)
    im = ax.pcolormesh((np.arange(nx)+.5)*dx, (np.arange(nx)+.5)*dx, Tg, shading='nearest')
    ax.contour((np.arange(nx)+.5)*dx, (np.arange(nx)+.5)*dx, G(a), [0.5], colors='w', linewidths=0.8)
    plt.colorbar(im, ax=ax); ax.set_aspect('equal'); ax.set_title(f'{os.path.basename(os.path.normpath(c))}, t = {t}')
fig.tight_layout(); fig.savefig(out, dpi=110); print('written', out)
