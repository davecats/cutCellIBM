#!/usr/bin/env python
"""max |U - u_b| over the fluid cells of a free-stream case at its latest time."""
import sys, os, numpy as np
from foamio import read_field, latest_time
case = sys.argv[1]; ub = float(sys.argv[2])
t = latest_time(case)
U, _ = read_field(os.path.join(case, t, 'U'))
a, _ = read_field(os.path.join(case, t, 'ibmAlpha'))
fluid = a > 0
print(f'{case}: t = {t}, max |u - u_b| = {np.max(np.abs(U[fluid, 0] - ub)):.3e}, max |v| = {np.max(np.abs(U[fluid, 1])):.3e}')
