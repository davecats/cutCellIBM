"""Taylor-Couette pressure error (mean removed over the fluid, at the value points) and its radial
profile in r bins of 0.1. Usage: taylorCouette_pressureError.py case... [--outer: outer cylinder rotating]."""
import sys, os
import numpy as np
sys.path.insert(0, "/home/ws/xt8786/Codes/openFoam.ibm/validation")
from foamio import read_field, latest_time
Ri, Ro, Om, cx, cy = 0.5, 1.5, 1.0, 1.6, 1.6
args = sys.argv[1:]
if "--outer" in args:
    args.remove("--outer"); A = Om*Ro**2/(Ro**2 - Ri**2); B = -Om*Ri**2*Ro**2/(Ro**2 - Ri**2)
else:
    A = -Om*Ri**2/(Ro**2 - Ri**2); B = Om*Ri**2*Ro**2/(Ro**2 - Ri**2)
pex = lambda r: A*A*r*r/2 + 2*A*B*np.log(r) - B*B/(2*r*r)
bins = np.linspace(Ri, Ro, 11)
print("%-26s %9s" % ("case", "L2") + "".join("%9.2f" % (0.5*(bins[i]+bins[i+1])) for i in range(10)))
for case in args:
    t = latest_time(case)
    C, _ = read_field(os.path.join(case, "0", "C")); p, _ = read_field(os.path.join(case, t, "p"))
    a, _ = read_field(os.path.join(case, t, "ibmAlpha"))
    f = os.path.join(case, t, "ibmValueOffset")
    if not os.path.exists(f): f = os.path.join(case, t, "ibmCentroidOffset")
    off = read_field(f)[0]
    x = C[:, 0] + off[:, 0] - cx; y = C[:, 1] + off[:, 1] - cy; r = np.hypot(x, y)
    fl = a > 0; e = p - pex(r); e -= np.sum(e[fl]*a[fl])/np.sum(a[fl])
    s = "%-26s %9.2e" % ("/".join(case.split("/")[-2:]), np.sqrt(np.sum(e[fl]**2*a[fl])/np.sum(a[fl])))
    for i in range(10):
        m = fl & (r >= bins[i]) & (r < bins[i+1]); s += "%9.1e" % np.sqrt(np.mean(e[m]**2))
    print(s)
