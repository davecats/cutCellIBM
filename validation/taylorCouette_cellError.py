"""Taylor-Couette errors by cell class (linking sweep, tutorials/order/taylorCouette with
mergeAlpha 0 / 0.5 and N = 40..320, cases named so<true|false>_m<mergeAlpha>_N<N>): all fluid cells, cut cells
(0 < alpha < 1), bulk (alpha = 1 and wall distance > 2h), and the wall-gradient error of the
cut cells: rms of |U - Uex|/d over the exact wall shear scale 2|B|/Ri^2. Values at the
fluid centroids when ibmCentroidOffset exists (at the value points, ibmValueOffset, when
that exists: group centroids of linked cells). Usage: taylorCouette_cellError.py [sweep dir]."""
import sys, os, glob, re
import numpy as np
sys.path.insert(0, "/home/ws/xt8786/Codes/openFoam.ibm/validation")
from foamio import read_field, latest_time

Ri, Ro, Om, cx, cy = 0.5, 1.5, 1.0, 1.6, 1.6
A = -Om*Ri**2/(Ro**2 - Ri**2); B = Om*Ri**2*Ro**2/(Ro**2 - Ri**2)


def errors(case):
    t = latest_time(case)
    C, _ = read_field(os.path.join(case, "0", "C"))
    U, _ = read_field(os.path.join(case, t, "U"))
    a, _ = read_field(os.path.join(case, t, "ibmAlpha"))
    f = os.path.join(case, t, "ibmValueOffset")
    if not os.path.exists(f): f = os.path.join(case, t, "ibmCentroidOffset")
    off = read_field(f)[0] if os.path.exists(f) else 0*C
    x = C[:, 0] + off[:, 0] - cx; y = C[:, 1] + off[:, 1] - cy; r = np.hypot(x, y)
    ut = A*r + B/r
    eU = np.hypot(U[:, 0] + ut*y/r, U[:, 1] - ut*x/r)
    n = int(round(np.sqrt(len(C)))); h = 3.2/n
    d = np.minimum(r - Ri, Ro - r)
    fl = a > 0; cut = fl & (a < 1); bulk = (a >= 1) & (d > 2*h)
    rms = lambda m: np.sqrt(np.sum(eU[m]**2*a[m])/np.sum(a[m]))
    tau = 2*abs(B)/Ri**2
    g = eU[cut]/np.maximum(d[cut], 1e-3*h)/tau
    return n, rms(fl), rms(cut), rms(bulk), np.sqrt(np.mean(g**2))


D = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
res = {}
for c in sorted(glob.glob(os.path.join(D, "so*_N*"))):
    if not os.path.exists(os.path.join(c, "0", "C")): continue
    so, m, n = re.match(r"so(\w+)_m([\d.]+)_N(\d+)", os.path.basename(c)).groups()
    res.setdefault((so, m), []).append(errors(c))
print("%-24s %4s %10s %5s %10s %5s %10s %5s %10s %5s" % ("variant", "N", "L2 all", "ord", "L2 cut", "ord",
                                                     "L2 bulk", "ord", "wall grad", "ord"))
for (so, m), rows in sorted(res.items()):
    rows.sort()
    for i, r in enumerate(rows):
        s = "%-24s %4d" % ("secondOrder %s, mergeAlpha %s" % (so, m) if i == 0 else "", r[0])
        for k in range(1, 5):
            o = np.log2(rows[i - 1][k]/r[k]) if i else np.nan
            s += " %10.3e %5s" % (r[k], "" if i == 0 else "%.2f" % o)
        print(s)
