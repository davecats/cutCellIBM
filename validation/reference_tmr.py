"""Reference data of the NASA Turbulence Modeling Resource (TMR) verification cases with the
Spalart-Allmaras model (cases V15 2DZP and V16 2DBC), https://tmbwg.github.io/turbmodels/
(flatplate_sa.html, bump_sa.html; data files copied unchanged into data/tmr, see the README).

CFL3D (cell-centred, structured) and FUN3D (node-centred) on the five nested TMR grids, M = 0.2
(compressible), SA as defined on the TMR page, free-stream nuTilda = 3 nu; incompressible results
exist only for C_f(0.97) of the flat plate (FUN3D incompressible, SC/Tetra). The TMR note on
OpenFOAM: its SA (2014) did not match the TMR equations, so the primary reference of the IBM runs
is body-fitted OpenFOAM with the same model on the TMR grid (as for the NACA case, V8), and the
TMR values check that body-fitted reference.

- 2DZP: Re = 5e6 per unit length, plate x in [0, 2] (symmetry plane upstream from x = -0.33333),
  top at y = 1; C_f at x = 0.97008, C_D of the plate with reference length 2; u+(y+) and u(y) at
  x = 0.97008 and 1.90334 (CFL3D, finest grid 545 x 385); nu_t/nu at x = 0.97.
- 2DBC: Re = 3e6 per unit length, wall x in [0, 1.5], bump y = 0.05 sin^4(pi x/0.9 - pi/3) for
  0.3 <= x <= 1.2, symmetry planes from x = -25 and to x = 26.5, top at y = 5; C_L, C_D (pressure
  and viscous) with reference length 1.5; C_f at x = 0.6321975, 0.75, 0.8678025; C_p(x), C_f(x)
  (finest grid 1409 x 641); u(y - y0) at x = 0.75 and 1.20148; nu_t/nu at x = 0.75.
"""
import os
import re
import numpy as np

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "tmr")

FLATPLATE = {"Re_per_length": 5e6, "x_cf": 0.970084071, "ref_length": 2.0, "nuTilda_inf/nu": 3.0}
BUMP = {"Re_per_length": 3e6, "x_cf": (0.6321975, 0.75, 0.8678025), "ref_length": 1.5,
        "nuTilda_inf/nu": 3.0}


def bump_y(x):
    x = np.asarray(x, float)
    return np.where((x >= 0.3) & (x <= 1.2), 0.05*np.sin(np.pi*x/0.9 - np.pi/3)**4, 0.0)


def read_zones(name):
    """Tecplot ASCII file of the TMR -> (variable names, dict zone title -> array)"""
    names, zones, cur = [], {}, None
    with open(os.path.join(DATA, name)) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            if s.lower().startswith("variables"):
                names = re.findall(r'"([^"]*)"', s)
            elif s.lower().startswith("zone"):
                cur = re.search(r'[tT]\s*=\s*"([^"]*)"', s).group(1).strip()
                zones[cur] = []
            else:
                zones[cur].append([float(v) for v in s.split()])
    return names, {k: np.array(v) for k, v in zones.items()}


def extrapolate(h, q):
    """Richardson extrapolation from the three finest grids (h ascending): q_h0, observed order"""
    i = np.argsort(h)[:3]
    h, q = np.asarray(h)[i], np.asarray(q)[i]
    r = h[1]/h[0]
    e21, e32 = q[1] - q[0], q[2] - q[1]
    p = np.log(abs(e32/e21))/np.log(r)
    return q[0] - e21/(r**p - 1), p


def flatplate_cf():
    """dict code -> (h, C_f(0.97)) for CFL3D, FUN3D (compressible) and the incompressible codes"""
    _, z = read_zones("FlatPlate_SA_cf_convergence.dat")
    _, zi = read_zones("FlatPlate_validation_cf_incomp_results_sa.dat")
    z.update(zi)
    return {k: (v[:, 2], v[:, 3]) for k, v in z.items()}


def flatplate_drag():
    _, z = read_zones("FlatPlate_SA_drag_convergence.dat")
    return {k: (v[:, 2], v[:, 3]) for k, v in z.items()}


def bump_forces():
    """dict code -> dict h, C_L, C_D, C_Dp, C_Dv"""
    _, z = read_zones("Bump_SA_force_convergence.dat")
    return {k: dict(zip(("h", "C_L", "C_D", "C_Dp", "C_Dv"), v[:, 2:].T)) for k, v in z.items()}


def bump_cf():
    """dict (code, x) -> (h, C_f)"""
    _, z = read_zones("Bump_SA_cf_convergence_sa.dat")
    out = {}
    for k, v in z.items():
        code, x = k.split(", x=")
        out[(code, float(x))] = (v[:, 2], v[:, 3])
    return out


def check():
    print("2DZP C_f(0.97), finest grid and Richardson extrapolation of the finest three:")
    for k, (h, q) in flatplate_cf().items():
        qe, p = extrapolate(h, q)
        print("  %-26s %.6f -> %.6f (order %.2f)" % (k, q[np.argmin(h)], qe, p))
    for k, (h, q) in flatplate_drag().items():
        qe, p = extrapolate(h, q)
        print("2DZP C_D %-8s %.6f -> %.6f (order %.2f)" % (k, q[np.argmin(h)], qe, p))
    for k, d in bump_forces().items():
        print("2DBC %-6s finest: " % k + " ".join("%s %.6f" % (n, d[n][0]) for n in
                                                     ("C_L", "C_D", "C_Dp", "C_Dv")))
    for (code, x), (h, q) in bump_cf().items():
        print("2DBC C_f(%.4f) %-6s %.6f" % (x, code, q[np.argmin(h)]))
    _, cp = read_zones("Bump_SA_cp_bump.dat")
    c = cp["CFL3D"]
    print("2DBC min C_p %.4f at x = %.4f (CFL3D); bump height at 0.75: %.4f"
          % (c[:, 1].min(), c[c[:, 1].argmin(), 0], bump_y(0.75)))


if __name__ == "__main__":
    check()
