"""Reference data for the circular cylinder rotating in a uniform stream (case ROT), Mittal & Kumar
2003 (JFM 476, 303-334, doi:10.1017/S0022112002002938): Re = 2 a U/nu = 200, rotation rate
alpha = a omega/U (counter-clockwise, flow in +x, so the lift is negative), stabilised FEM, outer
boundary at L = 100 D (C_D and St depend on L for L/D < 75), 160 points on the circumference,
dt U/a = 0.025.

C_D = F_x/(rho U^2 a), C_L = F_y/(rho U^2 a), C_M = M/(rho U^2 a^2), power coefficient
C_P = C_D + C_M alpha (section 2.2), so C_M = (C_P - C_D)/alpha: the torque. St = f 2a/U.

Figure 11 digitised (validation/data/mittalKumar2003, bitmapFigures_extract.py, ~0.3 % of each
axis range): time averages and r.m.s. of the periodic state, alpha = 0...4.3. C_M from the
difference of two digitised curves is uncertain by ~0.03/alpha (use alpha >= 1.5).
"""
import json
import os

import numpy as np

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "mittalKumar2003")

TEXT = {
    "alpha=0": {"St": 0.1934, "C_D": 1.316},  # Williamson 1991 (experiment): St 0.196
    "shedding stops": 1.91,                    # vortex shedding for alpha < 1.91 (DNS and linear stability)
    "second instability": (4.34, 4.70),        # one-sided shedding, unsteady again
}


def figure11():
    """dict of arrays alpha, C_L, C_D, C_P, C_Lrms, C_Drms, St (nan: not readable) and C_M"""
    t = np.loadtxt(os.path.join(DATA, "fig11_unsteady.dat"))
    d = dict(zip(["alpha", "C_L", "C_D", "C_P", "C_Lrms", "C_Drms", "St"], t.T))
    with np.errstate(invalid="ignore", divide="ignore"):
        d["C_M"] = np.where(d["alpha"] > 0, (d["C_P"] - d["C_D"])/d["alpha"], np.nan)
    return d


def bodyFitted(M=128):
    """the V12 problem body-fitted (suite/bodyfitted.py rot, decision D3): {alpha: row} with C_L,
    C_D, C_Lrms, C_Drms, St of the O-grid with M cells per quarter circle (64: 61k cells, 128: 243k;
    mean C_L at alpha = 1: -2.4984, -2.4994)"""
    with open(os.path.join(os.path.dirname(DATA), "rotatingCylinder_bodyFitted", "rot.json")) as f:
        rows = json.load(f)["rows"]
    return {r["alpha"]: r for r in rows if r["M"] == M}


def check():
    d = figure11()
    print("alpha = 0: C_D %.3f (text 1.316), C_P %.3f (= C_D), St %.4f" % (d["C_D"][0], d["C_P"][0], d["St"][0]))
    print("C_M/alpha (smooth if the digitisation is consistent):",
          " ".join("%.2f:%.3f" % (a, m/a) for a, m in zip(d["alpha"], d["C_M"]) if a >= 1 and np.isfinite(m)))
    print("for comparison, the steady vortex of a cylinder rotating in fluid at rest: C_M/alpha = 8 pi/Re = %.3f"
          % (8*np.pi/200))


if __name__ == "__main__":
    check()
