"""Reference data for the rigid sphere in uniform flow (case SPH): drag, wake geometry, loss of axial
symmetry and shedding, D = 1, U = 1, Re = U D/nu, C_D = F_x/(rho U^2 pi D^2/8).

- Johnson & Patel 1999 (JFM 378, 19-70, doi:10.1017/S0022112098003206), finite volume on a body-
  fitted O-O grid, ~1.7 order, figures 4, 5 and 10 digitised (validation/data/johnsonPatel1999,
  bitmapFigures_extract.py) and the numbers of the text (Re = 211, 212, 300).
- Magnaudet, Rivero & Fabre 1995 (JFM 284, 97-135) and Masliyah & Epstein 1970 (JFM 44, 493-512)
  as quoted by Blanco & Magnaudet 1995 (Phys. Fluids 7, 1265, table II, values in parentheses).
- Fornberg 1988 (JFM 190, 471-489), steady axisymmetric Navier-Stokes, second-order finite
  differences on successively refined grids, tabulated, errors "well below 1 %" (section 6).
  J&P give their steady C_D and wake only as figures (axis 0-3 for C_D: the digitised values
  carry +-0.01, ~1 % at Re 200); Fornberg's table runs ~2 % below them at Re 200.
"""
import os
import numpy as np

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "johnsonPatel1999")

# C_D of the solid sphere, Blanco & Magnaudet 1995 table II
CD_TABLE = {
    "Magnaudet et al. 1995": {1: 27.540, 10: 4.317, 20: 2.707, 50: 1.581, 100: 1.092},
    "Masliyah & Epstein 1970": {1: 26.74, 10: 4.34, 20: 2.75, 50: 1.59, 100: 1.08},
}

# Johnson & Patel 1999, text
REGIMES = {
    "separation onset": 20.0,          # theta_s = 180 deg, x_s = 0 at Re = 20 (figure 4)
    "loss of axial symmetry": 211.0,   # C_L -3e-5 at 211, -2.4e-2 at 212 (text, figure 10)
    "unsteady": 270.0,                 # steady non-axisymmetric for 210 < Re < 270 (text)
}
# Fornberg 1988: C_D from table 4 (finest grid), the wake length from table 2 (grid G1), given from
# the centre in radii L: the bubble behind the sphere x_s = (L - 1)/2 in diameters
FORNBERG = {100: {"C_D": 1.0852, "x_s": (2.747 - 1)/2},
            200: {"C_D": 0.7683, "x_s": (3.871 - 1)/2}}

RE300 = {  # periodic shedding, section 6.1
    "C_D mean": 0.656, "C_D amplitude": 3.5e-3,
    "C_L mean": -0.069, "C_L amplitude": 1.6e-2,  # sign of C_L arbitrary
    "St": 0.137,
    "others": "Tomboulides 1993: C_D 0.671, St 0.136; Roos & Willmarth 1971 (experiment): C_D 0.629",
}
# published computations at Re 300 (mean C_D, St): the spread is the reference band of V10.
# J&P and Tomboulides (spectral elements) from J&P section 6.1; Kang et al. and Angelidis et al. as
# quoted in the comparison table of arXiv 2411.00272 (2024)
RE300_PUBLISHED = {"Johnson & Patel 1999": (0.656, 0.137), "Tomboulides 1993": (0.671, 0.136),
                   "Kang et al.": (0.663, 0.134), "Angelidis et al.": (0.665, 0.132)}


def steady():
    """figures 4 and 5: dict of arrays Re, theta_s (deg), x_s, x_c, y_c (D), C_D"""
    t = np.loadtxt(os.path.join(DATA, "sphere_steady.dat"))
    return dict(zip(["Re", "theta_s", "x_s", "x_c", "y_c", "C_D"], t.T))


def lateral():
    """figure 10: Re, |C_L| of the steady non-axisymmetric wake"""
    t = np.loadtxt(os.path.join(DATA, "sphere_lateral.dat"))
    return t[:, 0], np.abs(t[:, 1])


def check():
    s = steady()
    for Re in (50, 100):
        i = list(s["Re"]).index(Re)
        print("Re %3d: C_D Johnson & Patel (digitised) %.3f, Magnaudet et al. %.3f, Masliyah & Epstein %.2f"
              % (Re, s["C_D"][i], CD_TABLE["Magnaudet et al. 1995"][Re], CD_TABLE["Masliyah & Epstein 1970"][Re]))
    Re, cl = lateral()
    print("|C_L| at Re 212: %.4f (text 0.024)" % cl[list(Re).index(212)])
    print("Re = 100: theta_s %.1f deg, x_s %.3f, (x_c, y_c) (%.3f, %.3f)" % tuple(
        s[k][list(s["Re"]).index(100)] for k in ("theta_s", "x_s", "x_c", "y_c")))


if __name__ == "__main__":
    check()
