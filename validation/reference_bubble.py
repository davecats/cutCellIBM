"""Reference data for the stress-free (clean, spherical or spheroidal) bubble of fixed shape in uniform
flow (case BUB): zero tangential stress and no penetration at the surface, U = 1,
Re = 2 U r_eq/nu with r_eq = (b^2 a)^(1/3) the radius of the sphere of equal volume,
C_D = F_x/(rho U^2 pi r_eq^2/2).

- Hadamard-Rybczynski (Re -> 0): C_D = 16/Re.
- Moore 1963 (JFM 16, 161) boundary layer theory, sphere ("M1"): C_D = 48/Re (1 - 2.21/Re^0.5),
  and Moore 1965 for oblate spheroids ("M2", only through table III below).
- Blanco & Magnaudet 1995 (Phys. Fluids 7, 1265-1274, doi:10.1063/1.868515): axisymmetric finite
  volume on orthogonal body-fitted grids, drag to < 0.5 % (their grid study at Re = 1000);
  table I (sphere), table III (oblate spheroids of aspect ratio chi = b/a).

Kinematic identity used as a pointwise check: on a stress-free sphere of radius a the surface
vorticity is omega = 2 u_theta/a at every point and every Re (zero shear: d u_theta/dr = u_theta/r;
u_r = 0). The free-slip condition d u_t/dn = 0 gives u_theta/a instead. For Re -> infinity
omega_max -> 3 U/a (potential flow, Blanco & Magnaudet).
"""
import numpy as np

RE_TABLE = (50, 100, 200, 300, 400, 600, 800, 1000)
SPHERE = {  # table I
    "this study": (0.664, 0.369, 0.200, 0.138, 0.104, 0.072, 0.055, 0.045),
    "Moore": (0.660, 0.374, 0.202, 0.140, 0.107, 0.073, 0.055, 0.045),
    "Ryskin & Leal 1984": (0.670, 0.380, 0.220, None, None, None, None, None),
}
SPHEROID = {  # table III: chi -> (Moore M2, this study)
    1.50: ((1.130, 0.617, 0.327, 0.223, 0.170, 0.115, 0.087, 0.070),
           (0.802, 0.516, 0.296, 0.209, 0.161, 0.112, 0.084, 0.068)),
    1.65: ((1.333, 0.715, 0.374, 0.255, 0.193, 0.131, 0.099, 0.079),
           (0.863, 0.549, 0.320, 0.227, 0.173, 0.122, 0.094, 0.076)),
    1.75: ((1.488, 0.787, 0.409, 0.277, 0.210, 0.142, 0.107, 0.086),
           (0.924, 0.609, 0.357, 0.251, 0.195, 0.132, 0.102, 0.083)),
    1.84: ((1.643, 0.859, 0.443, 0.299, 0.226, 0.152, 0.115, 0.092),
           (0.950, 0.594, 0.359, 0.259, 0.202, 0.141, 0.108, 0.088)),
    1.95: ((1.854, 0.954, 0.487, 0.328, 0.247, 0.166, 0.125, 0.100),
           (0.979, 0.627, 0.386, 0.277, 0.220, 0.152, 0.117, 0.095)),
}
# Dandy & Leal 1986 for the sphere at Re = 300: 0.143 (quoted by Blanco & Magnaudet, section V)


def hadamard(Re):
    return 16.0/np.asarray(Re, float)


def moore(Re):
    Re = np.asarray(Re, float)
    return 48.0/Re*(1.0 - 2.21/np.sqrt(Re))


def drag(Re, chi=1.0):
    """numerical C_D of Blanco & Magnaudet at a tabulated Re"""
    i = RE_TABLE.index(Re)
    return SPHERE["this study"][i] if chi == 1.0 else SPHEROID[chi][1][i]


def check():
    print("sphere: Moore formula against the Moore row of table I:",
          " ".join("%.3f/%.3f" % (moore(Re), m) for Re, m in zip(RE_TABLE, SPHERE["Moore"])))
    print("sphere: numerical/Moore - 1 [%]:",
          " ".join("%+.1f" % (100*(n/m - 1)) for n, m in zip(SPHERE["this study"], SPHERE["Moore"])))


if __name__ == "__main__":
    check()
