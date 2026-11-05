"""Reference data for the circular cylinder oscillating in-line in fluid at rest (case OSC), Dütsch,
Durst, Becker & Lienhart 1998 (JFM 360, 249-271, doi:10.1017/S002211209800860X), two-dimensional
finite volume, grid of concentric circles to 60 D, in the cylinder frame (body force -rho dU_s/dt,
outer boundary at -U_s).

Motion x(t) = -A sin(2 pi f t), U_max = 2 pi f A, Re = U_max D/nu, KC = U_max/(f D) = 2 pi A/D,
beta = Re/KC. In-line force per unit length in the inertial frame fitted with Morison's equation
(1.3):
    F(t) = -1/2 rho D c_d xdot |xdot| - pi/4 rho D^2 c_i xddot
(c_i: added mass; c_m = c_i + 1 includes the Froude-Krylov force that appears in an accelerated
frame or an oscillating flow). Paper: D = 10 mm, rho = 998.2 kg/m^3, mu = 1.00818 mPa s.
"""
import numpy as np

TABLE1 = {  # Re = 100, KC = 5: set -> (c_d, c_i, cells on the perimeter, cells total, dt/T)
    "A": (2.15, 1.42, 96, 6144, 1/720),
    "B1": (2.10, 1.45, 192, 24576, 1/720),
    "C": (2.09, 1.45, 384, 98304, 1/720),
    "B2": (2.15, 1.42, 192, 24576, 1/180),
    "B3": (2.17, 1.41, 192, 24576, 1/60),
}
# Re = 200, KC = 10 is not periodic (cycle-to-cycle variations): set D1 c_d 1.81-1.84, c_i 1.02-1.04

TABLE2_BETA35 = {  # KC -> (regime of Tatsuno & Bearman 1990, c_d, c_i); Re = 35 KC; ranges as (lo, hi)
    0.5: ("A*", 10.4, 1.37), 1.0: ("A*", 5.39, 1.36), 1.5: ("A*", 3.72, 1.35), 2.0: ("A*", 2.98, 1.34),
    3.0: ("A*", 2.28, 1.31), 4.0: ("A", 1.97, 1.30), 5.0: ("C", 1.82, 1.30), 6.0: ("E", 1.73, 1.17),
    7.0: ("E/F", (1.78, 1.79), (1.07, 1.15)), 8.0: ("F", (1.72, 1.73), (1.14, 1.15)),
    9.0: ("G", 1.80, 0.946), 15.0: ("G", (1.78, 1.81), (0.385, 0.511)),
}
REFERENCE = {"Re": 100, "KC": 5, "c_d": 2.09, "c_i": 1.45}  # set C, the finest


def morison_fit(t, F, f, A, D=1.0, rho=1.0):
    """least-squares c_d, c_i of Morison's equation from the in-line force F(t) per unit length
    (whole cycles of the periodic state) for x(t) = -A sin(2 pi f t)"""
    w = 2*np.pi*f
    xd = -A*w*np.cos(w*t)
    xdd = A*w*w*np.sin(w*t)
    M = np.column_stack([-0.5*rho*D*xd*np.abs(xd), -0.25*np.pi*rho*D*D*xdd])
    (cd, ci), *_ = np.linalg.lstsq(M, F, rcond=None)
    return cd, ci


def check():
    t = np.linspace(0, 3, 3001)
    A = 5/(2*np.pi)
    w = 2*np.pi
    xd, xdd = -A*w*np.cos(w*t), A*w*w*np.sin(w*t)
    F = -0.5*2.09*xd*np.abs(xd) - 0.25*np.pi*1.45*xdd
    print("Morison fit of a synthetic force (2.09, 1.45):", "%.4f %.4f" % morison_fit(t, F, 1.0, A))
    nu = 1.00818e-3/998.2
    print("paper units: Re 100, KC 5, D 10 mm: U_max %.4f m/s, A %.2f mm, f %.4f Hz"
          % (100*nu/0.01, 5/(2*np.pi)*10, 100*nu/0.01/(5*0.01)))


if __name__ == "__main__":
    check()
