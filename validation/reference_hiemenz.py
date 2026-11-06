"""Exact Navier-Stokes solutions of stagnation-point flow onto a plane wall (gate G9): planar
Hiemenz flow and axisymmetric Homann flow, strain rate a, kinematic viscosity nu, rho = 1.

Wall coordinates: s along the wall (2-D) or r = |s| in the wall plane (3-D), n normal into the
fluid, eta = n sqrt(a/nu).
- Hiemenz:  u_s = a s f'(eta), u_n = -sqrt(a nu) f(eta),  f''' + f f'' - f'^2 + 1 = 0,
            p = p0 - a^2 s^2/2 - a nu (f' + f^2/2),       f''(0) = 1.232588 (White, Viscous
            Fluid Flow), f -> eta - 0.6479.
- Homann:   u_r = a r f'(eta), u_n = -2 sqrt(a nu) f(eta), f''' + 2 f f'' - f'^2 + 1 = 0,
            p = p0 - a^2 r^2/2 - 2 a nu (f' + f^2),       f''(0) = 1.311938, f -> eta - 0.5689.
f(0) = f'(0) = 0, f'(inf) = 1. Wall shear stress (rho nu du_s/dn at n = 0) = a s sqrt(a nu) f''(0),
wall pressure p0 - a^2 s^2/2.

`exact(points, ...)` evaluates U and p in global coordinates for a wall through `origin` with unit
normal `normal` (pointing into the fluid), so that tilted walls (tan alpha = 1/3, 1) use the same
function. In 2-D the tangent is the normal rotated by -90 deg in the (x, y) plane.
"""
import numpy as np
from scipy.integrate import solve_bvp, solve_ivp

ETA_MAX = 12.0
LITERATURE = {"hiemenz": {"f''(0)": 1.232588, "delta*": 0.6479},
              "homann": {"f''(0)": 1.311938, "delta*": 0.5689}}
_SOL = {}


def solve(kind="hiemenz"):
    """returns (callable eta -> (f, f', f''), f''(0), displacement constant lim eta - f)"""
    if kind in _SOL:
        return _SOL[kind]
    c = {"hiemenz": 1.0, "homann": 2.0}[kind]

    def rhs(eta, y):
        return np.vstack([y[1], y[2], -c*y[0]*y[2] + y[1]**2 - 1.0])

    def bc(ya, yb):
        return np.array([ya[0], ya[1], yb[1] - 1.0])

    eta = np.linspace(0.0, ETA_MAX, 400)
    guess = np.vstack([eta - 1.0 + np.exp(-eta), 1.0 - np.exp(-eta), np.exp(-eta)])
    r = solve_bvp(rhs, bc, eta, guess, tol=1e-10, max_nodes=100000)
    assert r.success, r.message
    fpp0 = float(r.sol(0.0)[2])
    dstar = float(ETA_MAX - r.sol(ETA_MAX)[0])

    # below the wall (eta < 0): the same ODE integrated backwards from the wall values, a smooth
    # divergence-free continuation for boundary data on faces inside the solid
    back = solve_ivp(lambda e, y: rhs(e, y[:, None])[:, 0], (0.0, -3.0), [0.0, 0.0, fpp0],
                     dense_output=True, rtol=1e-12, atol=1e-14)

    def f(e):
        e = np.asarray(e, float)
        inner = r.sol(np.clip(e, 0.0, ETA_MAX))
        outer = np.array([e - dstar, np.ones_like(e), np.zeros_like(e)])
        below = back.sol(np.clip(e, -3.0, 0.0))
        return np.where(e < 0.0, below, np.where(e <= ETA_MAX, inner, outer))

    _SOL[kind] = (f, fpp0, dstar)
    return _SOL[kind]


def exact(points, a=1.0, nu=0.01, origin=(0, 0, 0), normal=(0, 1, 0), kind="hiemenz", p0=0.0):
    """U (N, 3) and p (N,) at points (N, 3); kind 'hiemenz' (2-D, in the plane z = const) or
    'homann' (3-D)"""
    f, _, _ = solve(kind)
    x = np.atleast_2d(np.asarray(points, float)) - np.asarray(origin, float)
    nhat = np.asarray(normal, float)/np.linalg.norm(normal)
    n = x @ nhat
    sv = x - np.outer(n, nhat)  # in-plane position vector
    if kind == "hiemenz":
        that = np.array([nhat[1], -nhat[0], 0.0])
        sv = np.outer(sv @ that, that)
    s2 = np.einsum("ij,ij->i", sv, sv)
    F, Fp, _ = f(n*np.sqrt(a/nu))
    c = 1.0 if kind == "hiemenz" else 2.0
    U = a*sv*Fp[:, None] - c*np.sqrt(a*nu)*np.outer(F, nhat)
    p = p0 - 0.5*a*a*s2 - c*a*nu*(Fp + 0.5*c*F*F)
    return U, p


def wall_shear(s, a=1.0, nu=0.01, kind="hiemenz"):
    """tangential wall shear stress (kinematic) at wall distance s from the stagnation point"""
    return a*np.asarray(s, float)*np.sqrt(a*nu)*solve(kind)[1]


def check():
    for kind, lit in LITERATURE.items():
        _, fpp0, dstar = solve(kind)
        print("%-8s f''(0) = %.6f (literature %.6f), delta* = %.4f (%.4f)"
              % (kind, fpp0, lit["f''(0)"], dstar, lit["delta*"]))
    # residual of the momentum equations by finite differences, tilted wall tan alpha = 1/3
    nrm = np.array([-1.0, 3.0, 0.0])/np.sqrt(10.0)
    for kind in LITERATURE:
        h, nu = 1e-4, 0.01
        x0 = np.array([0.3, 0.25, 0.1 if kind == "homann" else 0.0])
        dims = 3 if kind == "homann" else 2

        def field(x):
            U, p = exact(x[None, :], nu=nu, normal=nrm, kind=kind)
            return U[0], p[0]

        U0, _ = field(x0)
        res = np.zeros(3)
        div = 0.0
        for d in range(dims):
            e = np.zeros(3)
            e[d] = h
            Up, pp = field(x0 + e)
            Um, pm = field(x0 - e)
            res += U0[d]*(Up - Um)/(2*h) - nu*(Up - 2*U0 + Um)/h**2
            res[d] += (pp - pm)/(2*h)
            div += (Up[d] - Um[d])/(2*h)
        print("%-8s tilted wall, momentum residual %.1e, divergence %.1e (finite differences)"
              % (kind, np.abs(res).max(), abs(div)))


if __name__ == "__main__":
    check()
