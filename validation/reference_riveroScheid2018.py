"""Reference data of Rivero-Rodriguez & Scheid (2018), "Bubble dynamics in microchannels: inertial
and capillary migration forces", J. Fluid Mech. 842, 215-247, doi:10.1017/jfm.2018.78.

A sphere (bubble with a rigid or a stress-free interface) of diameter d in a periodic train in a
circular pipe of diameter d_h = 1, mean flow J = 1, transverse position eps (centre offset from the
axis), in the frame of the sphere (pipe wall moving at -V). Pure linear inertial regime (Re -> 0,
Ca = 0), appendix C:

    phi(d, eps) = sum_i sum_j phi_ij (eps/eps_*)^i d^j,   eps_* = (1 - d)/2,  (eps/eps_*)^0 = 1

    f1     = f/Re, the balanced transverse body force per unit bubble volume (= migration force
             L_y / V_B), in units mu J / d_h^2; f1 < 0 pushes the sphere towards the axis
    V0/VP  bubble velocity over the Poiseuille velocity at the centre, VP = 2 (1 - 4 eps^2)
    100 b0 pressure correction factor, Delta_s p_T / d_x p_P = 1 + beta alpha_G
    O0/8e  rotation rate over the Poiseuille one, Omega_P = 8 eps (rigid interface only)

The printed captions of tables 3 and 4 are swapped: table 3 (with Omega) is the rigid interface of
figure 7, table 4 the stress-free interface of figure 8. check() confirms this against the text
(centred position unstable for d < 0.73 rigid, d < 0.85 stress-free) and against figure 4.
Fits valid for 0.1 <= d <= 0.9, 0 <= eps <= 0.95 eps_*, error < 1 % (paper).
"""
import numpy as np

# rows j (power of d), columns i (power of eps/eps_*)
_RIGID = {
    "f1": ((1, 3, 5, 7), (1, 2, 3, 4, 5), [
        [13.28, -33.15, 2.20, -0.63],
        [-43.14, 115.27, 19.72, -28.76],
        [48.60, -104.32, -137.39, 122.41],
        [-20.20, -36.05, 255.14, -177.07],
        [0.57, 65.47, -150.87, 89.05]]),
    "V0/VP": ((0, 2, 4, 6), (0, 2, 4, 6, 8), [
        [1.00, -0.02, 0.10, -0.20],
        [-0.67, -0.63, -1.07, -0.47],
        [0.05, 2.42, 4.19, 5.87],
        [0.31, -3.47, -7.35, -14.90],
        [-0.22, 1.95, 4.88, 12.49]]),
    "100b0": ((0, 2, 4, 6), (0, 2, 4, 6, 8), [
        [-0.01, 0.14, -0.68, 0.45],
        [0.23, 3.65, 16.84, -9.38],
        [3.44, 60.06, -95.98, 116.73],
        [18.31, -130.40, 171.85, -146.87],
        [12.72, 105.94, -100.60, 50.42]]),
    "O0/8e": ((0, 2, 4, 6), (0, 2, 4, 6, 8), [
        [1.00, -0.02, 0.12, -0.20],
        [-0.14, -0.42, -0.93, -0.98],
        [-0.77, 1.30, 1.00, 10.63],
        [0.94, -2.40, 6.04, -31.08],
        [-0.62, 2.43, -10.59, 29.65]]),
}

_STRESSFREE = {
    "f1": ((1, 3, 5, 7), (1, 2, 3, 4, 5), [
        [3.55, -7.42, -0.10, -0.51],
        [-0.27, 16.33, 5.85, 1.27],
        [-4.67, 7.65, -7.98, -38.37],
        [0.45, -43.31, -14.05, 100.91],
        [-0.79, 28.12, 19.56, -68.44]]),
    "V0/VP": ((0, 2, 4, 6), (0, 2, 4, 6, 8), [
        [1.00, -0.01, 0.04, -0.07],
        [0.02, -0.17, 0.26, -1.56],
        [-0.32, -0.02, -2.70, 7.93],
        [-1.05, 0.32, 6.87, -15.38],
        [0.94, 0.07, -5.67, 11.28]]),
    "100b0": ((0, 2, 4, 6), (0, 2, 4, 6, 8), [
        [-0.03, 0.31, -1.11, 0.98],
        [0.95, -5.67, 31.40, -27.27],
        [-8.76, 73.27, -218.95, 205.73],
        [-6.43, -181.60, 511.37, -393.88],
        [18.21, 163.32, -396.80, 262.54]]),
}

FITS = {"rigid": _RIGID, "stressFree": _STRESSFREE}


def eps_star(d):
    return 0.5*(1.0 - np.asarray(d, float))


def fit(interface, quantity, d, eps):
    """phi(d, eps) of table 3 (interface 'rigid') or 4 ('stressFree')"""
    ipow, jpow, c = FITS[interface][quantity]
    d = np.asarray(d, float)
    s = np.asarray(eps, float)/eps_star(d)
    out = 0.0
    for jj, row in zip(jpow, c):
        for ii, cij in zip(ipow, row):
            out = out + cij*s**ii*d**jj
    return out


def f1(interface, d, eps):
    return fit(interface, "f1", d, eps)


def V0(interface, d, eps):
    return fit(interface, "V0/VP", d, eps)*2.0*(1.0 - 4.0*np.asarray(eps, float)**2)


def beta0(interface, d, eps):
    return fit(interface, "100b0", d, eps)/100.0


def Omega0(d, eps):
    return fit("rigid", "O0/8e", d, eps)*8.0*np.asarray(eps, float)


def slope_centre(interface, d):
    """d f1 / d eps at eps = 0 (> 0: centred position unstable)"""
    ipow, jpow, c = FITS[interface]["f1"]
    d = np.asarray(d, float)
    return sum(row[0]*d**jj for jj, row in zip(jpow, c))/eps_star(d)


def neutral(interface, d, n=2001):
    """equilibrium positions of the neutral sphere (f1 = 0, eps > 0) with their stability"""
    s = np.linspace(1e-3, 0.95, n)
    e = s*eps_star(d)
    f = f1(interface, d, e)
    out = []
    for k in np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0]:
        ek = e[k] - f[k]*(e[k + 1] - e[k])/(f[k + 1] - f[k])
        out.append((ek, "stable" if f[k + 1] < f[k] else "unstable"))
    return out


def check():
    """the fits against statements of the text"""
    from scipy.optimize import brentq
    for itf, dt in (("rigid", 0.73), ("stressFree", 0.85)):
        dc = brentq(lambda d: slope_centre(itf, d), 0.5, 0.95)
        print("%-10s centred position: unstable for d < %.3f (text %.2f)" % (itf, dc, dt))
    d = 0.4
    s = np.linspace(0, 0.95, 96)
    f = f1("rigid", d, s*eps_star(d))
    print("rigid d=0.4: max f/Re %.3f at eps/eps* %.2f (fig. 4a: ~0.3 at ~0.4), neutral %s" % (
        f.max(), s[np.argmax(f)], ["%.3f eps*=%.3f %s" % (e/eps_star(d), e, st) for e, st in neutral("rigid", d)]))
    print("rigid d=0.4 at eps=0: V0 %.3f (fig. 4b ~1.75), 100 beta0 %.3f, Omega0/8eps %.3f" % (
        V0("rigid", d, 0.0), 100*beta0("rigid", d, 0.0), fit("rigid", "O0/8e", d, 0.0)))
    for itf in ("rigid", "stressFree"):
        print("%-10s neutral positions:" % itf, "  ".join(
            "d=%.1f: %s" % (dd, ",".join("%.3f(%s)" % (e, st[0]) for e, st in neutral(itf, dd)))
            for dd in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)))
    print("small stress-free bubble: neutral eps -> sqrt(2)/4 = %.4f (text)" % (np.sqrt(2)/4))


if __name__ == "__main__":
    check()
