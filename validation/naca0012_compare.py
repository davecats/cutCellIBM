#!/usr/bin/env python
"""Lift and drag of NACA 0012 runs (tutorials/turbulence/naca0012/makeCase.py).

    naca0012_compare.py [--cv] case1 case2 ...

The angle of attack is read from the inflow direction (0/U), the forces from
postProcessing/forces/0/force.dat (body-fitted) or
postProcessing/ibmForces/<first time>/forces.dat
(IBM), per unit density for the depth
of the mesh (DZ = 0.01): C = F . e / (0.5 U^2 c DZ), e the lift or drag
direction. Prints the final coefficients, their pressure and viscous parts,
and the change of C_l and C_d over the last 500 iterations (convergence),
relative to the first case.

--cv adds C_l and C_d from the momentum balance of a control volume of whole
cells around the airfoil (cvforce.py, latest written time) for the box
-1 < x < 2, |y| < 1, and the spread of C_d over five boxes from
(-0.25, 1.25) x (-0.25, 0.25) to (-5, 10) x (-5, 5). The discrete fluxes
through the interior faces cancel, so this equals the force on the body when
the solution is in momentum balance; a spread larger than the iteration
change marks cells (usually thin, skewed ones in the wake) that are not.
"""
import math, os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

DZ = 0.01
CV_BOXES = [(-0.25, 1.25, -0.25, 0.25), (-0.5, 1.5, -0.5, 0.5), (-1, 2, -1, 1), (-2, 4, -2, 2), (-5, 10, -5, 5)]


def alpha_of(case):
    s = open(os.path.join(case, "0", "U")).read()
    ux, uy = map(float, re.search(r"internalField\s+uniform\s+\(\s*(\S+)\s+(\S+)", s).groups())
    return math.atan2(uy, ux)


def history(case):
    """Iterations x 4: pressure (x, y), viscous (x, y); restarts (one
    directory per start time) are joined."""
    d = os.path.join(case, "postProcessing", "forces")
    ibm = not os.path.isdir(d)
    if ibm:
        d = os.path.join(case, "postProcessing", "ibmForces")
    rows = []
    for t in sorted(os.listdir(d), key=float):
        f = os.path.join(d, t, "forces.dat" if ibm else "force.dat")
        for line in open(f):
            if line.startswith("#"):
                continue
            v = [float(x) for x in re.sub(r"[()]", " ", line).split()]
            rows.append([v[1], v[2], v[4], v[5]] if ibm else [v[4], v[5], v[7], v[8]])
    return np.array(rows)


def main():
    args = sys.argv[1:]
    cv = "--cv" in args
    if cv:
        from cvforce import cv_force
        args.remove("--cv")
    ref = None
    for case in args:
        a = alpha_of(case)
        eD, eL = np.array([math.cos(a), math.sin(a)]), np.array([-math.sin(a), math.cos(a)])
        h = history(case) / (0.5 * DZ)
        Fp, Fv = h[:, 0:2], h[:, 2:4]
        Cl_p, Cl_v = Fp @ eL, Fv @ eL
        Cd_p, Cd_v = Fp @ eD, Fv @ eD
        Cl, Cd = Cl_p + Cl_v, Cd_p + Cd_v
        k = min(500, len(Cl) - 1)
        line = ("%-22s alpha %4.1f  its %5d  Cl %.5f (p %.5f v %+.5f)  Cd %.6f (p %.6f v %.6f)"
                "  dCl/500 %.1e dCd/500 %.1e"
                % (os.path.basename(os.path.normpath(case)), math.degrees(a), len(Cl), Cl[-1], Cl_p[-1],
                   Cl_v[-1], Cd[-1], Cd_p[-1], Cd_v[-1], abs(Cl[-1] - Cl[-1 - k]), abs(Cd[-1] - Cd[-1 - k])))
        if ref is None:
            ref = (Cl[-1], Cd[-1])
        else:
            line += "  dCd %+.2f%%" % (100 * (Cd[-1] / ref[1] - 1))
            if abs(ref[0]) > 1e-3:
                line += " dCl %+.2f%%" % (100 * (Cl[-1] / ref[0] - 1))
        if cv:
            C = [(F[:2] @ eL / (0.5 * DZ), F[:2] @ eD / (0.5 * DZ)) for _, F, _ in cv_force(case, None, CV_BOXES)]
            cd = [c[1] for c in C]
            line += "  CV Cl %.5f Cd %.6f (spread %.1e)" % (C[2][0], C[2][1], max(cd) - min(cd))
        print(line)


if __name__ == "__main__":
    main()
