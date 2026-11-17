"""V14 DFG 3D-1Z (Re 20, steady) and 3D-3Z (U_m(t) = 2.25 sin(pi t/8), t in [0, 8]): cylinder
D = 0.1, axis z through (0.5, 0.2), across the square channel [0, 2.5] x [0, 0.41]^2 (Schaefer &
Turek 1996; references in validation/reference_dfg.py, specification doc/validationSuite.md §11).
Channel walls body-fitted, only the cylinder immersed: it meets the walls z = 0 and z = H at right
angles (cells with a wall face and a cut face, G1).

Mesh family n = 16, 24, 32, 48: background h0 = H/n in y and z, 2.5/round(2.5/h0) in x (aspect
ratio <= 1.005: 0.41 and 2.5 have no convenient common divisor), one level of refinement in the
box B1 and one more in the finest box B2, both fixed in space and symmetric about the channel
centre: the whole mesh scales with 1/n, D/h = 15.6, 23.4, 31.2, 46.8 at the cylinder. The boxes of
the spec (B2 = [0.4, 0.75] x [0.1, 0.3] x [0, H]) left 16-32 cells across the duct upstream, whose
error dominated C_D (order 3.6, extrapolated -0.17 %); with these C_D converges at order 1.96 to
-0.034 % of the reference (user, 2026-11-15). C_D and C_L are judged at D/h 46.8 (user,
2026-11-15): the lift, 0.15 % of the drag, is set by the cut pattern at D/h <= 31 (2-D DFG on the
same y-mesh: +114, -57, +29, +2.3, -4.8 % at D/h 15.6 ... 62.4), C_L within 10 %.
3D-3Z on n = 16, 24, 48 with a fixed time step from Courant <= 0.5 at |u| = 3.75 (user, 2026-11-15:
the spec's dt = 0.005 is Courant 2.3-3.5; |u| reaches 3.6 past the cylinder), and once with dt/2 on
n = 16. C_Dmax and C_Lmin judged at D/h 46.8 as in 3D-1Z (user, 2026-11-16: C_Dmax +4.65, +2.37 %
at D/h 15.6, 23.4, order 1.66; the lift unresolved below). C_Lmin is judged on the Richardson
extrapolation from D/h 23.4 and 46.8 with the formal order 2 (user, 2026-11-17: -0.0162, -0.0127
against -0.0110 at D/h 23.4, 46.8, the C_L(t) shape and the time of the minimum right at 46.8; the
lift at D/h 15.6 has the wrong sign, a D/h 94 run ~68 M cells), the finest mesh reported. The PIMPLE steady state differs from the
SIMPLEC one through the rAU in the Rhie-Chow flux (dt against relaxation): -0.5 % C_D at D/h 15.6,
order >= 3 in h in 2-D and 3-D, so the dt dependence at D/h 15.6 is a resolution effect."""
import os

import numpy as np

import common
import geometry
import meshes
import reference_dfg as dfg
from common import write_dict
from foamio import latest_time, read_field

ID, KIND, TIER = "V14", "case", 2
TITLE = "DFG 3D-1Z (Re 20) and 3D-3Z (Re(t) <= 100), cylinder across a square channel"
STRESSES, RELIES = ["K1", "K3", "K8"], ["V5", "G1"]
REGRESSION = (1e-4,)

LX, H, D, NU = 2.5, 0.41, 0.1, 1e-3
UM1 = 0.45                       # 3D-1Z: U_mean = 4/9 U_m = 0.2, Re 20
CENTRE = (0.5, 0.2, 0.0)
CYL = common.TUTORIALS / "cylinderRe1"
B2 = ((0.35, 0.05, 0.0), (0.85, 0.36, H))
B1 = ((0.25, 0.0, 0.0), (1.4, H, H))
P1, P2 = (0.45, 0.2, 0.205), (0.55, 0.2, 0.205)
NS = [16, 24, 32, 48]
NP = {16: 16, 24: 32, 32: 64, 48: 128}    # ranks per run (0.32, 1.06, 2.5, 8.5 M cells)
NS3 = [16, 24, 48]                   # 3D-3Z
UMAX3, CO3, END3 = 3.75, 0.5, 8.0  # |u| bound for the time step: U_m 2.25 accelerates to 3.6 past the cylinder (V14.1)


def setup(d, n, Um, dtfactor=None):
    h0 = H/n
    nx = int(round(LX/h0))
    patches = {"inlet": ("patch", ["xmin"]), "outlet": ("patch", ["xmax"]),
               "walls": ("wall", ["ymin", "ymax", "zmin", "zmax"])}
    ncells = meshes.refined(d, (0, 0, 0), (LX, H, H), None, [B1, B2], patches, two_d=False,
                            cells=(nx, n, n))
    geometry.ibm_properties(d, {"cylinder": {"type": "cylinder", "centre": CENTRE, "axis": (0, 0, 1),
                                             "radius": D/2, "wallType": "noSlip"}})
    for x in ("transportProperties", "turbulenceProperties"):
        (d / "constant" / x).write_text((CYL / "constant" / x).read_text())
    common.set_keys(d / "constant" / "transportProperties", nu=NU)
    pc = meshes.face_centres(d)
    y, z = pc["inlet"][:, 1], pc["inlet"][:, 2]
    Uin = np.zeros((len(y), 3))
    Uin[:, 0] = 16*Um*y*z*(H - y)*(H - z)/H**4
    if dtfactor:   # 3D-3Z: U_m(t) = Um sin(pi t/8)
        inlet = ('type exprFixedValue; value uniform (0 0 0); valueExpr "vector(16*%r*sin(pi()*time()/8)'
                 '*pos().y()*pos().z()*(%r - pos().y())*(%r - pos().z())/%r, 0, 0)";' % (Um, H, H, H**4))
    else:
        inlet = "type fixedValue; value %s;" % meshes.vector_list(Uin)
    meshes.write_fields(d, {
        "U": ("volVectorField", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)",
              {"inlet": inlet, "walls": "type noSlip;",
               "outlet": "type zeroGradient;"}),
        "p": ("volScalarField", "[0 2 -2 0 0 0 0]", "uniform 0",
              {"inlet": "type zeroGradient;", "walls": "type zeroGradient;",
               "outlet": "type fixedValue; value uniform 0;"})}, two_d=False)
    hmin = h0/4
    if dtfactor:
        return ncells, hmin, unsteady_controls(d, hmin, dtfactor)
    (d / "system" / "fvSchemes").write_text((CYL / "system" / "fvSchemes").read_text())
    # as V5 2D-1: SIMPLEC, GAMG-preconditioned CG behind linkedGroup
    write_dict(d / "system" / "fvSolution", "dictionary", "fvSolution", """
solvers { p { solver PCG; preconditioner { preconditioner linkedGroup; base { preconditioner GAMG;
                smoother DICGaussSeidel; nCellsInCoarsestLevel 10; tolerance 0; relTol 0.01; } }
              tolerance 1e-12; relTol 0.01; maxIter 1000; }
          U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-12; relTol 0.1; } }
SIMPLE { nNonOrthogonalCorrectors 0; consistent yes; residualControl { p 1e-9; U 1e-9; } }
relaxationFactors { equations { U 0.9; } fields { p 1.0; } }
""")
    geometry.control_dict(d, "ibmSimpleFoam", endTime=20000, deltaT=1, writeInterval=20000)
    return ncells, hmin


def unsteady_controls(d, hmin, dtfactor):
    """3D-3Z as V5 2D-2: backward, PIMPLE 3 outer x 2 correctors, DIC-PCG behind linkedGroup; a
    fixed dt = dtfactor * Co h/|u|, an integer number of steps to t = 8"""
    t = (CYL / "system" / "fvSchemes").read_text().replace("default steadyState;", "default backward;")
    (d / "system" / "fvSchemes").write_text(t)
    write_dict(d / "system" / "fvSolution", "dictionary", "fvSolution", """
solvers { p { solver PCG; preconditioner { preconditioner linkedGroup; base DIC; } tolerance 1e-9; relTol 0.01; }
          pFinal { $p; relTol 0; }
          U { solver PBiCGStab; preconditioner { preconditioner linkedGroup; base DILU; } tolerance 1e-10; relTol 0; }
          UFinal { $U; } }
PIMPLE { momentumPredictor yes; nOuterCorrectors 3; nCorrectors 2; nNonOrthogonalCorrectors 0;
         consistent yes; ddtCorr yes; }
""")
    steps = int(np.ceil(END3/(dtfactor*CO3*hmin/UMAX3)))
    dt = END3/steps
    geometry.control_dict(d, "ibmPimpleFoam", endTime=END3, deltaT=repr(float(dt)), writeInterval=steps)
    return dt


def decomposed(d, app, np_, times):
    write_dict(d / "system" / "decomposeParDict", "dictionary", "decomposeParDict",
               "numberOfSubdomains %d;\nmethod scotch;\n%s" % (np_, common.LINK_CONSTRAINT))
    common.link_zone(d)
    common.sh("decomposePar", d, "log.decomposePar")
    common.sh(app, d, "log." + app, np_=np_)
    common.sh("reconstructPar " + times, d, "log.reconstructPar")
    common.sh("ibmSetGeometry", d, "log.ibmSetGeometry")   # geometry fields at 0 for the sampling


def point_value(xc, a, f, x0, h, radius=4.0):
    """value at x0 by an alpha-weighted quadratic least-squares fit in x, y, z over the fluid
    cells within radius*h (10 coefficients; a wall point has ~130 cells within 4 h)"""
    dx = (xc - np.asarray(x0, float))/h
    m = (a > 0) & (np.linalg.norm(dx, axis=1) < radius)
    X = dx[m]
    cols = [np.ones(len(X))] + [X[:, i] for i in range(3)] + [X[:, i]*X[:, j] for i in range(3) for j in range(i, 3)]
    A = np.column_stack(cols)*np.sqrt(a[m])[:, None]
    c, *_ = np.linalg.lstsq(A, f[m]*np.sqrt(a[m]), rcond=None)
    return float(c[0])


def courant(log):
    import re
    return float(max(float(x) for x in re.findall(r"Courant Number mean: \S+ max: (\S+)", log)))


def unsteady(ctx, n, dtfactor):
    d = ctx.root / ("3z_n%d_dt%g" % (n, dtfactor))
    ncells, hmin, dt = setup(d, n, 2.25, dtfactor)
    decomposed(d, "ibmPimpleFoam", NP[n], "-latestTime")
    t, F = meshes.force_history(d, 1.0)
    c = 2/(1.0*D*H)                                   # U_mean = 4/9 2.25 = 1
    cd, cl = c*F[:, 0], c*F[:, 1]
    tr, cdr, clr = dfg.bench3d_3z()
    l2 = lambda a, b: float(np.sqrt(np.mean((np.interp(tr, t, a) - b)**2)/np.mean(b**2)))
    return {"n": n, "Dh": D/hmin, "dtfactor": dtfactor, "dt": dt, "cells": ncells,
            "Co": courant((d / "log.ibmPimpleFoam").read_text()),
            "C_Dmax": float(cd.max()), "t_C_Dmax": float(t[cd.argmax()]),
            "C_Lmin": float(cl.min()), "t_C_Lmin": float(t[cl.argmin()]), "C_Lmax": float(cl.max()),
            "L2_C_D": l2(cd, cdr), "L2_C_L": l2(cl, clr),
            "series": [t[::10].tolist(), cd[::10].tolist(), cl[::10].tolist()]}


def steady(ctx, n):
    d = ctx.root / ("re20_n%d" % n)
    ncells, hmin = setup(d, n, UM1)
    decomposed(d, "ibmSimpleFoam", NP[n], "-latestTime")
    log = (d / "log.ibmSimpleFoam").read_text()
    U = 4/9*UM1
    c = 2/(U*U*D*H)
    F = common.ibm_force(log, "cylinder")["total"]
    t, Fh = meshes.force_history(d, 1.0)
    cd = c*Fh[:, 0]
    t = latest_time(str(d))
    xc, a = meshes.fluid_cells(str(d), t)
    p, _ = read_field(os.path.join(d, t, "p"))
    dp = point_value(xc, a, p, P1, hmin) - point_value(xc, a, p, P2, hmin)
    return {"n": n, "Dh": D/hmin, "cells": ncells, "C_D": float(c*F[0]), "C_L": float(c*F[1]),
            "dp": dp, "iterations": int(float(t)),
            "C_Dchange": float(np.ptp(cd[-200:])/abs(cd[-1])),   # stationarity over the last 200
            "converged": "SIMPLE solution converged" in log}


def run(ctx, res):
    ns = NS[:1] if ctx.quick else NS
    jobs = [(steady, (n,)) for n in ns] + [(unsteady, (n, 1.0)) for n in (NS3[:1] if ctx.quick else NS3)]
    jobs += [] if ctx.quick else [(unsteady, (NS3[0], 0.5))]
    # 304 ranks in all: as many runs at a time as --jobs (cores) allows
    import threading
    free, cv = [max(ctx.jobs, max(NP[a[0]] for _, a in jobs))], threading.Condition()

    def one(j):
        need = NP[j[1][0]]
        with cv:
            cv.wait_for(lambda: free[0] >= need)
            free[0] -= need
        try:
            return j[0].__name__, j[0](ctx, *j[1])
        finally:
            with cv:
                free[0] += need
                cv.notify_all()
    jobs.sort(key=lambda j: (j[0] is steady, -NP[j[1][0]]))   # the long transient runs first
    out = common.pmap(one, jobs, len(jobs))
    st = sorted([r for k, r in out if k == "steady"], key=lambda r: r["n"])
    un = sorted([r for k, r in out if k == "unsteady"], key=lambda r: (r["n"], -r["dtfactor"]))
    res.data["steady"], res.data["unsteady"] = st, un
    res.data["metrics"] = {"re20_n%d_%s" % (r["n"], q): r[q] for r in st for q in ("C_D", "C_L", "dp")}
    res.data["metrics"].update({"3z_n%d_dt%g_%s" % (r["n"], r["dtfactor"], q): r[q] for r in un
                                for q in ("C_Dmax", "C_Lmin", "C_Lmax")})
    res.le("3D-1Z runs not converged", sum(not r["converged"] for r in st), 0,
           "relative change of C_D over the last 200 iterations %s"
           % " ".join("%.1e" % r["C_Dchange"] for r in st))
    seq = lambda q: " ".join("%.5g" % r[q] for r in st)
    B = dfg.BENCH3D_1Z
    f = st[-1]
    res.rel("3D-1Z C_D, D/h = %.1f" % f["Dh"], f["C_D"], B["C_D"], 0.01, seq("C_D"))
    res.rel("3D-1Z C_L, D/h = %.1f" % f["Dh"], f["C_L"], B["C_L"], 0.10, seq("C_L"))
    res.report("3D-1Z dp, D/h = %.1f (no sourced reference)" % f["Dh"], f["dp"], note=seq("dp"))
    if len(st) >= 3:
        for q in ("C_D", "C_L"):
            Ce, p = common.richardson([1/r["Dh"] for r in st], [r[q] for r in st])
            res.report("3D-1Z observed order of %s" % q, p,
                       note="extrapolated %.6g (reference %.6g, %+.3f %%)" % (Ce, B[q], 100*(Ce/B[q] - 1)))
    u1 = [r for r in un if r["dtfactor"] == 1.0]
    if not u1:
        return
    B3 = dfg.BENCH3D_3Z
    g = u1[-1]
    sq = lambda q: " ".join("%.5g" % r[q] for r in u1)
    res.le("3D-3Z max Courant number", max(r["Co"] for r in un), CO3,
           " ".join("%.3f" % r["Co"] for r in un))
    res.rel("3D-3Z C_Dmax, D/h = %.1f" % g["Dh"], g["C_Dmax"], B3["C_Dmax"], 0.01, sq("C_Dmax"))
    if len(u1) >= 2:   # user (2026-11-17): judged extrapolated, as V11
        a = u1[-2]
        r = g["Dh"]/a["Dh"]
        ext = lambda p: g["C_Lmin"] + (g["C_Lmin"] - a["C_Lmin"])/(r**p - 1)
        res.rel("3D-3Z C_Lmin, Richardson-extrapolated from D/h %.1f/%.1f with the formal order 2"
                % (a["Dh"], g["Dh"]), ext(2.0), B3["C_Lmin"], 0.10,
                sq("C_Lmin") + "; order 1.5 / 2.5: %.5g / %.5g; D/h %.1f %+.1f %%"
                % (ext(1.5), ext(2.5), g["Dh"], 100*(g["C_Lmin"]/B3["C_Lmin"] - 1)))
        res.report("3D-3Z C_Lmin, D/h = %.1f (finest mesh)" % g["Dh"], g["C_Lmin"], B3["C_Lmin"], sq("C_Lmin"))
    else:
        res.rel("3D-3Z C_Lmin, D/h = %.1f" % g["Dh"], g["C_Lmin"], B3["C_Lmin"], 0.10, sq("C_Lmin"))
    res.report("3D-3Z C_Lmax, D/h = %.1f" % g["Dh"], g["C_Lmax"], B3["C_Lmax"], sq("C_Lmax"))
    if len(u1) >= 3:
        Ce, p = common.richardson([1/r["Dh"] for r in u1], [r["C_Dmax"] for r in u1])
        res.report("3D-3Z observed order of C_Dmax", p,
                   note="extrapolated %.6g (reference %.6g, %+.3f %%)" % (Ce, B3["C_Dmax"], 100*(Ce/B3["C_Dmax"] - 1)))
    for q in ("C_D", "C_L"):
        res.report("3D-3Z relative L2 difference of %s(t) to FeatFlow level 4, D/h = %.1f" % (q, g["Dh"]),
                   g["L2_" + q], note=sq("L2_" + q))
    h = [r for r in un if r["n"] == NS3[0]]
    if len(h) == 2:
        a, b = h
        res.report("3D-3Z time step: dt/2 against dt at D/h = %.1f, C_Dmax, C_Lmin (rel.)" % a["Dh"],
                   abs(b["C_Dmax"]/a["C_Dmax"] - 1),
                   note="C_Dmax %+.2e, C_Lmin %+.2e" % (b["C_Dmax"]/a["C_Dmax"] - 1, b["C_Lmin"]/a["C_Lmin"] - 1))
