"""V5 DFG 2D-1 (Re 20, steady) and 2D-2 (Re 100, periodic): cylinder D = 0.1 at (0.2, 0.2) in
a channel [0, 2.2] x [0, 0.41] (Schäfer & Turek 1996; references in validation/reference_dfg.py,
specification doc/validationSuite.md §11). Channel walls body-fitted, only the cylinder
immersed. Mesh family k: background h0 = 0.01/2^k, one level more in [0.05, 1.2] x [0, 0.41]
and one more in [0.1, 0.45] x [0.1, 0.3]: D/h = 40, 80, 160 at the cylinder for k = 0, 1, 2."""
import os
import re

import numpy as np

import common
import meshes
import reference_dfg as dfg
from common import write_dict
from foamio import latest_time, read_field

ID, KIND, TIER = "V5", "case", 1
TITLE = "DFG 2D-1 (Re 20) and 2D-2 (Re 100), cylinder in a channel"
STRESSES, RELIES = ["K7", "K3", "K5", "K8"], ["V2", "G9"]
REGRESSION = (1e-4,)

H, D, NU = 0.41, 0.1, 1e-3
CYL = common.TUTORIALS / "cylinderRe1"
BOXES = [((0.05, 0.0, 0.0), (1.2, H, 0.0)), ((0.1, 0.1, 0.0), (0.45, 0.3, 0.0))]
P1, P2 = (0.15, 0.2), (0.25, 0.2)
END = 16.0       # 2D-2: from rest; statistics over the last 3 s


def setup(d, k, Um, unsteady):
    h0 = 0.01/2**k
    patches = {"inlet": ("patch", ["xmin"]), "outlet": ("patch", ["xmax"]), "walls": ("wall", ["ymin", "ymax"])}
    ncells = meshes.refined(d, (0, 0, 0), (2.2, H, 0), h0, BOXES, patches)
    import geometry
    geometry.ibm_properties(d, {"cylinder": {"type": "cylinder", "centre": (0.2, 0.2, 0), "axis": (0, 0, 1),
                                             "radius": D/2, "wallType": "noSlip"}})
    for x in ("transportProperties", "turbulenceProperties"):
        (d / "constant" / x).write_text((CYL / "constant" / x).read_text())
    common.set_keys(d / "constant" / "transportProperties", nu=NU)
    pc = meshes.face_centres(d)
    y = pc["inlet"][:, 1]
    Uin = np.zeros((len(y), 3))
    Uin[:, 0] = 4*Um*y*(H - y)/H**2
    meshes.write_fields(d, {
        "U": ("volVectorField", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)",
              {"inlet": "type fixedValue; value %s;" % meshes.vector_list(Uin), "walls": "type noSlip;",
               "outlet": "type zeroGradient;"}),
        "p": ("volScalarField", "[0 2 -2 0 0 0 0]", "uniform 0",
              {"inlet": "type zeroGradient;", "walls": "type zeroGradient;", "outlet": "type fixedValue; value uniform 0;"})})
    hmin = h0/4
    if not unsteady:
        (d / "system" / "fvSchemes").write_text((CYL / "system" / "fvSchemes").read_text())
        write_dict(d / "system" / "fvSolution", "dictionary", "fvSolution", """
solvers { p { solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0.01; }
          U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-12; relTol 0.1; } }
SIMPLE { nNonOrthogonalCorrectors 0; consistent yes; residualControl { p 1e-9; U 1e-9; } }
relaxationFactors { equations { U 0.9; } fields { p 1.0; } }
""")
        import geometry as g
        g.control_dict(d, "ibmSimpleFoam", endTime=40000, deltaT=1)
        return ncells, hmin, None
    t = (CYL / "system" / "fvSchemes").read_text().replace("default steadyState;", "default backward;")
    (d / "system" / "fvSchemes").write_text(t)
    write_dict(d / "system" / "fvSolution", "dictionary", "fvSolution", """
solvers { p { solver PCG; preconditioner DIC; tolerance 1e-9; relTol 0.01; }
          pFinal { $p; relTol 0; }
          U { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; } UFinal { $U; } }
PIMPLE { momentumPredictor yes; nOuterCorrectors 1; nCorrectors 2; nNonOrthogonalCorrectors 0;
         consistent yes; ddtCorr yes; }
""")
    # GAMG does not converge on the cut-cell pressure equation (stalls near 1e-6): PCG
    dt = 0.8*hmin/(1.5*Um)            # Courant <= 0.8 with |u| <= 1.5 U_m
    dt = 0.33/ (np.ceil(0.33/dt))     # an integer number of steps in 0.33
    end = END
    extra = """functions { pw { type writeObjects; libs (utilityFunctionObjects); objects (p);
    writeControl adjustableRunTime; writeInterval 0.004; timeStart %r; } }""" % (end - 1.0)
    import geometry as g
    g.control_dict(d, "ibmPimpleFoam", endTime=end, deltaT=repr(float(dt)), writeInterval=int(round(end/dt)), extra=extra)
    common.set_keys(d / "system" / "controlDict", writeControl="timeStep")
    return ncells, hmin, dt


def dp_at(d, t, hmin):
    xc, a = meshes.fluid_cells(str(d), t)
    p, _ = read_field(os.path.join(d, t, "p"))
    return meshes.point_value(xc, a, p, P1, hmin, quadratic=True) - meshes.point_value(xc, a, p, P2, hmin, quadratic=True)


def steady(ctx, k):
    d = ctx.root / ("re20_k%d" % k)
    ncells, hmin, _ = setup(d, k, 0.3, False)
    common.sh("ibmSetGeometry", d, "log.ibmSetGeometry")
    common.sh("ibmSimpleFoam", d, "log.ibmSimpleFoam")
    log = (d / "log.ibmSimpleFoam").read_text()
    F = common.ibm_force(log, "cylinder")["total"]/(4*hmin)   # per unit depth: thickness h0 = 4 hmin
    U = 0.2
    t = latest_time(str(d))
    xc, a = meshes.fluid_cells(str(d), t)
    Uf, _ = read_field(os.path.join(d, t, "U"))
    # recirculation: u on the axis y = 0.2 (a mesh line: the two rows of cells beside it)
    m = (a > 0) & (np.abs(xc[:, 1] - 0.2) < hmin) & (xc[:, 0] > 0.25)
    xs = np.unique(np.round(xc[m, 0], 9))
    us = np.array([Uf[m][np.isclose(np.round(xc[m, 0], 9), x), 0].mean() for x in xs])
    i = np.where((us[:-1] < 0) & (us[1:] >= 0))[0]
    xr = xs[i[0]] - us[i[0]]*(xs[i[0] + 1] - xs[i[0]])/(us[i[0] + 1] - us[i[0]]) if len(i) else np.nan
    return {"k": k, "Dh": D/hmin, "cells": ncells, "C_D": float(2*F[0]/(U*U*D)), "C_L": float(2*F[1]/(U*U*D)),
            "dp": dp_at(d, t, hmin), "L_a": float(xr - 0.25), "converged": "SIMPLE solution converged" in log}


NP = 8   # ranks of a 2D-2 run


def periodic(ctx, k):
    d = ctx.root / ("re100_k%d" % k)
    ncells, hmin, dt = setup(d, k, 1.5, True)
    write_dict(d / "system" / "decomposeParDict", "dictionary", "decomposeParDict",
               "numberOfSubdomains %d;\nmethod scotch;\n" % NP)
    common.sh("decomposePar", d, "log.decomposePar")
    common.sh("ibmPimpleFoam", d, "log.ibmPimpleFoam", np_=NP)
    common.sh("reconstructPar -newTimes", d, "log.reconstructPar")
    common.sh("ibmSetGeometry", d, "log.ibmSetGeometry")   # geometry fields at 0 for the sampling
    f = sorted((d / "postProcessing" / "ibmForces").glob("*/forces.dat"))[-1]
    h = np.loadtxt(f)
    t, F = h[:, 0], (h[:, 1:4] + h[:, 4:7])/(4*hmin)
    cd, cl = 2*F[:, 0]/(D), 2*F[:, 1]/(D)          # U_mean = 1
    s = t > t[-1] - 3.0
    r = dfg.cycle_statistics(t[s], cd[s], cl[s])
    T = D/r["St"]
    times = sorted((x for x in os.listdir(d) if re.fullmatch(r"[0-9.e+-]+", x) and (d / x / "p").exists()
                    and float(x) > t[-1] - 1.0), key=float)
    tt = np.array([float(x) for x in times])
    j = np.argsort(np.abs(tt - (r["t0"] + T/2)))[:2]
    dps = [dp_at(d, times[i], hmin) for i in j]
    dp = float(np.interp(r["t0"] + T/2, tt[j][np.argsort(tt[j])], np.array(dps)[np.argsort(tt[j])]))
    return {"k": k, "Dh": D/hmin, "cells": ncells, "dt": dt, "St": r["St"], "C_Dmax": float(r["C_Dmax"]),
            "C_Dmin": float(r["C_Dmin"]), "C_Lmax": float(r["C_Lmax"]), "C_Lmin": float(r["C_Lmin"]), "dp": dp}


def run(ctx, res):
    ks = [0] if ctx.quick else [0, 1]
    jobs = [(steady, k) for k in ks + ([] if ctx.quick else [2])] + [(periodic, k) for k in ks]
    out = common.pmap(lambda j: (j[0].__name__, j[0](ctx, j[1])), jobs, max(1, ctx.jobs // NP))
    st = sorted([r for n, r in out if n == "steady"], key=lambda r: r["k"])
    pe = sorted([r for n, r in out if n == "periodic"], key=lambda r: r["k"])
    res.data["steady"], res.data["periodic"] = st, pe
    res.data["metrics"] = {"re20_k%d_%s" % (r["k"], q): r[q] for r in st for q in ("C_D", "C_L", "dp", "L_a")}
    res.data["metrics"].update({"re100_k%d_%s" % (r["k"], q): r[q] for r in pe
                                for q in ("St", "C_Dmax", "C_Lmax", "C_Lmin", "dp")})
    res.le("2D-1 runs not converged", sum(not r["converged"] for r in st), 0)
    f = [r for r in st if r["Dh"] == 80][0] if any(r["Dh"] == 80 for r in st) else st[-1]
    seq = lambda q, rs: " ".join("%.5g" % r[q] for r in rs)
    B = dfg.BENCH1
    res.rel("2D-1 C_D, D/h = %d" % f["Dh"], f["C_D"], B["C_D"], 0.005, seq("C_D", st))
    res.rel("2D-1 C_L, D/h = %d" % f["Dh"], f["C_L"], B["C_L"], 0.05, seq("C_L", st))
    res.rel("2D-1 dp, D/h = %d" % f["Dh"], f["dp"], B["dp"], 0.005, seq("dp", st))
    lo, hi = dfg.BENCH1_INTERVAL["L_a"]
    res.absolute("2D-1 L_a inside the 1996 interval, D/h = %d" % f["Dh"], f["L_a"], (lo + hi)/2, (hi - lo)/2,
                 seq("L_a", st))
    if len(st) >= 3:
        Ce, p = common.richardson([1/r["Dh"] for r in st], [r["C_D"] for r in st])
        res.ge("2D-1 observed order of C_D (D/h = 40, 80, 160)", p, 1.5, "extrapolated %.5f" % Ce)
        g = st[-1]
        for q, tol in (("C_D", 0.002), ("dp", 0.002)):
            res.rel("2D-1 %s, D/h = %d (tier 2 tolerance)" % (q, g["Dh"]), g[q], B[q], tol)
    R = dfg.bench2_summary()[6]
    g = pe[-1]
    res.rel("2D-2 St, D/h = %d" % g["Dh"], g["St"], R["St"], 0.005, seq("St", pe))
    for q, tol in (("C_Dmax", 0.005), ("C_Dmin", 0.005), ("C_Lmax", 0.01), ("C_Lmin", 0.01), ("dp", 0.01)):
        res.rel("2D-2 %s, D/h = %d" % (q, g["Dh"]), g[q], R[q], tol, seq(q, pe))
