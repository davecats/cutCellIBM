"""V12 ROT: cylinder rotating in a uniform stream, Re = 200, alpha = a omega/U (Mittal & Kumar
2003; reference_rotatingCylinder.py; doc/validationSuite.md §10). D = 1, U = 1, nu = 1/200, square
domain of half width 102.4 (their L = 100 D), inflow U, outflow p = 0, slip sides. Background
h0 = 3.2 with 7 (D/h = 40) or 8 (D/h = 80) nested levels. The rotation is a tangential wall
velocity of an invariant geometry (motion none, angularVelocity 2 alpha), ibmPimpleFoam; statistics
over the whole cycles in the last 50 time units of 150. Tier 1: alpha = 0, 1, 2; alpha = 2 also at
D/h = 160 (9 levels, 4 NP ranks; decision D4: -2.2 % at D/h = 80, converging); torque (C_M) needs C1."""
import numpy as np

import common
import geometry
import meshes
import reference_rotatingCylinder as rot
from common import write_dict

ID, KIND, TIER = "V12", "case", 1
TITLE = "ROT: rotating cylinder in a stream, Re 200 (Mittal & Kumar 2003)"
STRESSES, RELIES = ["K3", "K7", "K8"], ["V2", "V5"]
REGRESSION = (1e-4,)

L, H0, NU = 102.4, 3.2, 1/200.0
CENTRE = (0.0123, 0.0071, 0.0)   # off the mesh lines
END, WINDOW = 150.0, 50.0
NP = 8
CYL = common.TUTORIALS / "cylinderRe1"


def ranks(levels):
    return NP if levels < 9 else 4*NP


def setup(d, alpha, levels):
    boxes = meshes.nested_boxes((-0.8, -0.8, 0), (3.0, 0.8, 0), H0, levels, (-L, -L, 0), (L, L, 0), wake=2.0)
    patches = {"inlet": ("patch", ["xmin"]), "outlet": ("patch", ["xmax"]), "sides": ("patch", ["ymin", "ymax"])}
    ncells = meshes.refined(d, (-L, -L, 0), (L, L, 0), H0, boxes, patches, thickness=0.1)
    geometry.ibm_properties(d, {"cylinder": {"type": "cylinder", "centre": CENTRE, "axis": (0, 0, 1), "radius": 0.5,
                                             "wallType": "noSlip",
                                             "motion": {"type": "none", "angularVelocity": (0, 0, 2.0*alpha),
                                                        "centreOfRotation": CENTRE}}})
    for x in ("transportProperties", "turbulenceProperties"):
        (d / "constant" / x).write_text((CYL / "constant" / x).read_text())
    common.set_keys(d / "constant" / "transportProperties", nu=NU)
    meshes.write_fields(d, {
        "U": ("volVectorField", "[0 1 -1 0 0 0 0]", "uniform (1 0 0)",
              {"inlet": "type fixedValue; value uniform (1 0 0);", "outlet": "type zeroGradient;", "sides": "type slip;"}),
        "p": ("volScalarField", "[0 2 -2 0 0 0 0]", "uniform 0",
              {"inlet": "type zeroGradient;", "outlet": "type fixedValue; value uniform 0;", "sides": "type slip;"})})
    t = (CYL / "system" / "fvSchemes").read_text().replace("default steadyState;", "default backward;")
    (d / "system" / "fvSchemes").write_text(t)
    write_dict(d / "system" / "fvSolution", "dictionary", "fvSolution", """
solvers { p { solver PCG; preconditioner DIC; tolerance 1e-8; relTol 0.01; } pFinal { $p; relTol 0; }
          U { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; } UFinal { $U; } }
PIMPLE { momentumPredictor yes; nOuterCorrectors 3; nCorrectors 2; nNonOrthogonalCorrectors 0;
         consistent yes; ddtCorr yes; }
""")
    # 3 outer correctors: one PISO pass per step (convection with the old flux) is first order in
    # time: DFG 2D-2, D/h = 40, C_L amplitude +12.8 %, +6.7 % at dt/2, +0.8 % with 3 (doc §12)
    h = H0/2**levels
    dt = 0.5*h/(1.5 + alpha)                  # Courant 0.5 with |u| <= 1.5 + alpha
    geometry.control_dict(d, "ibmPimpleFoam", endTime=END, deltaT=repr(float(dt)), writeInterval=int(round(END/dt)))
    write_dict(d / "system" / "decomposeParDict", "dictionary", "decomposeParDict",
               "numberOfSubdomains %d;\nmethod scotch;\n" % ranks(levels))
    return ncells, h


def one(ctx, alpha, levels):
    d = ctx.root / ("alpha%g_Dh%d" % (alpha, round(1/(H0/2**levels))))
    if not (common.REUSE and (d / "system").exists()):
        ncells, h = setup(d, alpha, levels)
    else:
        ncells, h = None, H0/2**levels
    common.sh("decomposePar", d, "log.decomposePar")
    common.sh("ibmPimpleFoam", d, "log.ibmPimpleFoam", np_=ranks(levels))
    t, F = meshes.force_history(d, 0.1)
    st = meshes.periodic_statistics(t, 2*F[:, 1], 2*F[:, 0], WINDOW)
    return {"alpha": alpha, "Dh": round(1/h), "cells": ncells, "St": st.pop("f"), **st}


def run(ctx, res):
    jobs = [(a, lv) for a in (0, 1, 2) for lv in ((7,) if ctx.quick else (7, 8))]
    jobs = ([] if ctx.quick else [(2, 9)]) + jobs      # longest first
    rows = common.pmap(lambda j: one(ctx, *j), jobs, max(1, ctx.jobs // NP))
    res.data["rows"] = rows
    res.data["metrics"] = {"a%g_Dh%d_%s" % (r["alpha"], r["Dh"], q): r[q] for r in rows
                           for q in ("C_L", "C_D", "C_Lrms", "St") if r[q] is not None and np.isfinite(r[q])}
    ref = rot.figure11()
    for a in (0, 1, 2):
        s = sorted([r for r in rows if r["alpha"] == a], key=lambda r: r["Dh"])
        g = s[-1]
        i = int(np.argmin(np.abs(ref["alpha"] - a)))
        seq = lambda q: " ".join("%.4g" % r[q] for r in s)
        for q in ("C_L", "C_D"):
            R = ref[q][i]
            res.absolute("alpha %g: mean %s, D/h = %d" % (a, q, g["Dh"]), g[q], R, max(0.02*abs(R), 0.02), seq(q))
        for q in ("C_Lrms", "C_Drms"):
            R = ref[q][i]
            res.absolute("alpha %g: %s, D/h = %d" % (a, q, g["Dh"]), g[q], R, max(0.05*abs(R), 0.01), seq(q))
        if a < 1.9:
            R = rot.TEXT["alpha=0"]["St"] if a == 0 else ref["St"][i]
            res.rel("alpha %g: St, D/h = %d" % (a, g["Dh"]), g["St"], R, 0.01, seq("St"))
    res.report("torque C_M", None, note="needs the torque output (C1)")
