"""V13 OSC: cylinder oscillating in-line in fluid at rest, Re = U_max D/nu = 100, KC = 2 pi A/D = 5
(Dütsch et al. 1998; reference_oscillatingCylinder.py; doc/validationSuite.md §10). D = 1, U_max = 1,
nu = 0.01, x(t) = -A sin(2 pi f t), A = 5/(2 pi), f = 0.2 (T = 5). Inertial frame, closed square of
half width 32 with slip walls, background h0 = 3.2 and 7, 8, 9 nested levels (D/h = 40, 80, 160)
around the stroke. Morison fit of the in-line force over cycles 5-10: c_d 2.09, c_i 1.45 (set C)."""
import math

import numpy as np

import common
import geometry
import meshes
import reference_oscillatingCylinder as osc
from common import write_dict

ID, KIND, TIER = "V13", "case", 1
TITLE = "OSC: cylinder oscillating in fluid at rest, Re 100, KC 5 (Dütsch et al. 1998)"
STRESSES, RELIES = ["K9", "K8"], ["V3"]
REGRESSION = (1e-4,)

L, H0, NU = 32.0, 3.2, 0.01
A, F = 5/(2*math.pi), 0.2
T = 1/F
NP = 8
CYL = common.TUTORIALS / "cylinderRe1"


def setup(d, levels):
    boxes = meshes.nested_boxes((-A - 1.2, -1.2, 0), (A + 1.2, 1.2, 0), H0, levels, (-L, -L, 0), (L, L, 0))
    ncells = meshes.refined(d, (-L, -L, 0), (L, L, 0), H0, boxes, {"walls": ("patch", ["xmin", "xmax", "ymin", "ymax"])},
                            thickness=0.1)
    geometry.ibm_properties(d, {"cylinder": {"type": "cylinder", "centre": (0.0123, 0.0071, 0), "axis": (0, 0, 1),
                                             "radius": 0.5, "wallType": "noSlip",
                                             "motion": {"type": "oscillating", "amplitude": (-A, 0, 0),
                                                        "frequency": F}}})
    for x in ("transportProperties", "turbulenceProperties"):
        (d / "constant" / x).write_text((CYL / "constant" / x).read_text())
    common.set_keys(d / "constant" / "transportProperties", nu=NU)
    meshes.write_fields(d, {
        "U": ("volVectorField", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)", {"walls": "type slip;"}),
        "p": ("volScalarField", "[0 2 -2 0 0 0 0]", "uniform 0", {"walls": "type slip;"})})
    (d / "system" / "fvSchemes").write_text(
        (CYL / "system" / "fvSchemes").read_text().replace("default steadyState;", "default Euler;"))
    # GAMG as the base: DIC-PCG runs into maxIter 1000 at D/h = 160 (575 iterations per solve,
    # with or without split groups). One V-cycle and the coarsest level to relTol 0.01: 0.17,
    # 0.45, 1.76 s per step at D/h = 40, 80, 160 against 0.10, 0.42 and 2.04 for DIC (maxIter
    # 5000 at D/h = 160); 32, 43, 51 iterations per solve. linkedGroup: the body splits 0 to 20
    # linked groups over the processors as it moves
    write_dict(d / "system" / "fvSolution", "dictionary", "fvSolution", """
solvers { p { solver PCG; preconditioner { preconditioner linkedGroup; base { preconditioner GAMG;
                smoother DICGaussSeidel; tolerance 0; relTol 0.01; nVcycles 1; } }
              tolerance 1e-9; relTol 0.01; }
          pFinal { $p; relTol 0; }
          U { solver PBiCGStab; preconditioner { preconditioner linkedGroup; base DILU; } tolerance 1e-11; relTol 0; }
          UFinal { $U; } }
PIMPLE { momentumPredictor yes; nOuterCorrectors 2; nCorrectors 2; nNonOrthogonalCorrectors 0;
         consistent yes; ddtCorr yes; pRefPoint (20 20 0.05); pRefValue 0; }
""")
    h = H0/2**levels
    n = max(720, int(math.ceil(T/(0.8*h/1.2))))   # T/720 (the paper) or Courant 0.8
    dt = T/n
    geometry.control_dict(d, "ibmPimpleFoam", endTime=10*T, deltaT=repr(float(dt)), writeInterval=10*n)
    write_dict(d / "system" / "decomposeParDict", "dictionary", "decomposeParDict",
               "numberOfSubdomains %d;\nmethod scotch;\n" % NP)
    return ncells, h, dt


def one(ctx, levels):
    d = ctx.root / ("Dh%d" % round(1/(H0/2**levels)))
    if not (common.REUSE and (d / "system").exists()):
        ncells, h, dt = setup(d, levels)
    else:
        ncells, h, dt = None, H0/2**levels, None
    common.sh("decomposePar", d, "log.decomposePar")
    common.sh("ibmPimpleFoam", d, "log.ibmPimpleFoam", np_=NP)
    t, Fv = meshes.force_history(d, 0.1)
    s = t >= 5*T - 1e-9
    cd, ci = osc.morison_fit(t[s], Fv[s, 0], F, A)
    return {"Dh": round(1/h), "cells": ncells, "dt": dt, "c_d": float(cd), "c_i": float(ci)}


def run(ctx, res):
    lv = (7, 8) if ctx.quick else (7, 8, 9)
    rows = common.pmap(lambda j: one(ctx, j), lv, max(1, ctx.jobs // NP))
    rows.sort(key=lambda r: r["Dh"])
    res.data["rows"] = rows
    res.data["metrics"] = {"Dh%d_%s" % (r["Dh"], q): r[q] for r in rows for q in ("c_d", "c_i")}
    R = osc.REFERENCE
    g = rows[-1]
    for q in ("c_d", "c_i"):
        res.rel("%s, D/h = %d (Morison fit, cycles 5-10)" % (q, g["Dh"]), g[q], R[q], 0.02,
                " ".join("%.4f" % r[q] for r in rows))
    if len(rows) >= 3:
        for q in ("c_d", "c_i"):
            e, p = common.richardson([1/r["Dh"] for r in rows], [r[q] for r in rows])
            res.report("%s Richardson" % q, e, R[q], "observed order %.2f" % p)
