"""V12 ROT: cylinder rotating in a uniform stream, Re = 200, alpha = a omega/U (Mittal & Kumar
2003; reference_rotatingCylinder.py; doc/validationSuite.md §10). D = 1, U = 1, nu = 1/200, square
domain of half width 102.4 (their L = 100 D), inflow U, outflow p = 0, slip sides. Background
h0 = 3.2 with 7 (D/h = 40) or 8 (D/h = 80) nested levels. The rotation is a tangential wall
velocity of an invariant geometry (motion none, angularVelocity 2 alpha), ibmPimpleFoam; statistics
over the whole cycles in the last 50 time units of 150. Tier 1: alpha = 0, 1, 2; alpha = 2 also at
D/h = 160 (9 levels, 6 NP ranks; decision D4: -2.2 % at D/h = 80, converging). Torque C_M = M/(U^2 a^2)
(positive against the rotation; postProcessing/ibmForces moments.dat, C1) against the body-fitted
solution, figure 11 reported."""
import numpy as np

import common
import geometry
import meshes
import reference_rotatingCylinder as rot
from common import write_dict

ID, KIND, TIER = "V12", "case", 1
TITLE = "ROT: rotating cylinder in a stream, Re 200 (Mittal & Kumar 2003)"
STRESSES, RELIES = ["K3", "K7", "K8"], ["V2", "V5"]
RTOL = 1e-4


def load_atol(b):
    """absolute tolerance of the force and moment coefficients: RTOL times the load scale of their
    run, max |C_D|, |C_L|, |C_M| of the baseline (like G7). The alpha = 0 mean C_L and C_M are zero
    by symmetry; a relative tolerance on them measures round-off. St keeps the relative tolerance."""
    scale = {}
    for k, v in b.items():
        if v is not None and k.endswith(("_C_L", "_C_D", "_C_M")):
            run = k[:k.index("_C_")]
            scale[run] = max(scale.get(run, 0.0), abs(v))
    return {k: RTOL*scale.get(k[:k.index("_C_")], 0.0) for k in b if "_C_" in k}


REGRESSION = (RTOL, None, load_atol)

L, H0, NU = 102.4, 3.2, 1/200.0
CENTRE = (0.0123, 0.0071, 0.0)   # off the mesh lines
END, WINDOW = 150.0, 50.0
SHEDDING = 1.9   # alpha below which the wake sheds (Mittal & Kumar 2003): St defined
NP = 8
CYL = common.TUTORIALS / "cylinderRe1"


def ranks(levels):
    return NP if levels < 9 else 6*NP


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
solvers { p { solver PCG; preconditioner { preconditioner linkedGroup; base DIC; } tolerance 1e-8; relTol 0.01; }
          pFinal { $p; relTol 0; }
          U { solver PBiCGStab; preconditioner { preconditioner linkedGroup; base DILU; } tolerance 1e-10; relTol 0; }
          UFinal { $U; } }
PIMPLE { momentumPredictor yes; nOuterCorrectors 3; nCorrectors 2; nNonOrthogonalCorrectors 0;
         consistent yes; ddtCorr yes; }
""")
    # 3 outer correctors: one PISO pass per step (convection with the old flux) is first order in
    # time: DFG 2D-2, D/h = 40, C_L amplitude +12.8 %, +6.7 % at dt/2, +0.8 % with 3 (doc §12)
    h = H0/2**levels
    dt = 0.5*h/(1.5 + alpha)                  # Courant 0.5 with |u| <= 1.5 + alpha
    geometry.control_dict(d, "ibmPimpleFoam", endTime=END, deltaT=repr(float(dt)), writeInterval=int(round(END/dt)))
    write_dict(d / "system" / "decomposeParDict", "dictionary", "decomposeParDict",
               "numberOfSubdomains %d;\nmethod scotch;\n%s" % (ranks(levels), common.LINK_CONSTRAINT))
    return ncells, h


def one(ctx, alpha, levels):
    d = ctx.root / ("alpha%g_Dh%d" % (alpha, round(1/(H0/2**levels))))
    if not (common.REUSE and (d / "system").exists()):
        ncells, h = setup(d, alpha, levels)
    else:
        ncells, h = None, H0/2**levels
    common.link_zone(d)      # the cylinder rotates in place: the links do not change
    common.sh("decomposePar", d, "log.decomposePar")
    common.sh("ibmPimpleFoam", d, "log.ibmPimpleFoam", np_=ranks(levels))
    t, F = meshes.force_history(d, 0.1)
    tm, M = meshes.moment_history(d, 0.1)
    assert np.array_equal(t, tm), "forces.dat and moments.dat times differ"
    st = meshes.periodic_statistics(t, 2*F[:, 1], 2*F[:, 0], WINDOW, -M[:, 2]/0.5**2)
    return {"alpha": alpha, "Dh": round(1/h), "cells": ncells, "St": st.pop("f"), **st}


def run(ctx, res):
    jobs = [(a, lv) for a in (0, 1, 2) for lv in ((7,) if ctx.quick else (7, 8))]
    jobs = ([] if ctx.quick else [(2, 9)]) + jobs      # longest first
    rows = common.pmap(lambda j: one(ctx, *j), jobs, max(1, ctx.jobs // NP))
    res.data["rows"] = rows
    # no St without shedding: the dominant frequency of a decaying lift is not a quantity to regress
    res.data["metrics"] = {"a%g_Dh%d_%s" % (r["alpha"], r["Dh"], q): r[q] for r in rows
                           for q in ("C_L", "C_D", "C_Lrms", "St", "C_M") if r.get(q) is not None and np.isfinite(r[q])
                           and (q != "St" or r["alpha"] < SHEDDING)}
    ref = rot.figure11()
    bf = rot.bodyFitted()
    for a in (0, 1, 2):
        s = sorted([r for r in rows if r["alpha"] == a], key=lambda r: r["Dh"])
        g = s[-1]
        i = int(np.argmin(np.abs(ref["alpha"] - a)))
        seq = lambda q: " ".join("%.4g" % r[q] for r in s)
        for q in ("C_L", "C_D"):
            R = ref[q][i]
            tol = max(0.02*abs(R), 0.02)
            if a == 1 and q == "C_L":
                # the digitised figure 11 (-2.265) is 10 % off the body-fitted solution of the same
                # problem (-2.499, mesh-converged), to which the IBM converges (finding 7, D3)
                res.report("alpha 1: mean C_L, D/h = %d, against figure 11 (digitised)" % g["Dh"], g[q],
                           ref=R, note=seq(q))
                R = bf[a][q]
                res.absolute("alpha 1: mean C_L, D/h = %d, against body-fitted (D3)" % g["Dh"], g[q], R,
                             max(0.02*abs(R), 0.02), seq(q))
                continue
            res.absolute("alpha %g: mean %s, D/h = %d" % (a, q, g["Dh"]), g[q], R, tol, seq(q))
        for q in ("C_L", "C_D", "C_Lrms", "St"):
            res.report("alpha %g: %s, D/h = %d, vs body-fitted M 128 (rel.; abs. where |ref| < 0.01)" % (a, q, g["Dh"]),
                       g[q]/bf[a][q] - 1 if abs(bf[a][q]) > 0.01 else g[q] - bf[a][q], note=seq(q))
        for q in ("C_Lrms", "C_Drms"):
            R = ref[q][i]
            res.absolute("alpha %g: %s, D/h = %d" % (a, q, g["Dh"]), g[q], R, max(0.05*abs(R), 0.01), seq(q))
        if a < SHEDDING:
            R = rot.TEXT["alpha=0"]["St"] if a == 0 else ref["St"][i]
            res.rel("alpha %g: St, D/h = %d" % (a, g["Dh"]), g["St"], R, 0.01, seq("St"))
        # torque: body-fitted (OpenFOAM forces, full stress) mesh-converged to 4 digits; figure 11
        # gives it only as (C_P - C_D)/alpha, uncertain by ~0.03/alpha
        R = bf[a]["C_M"]
        res.absolute("alpha %g: mean C_M, D/h = %d, against body-fitted" % (a, g["Dh"]), g["C_M"], R,
                     max(0.02*abs(R), 0.005), seq("C_M"))
        if a >= 1:
            res.report("alpha %g: mean C_M, D/h = %d, against figure 11 ((C_P - C_D)/alpha)" % (a, g["Dh"]),
                       g["C_M"], ref=ref["C_M"][i], note=seq("C_M"))
