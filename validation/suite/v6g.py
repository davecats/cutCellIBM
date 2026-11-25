"""V6G: V6 MIG in the laboratory frame (Galilean invariance). The pipe is still and the sphere moves
freely along the axis and turns about z under the force and torque of the flow (motion type free,
src/cutCellIBM/ibmBody/freeMotion.H), its transverse position held, as in V6. In the frame of the
sphere the flow is the steady V6 flow, so the time-averaged state of a free sphere that has
settled must reproduce the V6 balance (V, Omega, f/Re, beta) at the same point and mesh; the
difference is the error of the moving-body discretisation (the sphere crosses the cells) and of
the free-body coupling. Only the IBM runs of V6 are cloned, at a few points.

Gate V6G-G1, free sphere in the periodic Stokes box (V6 MIG-G2, rigid, L/d = 4, d/h = 16, Re = 0.01):
a constant external force F_ext along x on the sphere, zero total flux (fluid and sphere,
ibmMeanVelocityForce includeMovingBodies), mean pressure gradient on the sphere volume. The
terminal velocity U_t must satisfy F_ext = K 6 pi mu R U_t with K the array drag of the static
sphere on the same mesh (the G2 configuration: in the frame of the sphere the superficial velocity
is U_t). F_ext is set to the drag of the static run at U = 1, so U_t = 1. Runs: implicit coupling
at rho_s/rho = 1, 0.5, 0.1 (below the added mass) and explicit at 10 (explicit coupling in creeping
flow needs a heavy body: the unsteady viscous load of a step acts as an added mass growing like
dt^-1/2).

Cases: rigid d = 0.4 at d/h = 32 (V6 h80 meshes, refined over the whole length because the sphere
travels): Re 2 at eps/eps* = 0.4 and 0.9, Re 128 at 0.4; neutrally buoyant (rho_s/rho = 1,
implicit); at Re 2, 0.4 also rho_s/rho = 0.5 implicit and 2 explicit (the settled state must not
depend on the density or the coupling). The flow starts as Poiseuille flow, the sphere with the
velocity and half the vorticity of the undisturbed flow at its centre. Averages over the last
AVG of the run; settled if the averages of its two halves agree to SETTLE (V, Omega) and the
force balance holds.

Units rho = J = d_h = 1 (V6). Forces and volumes of the half domain (z >= 0) are doubled.
Tolerances (provisional, to be decided): gate U_t within 1 %; V and Omega within 0.5 %, f/Re within
0.01 absolute and beta within 2 % of V6 at the same point and mesh.
"""
import json
import math
import re
import shutil
from pathlib import Path

import numpy as np

import common
import geometry
import meshes
import mig
import v6
import reference_riveroScheid2018 as ref
from common import write_dict
from foamio import read_field

ID, KIND, TIER = "V6G", "case", 2
TITLE = "MIG in the laboratory frame: free sphere in a still pipe (Galilean check of V6)"
STRESSES, RELIES = ["K9", "K8", "K10", "K1", "K3", "K7"], ["V6", "V13", "G4"]
REGRESSION = None              # no baseline yet

G_L, G_DH, G_RE = 4.0, 16, 0.01
GATE_DT = 5e-4                 # the checks: U_t converges ~dt^2 to within 2e-4 of the static drag
GATE_RUNS = [("implicit", 1.0), ("implicit", 0.5), ("implicit", 0.1), ("explicit", 10.0)]
GATE_DT_SCAN = [5e-3, 2.5e-3, 1e-3]   # implicit rho_s/rho 1, REPORT (first runs, 2026-11-24:
                                      # +1.33, +0.32, +0.02 %, -0.02 % at 5e-4)
GATE_TOL = 0.01

D, DH = 0.4, 32
POINTS = [(2.0, 0.4), (2.0, 0.9), (128.0, 0.4)]
VARIANTS = [("implicit", 1.0), ("implicit", 0.5), ("explicit", 2.0)]   # at the first point
AVG = 1/3                      # averaging window: the last third of the run
SETTLE = 2e-3                  # V, Omega: the two halves of the window agree to this (relative)
TOL = {"V": 0.005, "Om": 0.005, "fRe": 0.01, "beta": 0.02}
NOUTER = 3                     # PIMPLE outer correctors = coupling iterations (implicit)
NOUTER_CAP = 60                # cap of the outer correctors with p residual control (ptol)
NP_GATE = 4
NP_CASE = 32                   # ~13k cells per rank: the geometry is rebuilt NOUTER times a step


# ----------------------------------------------------------------------------------- common

def free_lines(log, body="sphere"):
    """per time step (the 'final' line of the step): t, velocity, angularVelocity, force, moment"""
    v = r"\(([^)]*)\)"
    pat = (r"IBM free body %s final: t (\S+) displacement %s velocity %s angularVelocity %s force %s "
           r"moment %s" % (re.escape(body), v, v, v, v, v))
    rows = []
    for m in re.findall(pat, log):
        rows.append([float(m[0])] + [float(x) for g in m[1:] for x in g.split()])
    a = np.array(rows)
    return {"t": a[:, 0], "d": a[:, 1:4], "U": a[:, 4:7], "Om": a[:, 7:10], "F": a[:, 10:13], "M": a[:, 13:16]}


def gradients(log):
    g, _ = common.driving_gradient(log)
    return g


def window(t, q, frac=AVG):
    """mean over the last frac of the run, and the relative difference of the means of its halves"""
    t0 = t[-1] - frac*(t[-1] - t[0])
    m = t >= t0
    q = np.asarray(q)[m]
    h = len(q)//2
    a, b = q[:h].mean(axis=0), q[h:].mean(axis=0)
    return q.mean(axis=0), (b - a), q.std(axis=0)


def pimple_solution(pref, nouter, tight=True, ptol=None):
    """ptol: outer-loop residual control on p (PIMPLE residualControl, the initial residual of the
    latest p solve below ptol, then the final corrector; nouter the cap). The U residual is not used:
    OpenFOAM takes the largest component, the small cross-flow one in these cases.
    tight (the gate): U to 1e-12 in every outer corrector; else (the cases) relTol 0.01 in the
    non-final correctors and 1e-9 in the final one (U to 1e-12 costs 100-400 PBiCGStab iterations
    per component and corrector on the d/h 32 pipe)"""
    tu, ru, tp = (1e-12, 0, 1e-10) if tight else (1e-9, 0.01, 1e-9)
    return """
solvers { p { solver PCG; preconditioner { preconditioner linkedGroup; base { preconditioner GAMG;
                smoother DICGaussSeidel; tolerance 0; relTol 0.01; } }
              tolerance %g; relTol 0.01; }
          pFinal { $p; relTol 0; }
          U { solver PBiCGStab; preconditioner { preconditioner linkedGroup; base DILU; } tolerance %g; relTol %g; }
          UFinal { $U; relTol 0; } }
PIMPLE { momentumPredictor yes; nOuterCorrectors %d; nCorrectors 2; nNonOrthogonalCorrectors 0;
         consistent yes; ddtCorr yes; pRefPoint (%r %r %r); pRefValue 0;%s }
""" % ((tp, tu, ru, nouter) + tuple(pref) +
       (" residualControl { p { tolerance %g; relTol 0; } }" % ptol if ptol else "",))


def run_parallel(d, app, np_):
    write_dict(d / "system" / "decomposeParDict", "dictionary", "decomposeParDict",
               "numberOfSubdomains %d;\nmethod scotch;\n" % np_)
    common.sh("decomposePar -force", d, "log.decomposePar")
    common.sh(app, d, "log." + app, np_=np_)


# ------------------------------------------------------------------------------------- gate

def gate_mesh(d):
    """the G2 mesh (v6.g2_mesh) with the nested boxes over the whole length in x"""
    hf = 1.0/G_DH
    n = int(round(np.log2(G_DH/4.0)))
    h0 = hf*2**n
    boxes = []
    for k in range(1, n + 1):
        w = 0.5 + 8*h0/2**k
        boxes.append(((-G_L, -w, -1.0), (G_L, w, w)))
    patches = {"left": ("cyclic", ["xmin"], "neighbourPatch right;"),
               "right": ("cyclic", ["xmax"], "neighbourPatch left;"),
               "front": ("cyclic", ["ymin"], "neighbourPatch back;"),
               "back": ("cyclic", ["ymax"], "neighbourPatch front;"),
               "sym": ("symmetryPlane", ["zmin"]), "top": ("symmetryPlane", ["zmax"])}
    return meshes.refined(d, (-G_L/2, -G_L/2, 0.0), (G_L/2, G_L/2, G_L/2), h0, boxes, patches, two_d=False)


def gate_static(ctx):
    """the static sphere of G2 (rigid, U = 1) on the gate mesh: K and the drag at U = 1"""
    d = ctx.root / "g1_static"
    f = d / "result.json"
    if not f.exists():
        d.mkdir(parents=True, exist_ok=True)
        if not (d / "constant" / "polyMesh" / "owner").exists():
            gate_mesh(d)
        v6.g2_setup(d, G_L, G_DH, "rigid")
        Vf = mig.geometry_only(d)
        r = mig.solve(d)
        r["Vf"] = 2*Vf
        f.write_text(json.dumps(r))
    r = json.loads(f.read_text())
    mu, R = 1.0/G_RE, 0.5
    U = r["Vf"]/G_L**3                                   # superficial velocity
    r["drag"] = 2*r["F"][0] + r["G"]*4*np.pi/3*R**3      # force of the array on one sphere
    r["K"] = r["drag"]/(6*np.pi*mu*R*U)
    r["F1"] = r["K"]*6*np.pi*mu*R                        # the drag at U = 1
    return r


def gate_free(ctx, coupling, rho, st, dt=None):
    d = ctx.root / ("g1_free_%s_rho%g%s" % (coupling, rho, "_dt%g" % dt if dt else ""))
    f = d / "result.json"
    if f.exists():
        return json.loads(f.read_text())
    d.mkdir(parents=True, exist_ok=True)
    src = ctx.root / "g1_static"
    if not (d / "constant" / "polyMesh").exists():
        shutil.copytree(src / "constant" / "polyMesh", d / "constant" / "polyMesh")
        shutil.copytree(src / "system", d / "system")
        shutil.copy(src / "log.checkMesh", d / "log.checkMesh")
    R = 0.5
    Vh = 2*np.pi/3*R**3                                  # half sphere in the domain
    s = {"type": "sphere", "centre": (0, 0, 0), "radius": R, "wallType": "noSlip",
         "motion": {"type": "free", "density": rho, "translation": (1, 0, 0), "rotation": (0, 0, 0),
                    "coupling": coupling, "addedMass": 0.5, "volume": Vh,
                    "externalForce": (st["F1"]/2, 0, 0)}}
    geometry.ibm_properties(d, {"sphere": s}, periodicLengths=(G_L, G_L, 0))
    for x in ("transportProperties", "turbulenceProperties"):
        (d / "constant" / x).write_text((mig.CYL / "constant" / x).read_text())
    common.set_keys(d / "constant" / "transportProperties", nu=1.0/G_RE)
    write_dict(d / "constant" / "fvOptions", "dictionary", "fvOptions", """
momentumSource { type ibmMeanVelocityForce; selectionMode all; fields (U); Ubar (0 0 0); flowDir (1 0 0);
                 relaxation 1.0; includeMovingBodies yes; }
""")
    bc = ("left { type cyclic; } right { type cyclic; } front { type cyclic; } back { type cyclic; } "
          "sym { type symmetryPlane; } top { type symmetryPlane; }")
    for name, cls, dims, val in (("U", "volVectorField", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)"),
                                 ("p", "volScalarField", "[0 2 -2 0 0 0 0]", "uniform 0")):
        write_dict(d / "0" / name, cls, name, "dimensions %s;\ninternalField %s;\nboundaryField { %s }\n"
                   % (dims, val, bc))
    (d / "system" / "fvSchemes").write_text(
        (mig.CYL / "system" / "fvSchemes").read_text().replace("default steadyState;", "default Euler;"))
    write_dict(d / "system" / "fvSolution", "dictionary", "fvSolution",
               pimple_solution((0.45*G_L, 0.45*G_L, 0.2*G_L), NOUTER))
    # body relaxation time (m + M_a)/D ~ 5e-4 at rho 1; explicit: dt D/m < 2
    dt0, T = (0.005, 1.0) if coupling == "implicit" else (0.0005, 0.25)
    dt = dt or dt0
    T = min(T, 500*dt) if dt < dt0 else T
    geometry.control_dict(d, "ibmPimpleFoam", endTime=T, deltaT=dt, writeInterval=int(round(T/dt)))
    run_parallel(d, "ibmPimpleFoam", NP_GATE)
    h = free_lines((d / "log.ibmPimpleFoam").read_text())
    U, dU, sU = window(h["t"], h["U"][:, 0], 0.5)
    r = {"coupling": coupling, "rho": rho, "dt": dt, "Ut": float(U), "dUt": float(dU), "sUt": float(sU),
         "Fbal": float(np.abs(h["F"][-20:, 0]).max()/st["F1"]*2), "steps": len(h["t"])}
    f.write_text(json.dumps(r))
    return r


def gate_prescribed(ctx, Om, dt=5e-4, T=0.1):
    """diagnostic: the sphere of the gate box translating at U = 1 along x and turning at Om about z
    (prescribed, moving geometry), zero total flux, Stokes: the transverse force vanishes exactly
    (linearity, the symmetry of the box), so any F_y is an error of the moving discretisation"""
    d = ctx.root / ("g1_prescribed_Om%g_dt%g" % (Om, dt))
    f = d / "result.json"
    if f.exists():
        return json.loads(f.read_text())
    d.mkdir(parents=True, exist_ok=True)
    src = ctx.root / "g1_static"
    if not (d / "constant" / "polyMesh").exists():
        shutil.copytree(src / "constant" / "polyMesh", d / "constant" / "polyMesh")
        shutil.copytree(src / "system", d / "system")
        shutil.copy(src / "log.checkMesh", d / "log.checkMesh")
    s = {"type": "sphere", "centre": (0, 0, 0), "radius": 0.5, "wallType": "noSlip",
         "motion": {"type": "translating", "velocity": (1.0, 0, 0), "angularVelocity": (0, 0, Om)}}
    geometry.ibm_properties(d, {"sphere": s}, periodicLengths=(G_L, G_L, 0))
    for x in ("transportProperties", "turbulenceProperties"):
        (d / "constant" / x).write_text((mig.CYL / "constant" / x).read_text())
    common.set_keys(d / "constant" / "transportProperties", nu=1.0/G_RE)
    write_dict(d / "constant" / "fvOptions", "dictionary", "fvOptions", """
momentumSource { type ibmMeanVelocityForce; selectionMode all; fields (U); Ubar (0 0 0); flowDir (1 0 0);
                 relaxation 1.0; includeMovingBodies yes; }
""")
    bc = ("left { type cyclic; } right { type cyclic; } front { type cyclic; } back { type cyclic; } "
          "sym { type symmetryPlane; } top { type symmetryPlane; }")
    for name, cls, dims, val in (("U", "volVectorField", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)"),
                                 ("p", "volScalarField", "[0 2 -2 0 0 0 0]", "uniform 0")):
        write_dict(d / "0" / name, cls, name, "dimensions %s;\ninternalField %s;\nboundaryField { %s }\n"
                   % (dims, val, bc))
    (d / "system" / "fvSchemes").write_text(
        (mig.CYL / "system" / "fvSchemes").read_text().replace("default steadyState;", "default Euler;"))
    write_dict(d / "system" / "fvSolution", "dictionary", "fvSolution",
               pimple_solution((0.45*G_L, 0.45*G_L, 0.2*G_L), NOUTER))
    geometry.control_dict(d, "ibmPimpleFoam", endTime=T, deltaT=dt, writeInterval=int(round(T/dt)))
    run_parallel(d, "ibmPimpleFoam", NP_GATE)
    log = (d / "log.ibmPimpleFoam").read_text()
    m = re.findall(r"IBM force on sphere: pressure \(([^)]*)\) viscous \(([^)]*)\)", log)
    a = np.array([[float(x) for x in (p + " " + v).split()] for p, v in m])[len(m)//2:]
    r = {"Om": Om, "dt": dt, "Fp": a[:, :3].mean(axis=0).tolist(), "Fv": a[:, 3:].mean(axis=0).tolist(),
         "Fy_abs": float(np.abs(a[:, 1] + a[:, 4]).mean())}
    f.write_text(json.dumps(r))
    return r


def gate_galilean(ctx, Re, Om, dt, frame, T=2.0, np_=4, tag=""):
    """diagnostic: the sphere of the gate box at Re = U D/nu with inertia, turning at Om about z,
    either translating at U = 1 through fluid at rest (frame "lab", moving geometry) or held still
    in fluid moving at -1 (frame "body", static geometry, invariant rotation), both started
    impulsively and without a mean-flow source (the fluid exchanges momentum with the sphere
    only; ibmMeanVelocityForce is unstable in this wall-less box at Re 10, growing ~2x per step
    in either frame). The two are Galilean equivalent at every instant, so the difference of the
    force histories is the error of the moving-geometry discretisation alone (lift by the spin)"""
    d = ctx.root / ("g1_galilean_%s_Re%g_Om%g_dt%g_T%g%s" % (frame, Re, Om, dt, T, tag))
    f = d / "result.json"
    if f.exists():
        return json.loads(f.read_text())
    gate_static(ctx)
    d.mkdir(parents=True, exist_ok=True)
    src = ctx.root / "g1_static"
    if not (d / "constant" / "polyMesh").exists():
        shutil.copytree(src / "constant" / "polyMesh", d / "constant" / "polyMesh")
        shutil.copytree(src / "system", d / "system")
        shutil.copy(src / "log.checkMesh", d / "log.checkMesh")
    s = {"type": "sphere", "centre": (0, 0, 0), "radius": 0.5, "wallType": "noSlip"}
    if frame == "lab":
        s["motion"] = {"type": "translating", "velocity": (1.0, 0, 0), "angularVelocity": (0, 0, Om)}
        u0 = 0.0
    else:
        s["motion"] = {"type": "none", "angularVelocity": (0, 0, Om), "invariant": True}
        u0 = -1.0
    geometry.ibm_properties(d, {"sphere": s}, periodicLengths=(G_L, G_L, 0))
    for x in ("transportProperties", "turbulenceProperties"):
        (d / "constant" / x).write_text((mig.CYL / "constant" / x).read_text())
    common.set_keys(d / "constant" / "transportProperties", nu=1.0/Re)
    bc = ("left { type cyclic; } right { type cyclic; } front { type cyclic; } back { type cyclic; } "
          "sym { type symmetryPlane; } top { type symmetryPlane; }")
    for name, cls, dims, val in (("U", "volVectorField", "[0 1 -1 0 0 0 0]", "uniform (%r 0 0)" % u0),
                                 ("p", "volScalarField", "[0 2 -2 0 0 0 0]", "uniform 0")):
        write_dict(d / "0" / name, cls, name, "dimensions %s;\ninternalField %s;\nboundaryField { %s }\n"
                   % (dims, val, bc))
    (d / "system" / "fvSchemes").write_text(
        (mig.CYL / "system" / "fvSchemes").read_text().replace("default steadyState;", "default Euler;"))
    write_dict(d / "system" / "fvSolution", "dictionary", "fvSolution",
               pimple_solution((0.45*G_L, 0.45*G_L, 0.2*G_L), NOUTER))
    n = int(round(T/dt))
    geometry.control_dict(d, "ibmPimpleFoam", endTime=T, deltaT=dt, writeInterval=n)
    run_parallel(d, "ibmPimpleFoam", np_)
    log = (d / "log.ibmPimpleFoam").read_text()
    t = [float(x) for x in re.findall(r"^Time = (\S+)", log, re.M)]
    m = re.findall(r"IBM force on sphere: pressure \(([^)]*)\) viscous \(([^)]*)\)", log)
    a = np.array([[float(x) for x in (p + " " + v).split()] for p, v in m])
    t = np.array(t[:len(a)])
    np.save(d / "forces.npy", np.column_stack([t, a]))
    w = t > T/2
    F = 2*(a[:, :3] + a[:, 3:])                    # both halves of the box (symmetry plane z = 0)
    r = {"frame": frame, "Re": Re, "Om": Om, "dt": dt, "T": T, "steps": len(a),
         "F": F[w].mean(axis=0).tolist(), "Fsd": F[w].std(axis=0).tolist()}
    f.write_text(json.dumps(r))
    return r


# ------------------------------------------------------------------------------------ cases

# ------------------------------------------------- six degrees of freedom, several bodies

C_L, C_H0, C_LEVELS, C_RHO = 4.0, 0.25, 2, 1000.0   # cube side, background cell, levels (h 1/16), rho_s/rho
C_T, C_DT = 0.3, 0.005


def cube_case(ctx, name, bodies, walls=False, np_=8):
    """a cube of side C_L, periodic in x and z, and in y unless walls (then no-slip walls at y = -+C_L/2
    moving at -+1 along x: simple shear 2/C_L), refined C_LEVELS times in one box around all bodies
    (symmetric about the planes through the origin if the bodies are); laminar, Re 0.01 (nu 1/G_RE);
    bodies: name -> ibm body dict. Runs C_T and returns the free-body histories per body."""
    d = ctx.root / name
    f = d / "history.npz"
    if f.exists():
        z = np.load(f)
        return {k: {q: z["%s_%s" % (k, q)] for q in ("t", "U", "Om", "F", "M")} for k in bodies}
    d.mkdir(parents=True, exist_ok=True)
    R = max(b["radius"] for b in bodies.values())
    c = np.array([b["centre"] for b in bodies.values()], float)
    if not (d / "constant" / "polyMesh" / "owner").exists():
        cyc = lambda a, b: {a: ("cyclic", [a[0] + "min"], "neighbourPatch %s;" % b),
                            b: ("cyclic", [a[0] + "max"], "neighbourPatch %s;" % a)}
        patches = {}
        patches.update({"xlo": ("cyclic", ["xmin"], "neighbourPatch xhi;"),
                        "xhi": ("cyclic", ["xmax"], "neighbourPatch xlo;")})
        patches.update({"zlo": ("cyclic", ["zmin"], "neighbourPatch zhi;"),
                        "zhi": ("cyclic", ["zmax"], "neighbourPatch zlo;")})
        if walls:
            patches.update({"ylo": ("wall", ["ymin"]), "yhi": ("wall", ["ymax"])})
        else:
            patches.update({"ylo": ("cyclic", ["ymin"], "neighbourPatch yhi;"),
                            "yhi": ("cyclic", ["ymax"], "neighbourPatch ylo;")})
        boxes = []
        for k in range(1, C_LEVELS + 1):
            m = R + 8*C_H0/2**k
            boxes.append((tuple(np.minimum(c.min(axis=0) - m, -m)), tuple(np.maximum(c.max(axis=0) + m, m))))
        meshes.refined(d, (-C_L/2,)*3, (C_L/2,)*3, C_H0, boxes, patches, two_d=False)
        common.sh("checkMesh", d, "log.checkMesh")
    geometry.ibm_properties(d, bodies, periodicLengths=(C_L, 0 if walls else C_L, C_L))
    for x in ("transportProperties", "turbulenceProperties"):
        (d / "constant" / x).write_text((mig.CYL / "constant" / x).read_text())
    common.set_keys(d / "constant" / "transportProperties", nu=1.0/G_RE)
    fvo = d / "constant" / "fvOptions"
    if fvo.exists():
        fvo.unlink()
    cy = "type cyclic;"
    bU = "xlo { %s } xhi { %s } zlo { %s } zhi { %s } " % ((cy,)*4)
    bp = bU
    if walls:
        bU += "ylo { type fixedValue; value uniform (-1 0 0); } yhi { type fixedValue; value uniform (1 0 0); }"
        bp += "ylo { type zeroGradient; } yhi { type zeroGradient; }"
    else:
        bU += "ylo { %s } yhi { %s }" % (cy, cy)
        bp += "ylo { %s } yhi { %s }" % (cy, cy)
    (d / "0").mkdir(exist_ok=True)
    U0 = "uniform (0 0 0)"
    if walls:
        common.sh("postProcess -func writeCellCentres -time 0", d, "log.writeCellCentres")
        C, _ = read_field(str(d / "0" / "C"))
        U0 = "nonuniform List<vector>\n%d\n(\n%s\n)" % (len(C), "\n".join("(%.12g 0 0)" % (2*y/C_L)
                                                                        for y in C[:, 1]))
    write_dict(d / "0" / "U", "volVectorField", "U",
               "dimensions [0 1 -1 0 0 0 0];\ninternalField %s;\nboundaryField { %s }\n" % (U0, bU))
    write_dict(d / "0" / "p", "volScalarField", "p",
               "dimensions [0 2 -2 0 0 0 0];\ninternalField uniform 0;\nboundaryField { %s }\n" % bp)
    (d / "system" / "fvSchemes").write_text(
        (mig.CYL / "system" / "fvSchemes").read_text().replace("default steadyState;", "default Euler;"))
    write_dict(d / "system" / "fvSolution", "dictionary", "fvSolution",
               pimple_solution((0.45*C_L, 0.45*C_L, 0.45*C_L), NOUTER))
    geometry.control_dict(d, "ibmPimpleFoam", endTime=C_T, deltaT=C_DT, writeInterval=int(round(C_T/C_DT)))
    run_parallel(d, "ibmPimpleFoam", np_)
    log = (d / "log.ibmPimpleFoam").read_text()
    out, z = {}, {}
    for k in bodies:
        h = free_lines(log, k)
        out[k] = {q: h[q] for q in ("t", "U", "Om", "F", "M")}
        z.update({"%s_%s" % (k, q): h[q] for q in ("t", "U", "Om", "F", "M")})
    np.savez(f, **z)
    return out


def free_sphere(centre, U0=(0, 0, 0), Om0=(0, 0, 0), trans=(1, 1, 1), rot=(1, 1, 1), R=0.5, rho=C_RHO):
    return {"type": "sphere", "centre": tuple(centre), "radius": R, "wallType": "noSlip",
            "motion": {"type": "free", "density": rho, "translation": tuple(trans), "rotation": tuple(rot),
                       "velocity": tuple(U0), "angularVelocity": tuple(Om0), "coupling": "implicit",
                       "addedMass": 0.5, "meanPressureGradient": False}}


def gate_sixdof(ctx):
    """a heavy sphere (rho_s/rho 1000) at the centre of the periodic cube, free in all six directions,
    released with a velocity or a spin along x or along (1 1 1)/sqrt 3: by the cubic symmetry of the cell
    the decay along the diagonal is that along x (to O(Re)), the other components stay zero and
    translation and rotation do not couple (the sphere sits at a centre of inversion)"""
    e = np.ones(3)/np.sqrt(3)
    runs = {"U_x": ((1, 0, 0), (0, 0, 0)), "U_diag": (tuple(e), (0, 0, 0)),
            "Om_x": ((0, 0, 0), (1, 0, 0)), "Om_diag": ((0, 0, 0), tuple(e))}
    h = {k: cube_case(ctx, "c6_%s" % k, {"sphere": free_sphere((0, 0, 0), U0, Om0)})["sphere"]
         for k, (U0, Om0) in runs.items()}
    r = {}
    for q, other in (("U", "Om"), ("Om", "U")):
        ax, dg = h["%s_x" % q], h["%s_diag" % q]
        a0 = np.abs(ax[q][:, 0]).max()
        r["%s_x_decay" % q] = float(ax[q][-1, 0]/a0)
        # along the diagonal: the magnitude of the x run, the direction kept
        par = dg[q] @ e
        perp = np.linalg.norm(dg[q] - np.outer(par, e), axis=1)
        r["%s_diag_vs_x" % q] = float(np.abs(par - ax[q][:, 0]).max()/a0)
        r["%s_diag_perp" % q] = float(perp.max()/a0)
        r["%s_x_cross" % q] = float(np.abs(ax[q][:, 1:]).max()/a0)
        r["%s_coupling" % q] = float(max(np.abs(ax[other]).max(), np.abs(dg[other]).max())
                                     /(a0*(1 if q == "U" else 0.5)))
    return r


def gate_pair(ctx):
    """two heavy spheres at y = -+0.75 in the periodic cube released with U = (1 0 0): mirror images
    about y = 0 (U_x, U_z, Om_y equal; U_y, Om_x, Om_z opposite); several free bodies at once"""
    b = {"a": free_sphere((0, -0.75, 0), (1, 0, 0)), "b": free_sphere((0, 0.75, 0), (1, 0, 0))}
    h = cube_case(ctx, "c6_pair", b)
    a, bb = h["a"], h["b"]
    sU, sO = np.array([1, -1, 1]), np.array([-1, 1, -1])
    scale = np.abs(a["U"][:, 0]).max()
    return {"Ux_end": float(a["U"][-1, 0]), "Uy_end": float(a["U"][-1, 1]), "Omz_end": float(a["Om"][-1, 2]),
            "mirror_U": float(np.abs(a["U"] - sU*bb["U"]).max()/scale),
            "mirror_Om": float(np.abs(a["Om"] - sO*bb["Om"]).max()/scale)}


def gate_shear(ctx):
    """a sphere (rho_s/rho 10: spun up within a few steps) at the centre of the wall-bounded shear (walls at
    y = -+2 moving at -+1, rate 0.5),
    free to turn only: about z alone and about all three axes; the torque-free rotation -0.25 (Stokes,
    unbounded) up to the confinement; the two must agree and the free x, y spin stay zero"""
    out = {}
    for tag, rot in (("z", (0, 0, 1)), ("xyz", (1, 1, 1))):
        h = cube_case(ctx, "c6_shear_%s" % tag, {"sphere": free_sphere((0, 0, 0), trans=(0, 0, 0), rot=rot, rho=10.0)},
                      walls=True)["sphere"]
        out[tag] = h
    z, x = out["z"], out["xyz"]
    return {"Omz_end_z": float(z["Om"][-1, 2]), "Omz_end_xyz": float(x["Om"][-1, 2]),
            "Omz_stokes": -0.25, "xyz_vs_z": float(np.abs(x["Om"][:, 2] - z["Om"][:, 2]).max()/0.25),
            "xyz_xy": float(np.abs(x["Om"][:, :2]).max()/0.25)}


def mesh_full(d, hf, e, L=3.0, h0=None):
    """mig.mesh with the refinement box over the whole length (the sphere travels)"""
    h0 = h0 or 2*hf
    W = mig.RP + 2*h0
    R = D/2
    n = int(round(np.log2(h0/hf)))
    boxes = []
    for k in range(1, n + 1):
        m = 8*h0/2**k
        boxes.append(((-L, max(-W, e - R - m), -1.0), (L, W + 1.0, R + m)))
    patches = {"left": ("cyclic", ["xmin"], "neighbourPatch right;"),
               "right": ("cyclic", ["xmax"], "neighbourPatch left;"),
               "sym": ("symmetryPlane", ["zmin"]), "outer": ("wall", ["ymin", "ymax", "zmax"])}
    return meshes.refined(d, (-L/2, -W, 0.0), (L/2, W, W), h0, boxes, patches, two_d=False)


def empty_pipe(ctx, Re, hf, e, h0):
    """G0 and Vf0 of the empty pipe on the case mesh"""
    d = ctx.root / ("pipe_e%.4f_Re%g" % (e, Re))
    f = d / "result.json"
    if f.exists():
        return json.loads(f.read_text())
    d.mkdir(parents=True, exist_ok=True)
    if not (d / "constant" / "polyMesh" / "owner").exists():
        mesh_full(d, hf, e, h0=h0)
    mig.setup(d, Re, hf, D, e, sphere=False, h0=h0)
    # judged by stationarity (G to 1e-10 over 200 iterations, mig.extract): on the mesh refined
    # over the whole length the residual control of mig.setup is not reached
    common.set_keys(d / "system" / "controlDict", endTime=6000, writeInterval=6000)
    Vf0 = mig.geometry_only(d)
    r = dict(mig.solve(d), Vf0=Vf0)
    f.write_text(json.dumps(r))
    return r


def poiseuille(d):
    """u = 2 J (1 - 4 r^2) inside the pipe at the cell centres (0/C)"""
    C, _ = read_field(str(d / "0" / "C"))
    r2 = C[:, 1]**2 + C[:, 2]**2
    u = np.where(r2 < mig.RP**2, 2*mig.J*(1 - 4*r2), 0.0)
    return "nonuniform List<vector>\n%d\n(\n%s\n)" % (len(u), "\n".join("(%.12g 0 0)" % x for x in u))


def free_case(ctx, Re, er, coupling, rho, T=None, dtfac=1.0, np_=None, ptol=None):
    """the free sphere in the still pipe; coupling "prescribed" (diagnostic): the sphere translating
    and turning at the V6 values V, Omega instead (moving geometry without the coupling).
    ptol: outer correctors until the p residual is below ptol (cap NOUTER_CAP), else NOUTER: at Re 2
    the step is viscous-dominated (nu dt/h^2 ~ 8) and the outer loop contracts only ~0.93 per
    corrector (stock pimpleFoam alike), 3 correctors leave f/Re wrong by a factor 2 in the lab frame.
    The alternative, noted for later (user, 2026-11-25): an accelerated outer iteration
    (Aitken/Anderson) or a coupled p-U solve (not in OpenFOAM v2412/v2506; foam-extend pUCoupledFoam)"""
    hf, h0 = v6.mesh_of(D, DH)
    es = float(ref.eps_star(D))
    e = er*es
    name = "rigid_d%g_e%.4f_Re%g_h%d_%s_rho%g%s%s" % (D, e, Re, round(1/hf), coupling, rho,
                                                     "_dtx%g" % dtfac if dtfac != 1.0 else "",
                                                     "_T%g" % T if T else "") + ("_ptol%g" % ptol if ptol else "")
    d = ctx.root / name
    f = d / "result.json"
    if f.exists():
        return json.loads(f.read_text())
    pipe = empty_pipe(ctx, Re, hf, e, None if abs(h0 - 2*hf) < 1e-12 else h0)
    d.mkdir(parents=True, exist_ok=True)
    src = ctx.root / ("pipe_e%.4f_Re%g" % (e, Re))
    if not (d / "constant" / "polyMesh").exists():
        shutil.copytree(src / "constant" / "polyMesh", d / "constant" / "polyMesh")
        shutil.copytree(src / "system", d / "system")
        shutil.copy(src / "log.checkMesh", d / "log.checkMesh")
    L = 3.0
    R = D/2
    V0 = 2*mig.J*(1 - 4*e**2)                 # undisturbed flow at the centre
    Om0 = 8*mig.J*e                           # half its vorticity (-du/dy/2)
    pipe_b = {"type": "cylinder", "centre": (0, 0, 0), "axis": (1, 0, 0), "radius": mig.RP, "invert": True}
    s = {"type": "sphere", "centre": (0, e, 0), "radius": R, "wallType": "noSlip",
         "motion": {"type": "free", "density": rho, "translation": (1, 0, 0), "rotation": (0, 0, 1),
                    "velocity": (V0, 0, 0), "angularVelocity": (0, 0, Om0), "coupling": coupling,
                    "addedMass": 0.5, "volume": np.pi*D**3/12}}
    if coupling == "prescribed":
        p6 = v6_point(Re, er)
        s["motion"] = {"type": "translating", "velocity": (p6["V"], 0, 0), "angularVelocity": (0, 0, p6["Om"])}
    geometry.ibm_properties(d, {"pipe": pipe_b, "sphere": s}, periodicLengths=(L, 0, 0))
    for x in ("transportProperties", "turbulenceProperties"):
        (d / "constant" / x).write_text((mig.CYL / "constant" / x).read_text())
    common.set_keys(d / "constant" / "transportProperties", nu=1.0/Re)
    write_dict(d / "constant" / "fvOptions", "dictionary", "fvOptions", """
momentumSource { type ibmMeanVelocityForce; selectionMode all; fields (U); Ubar (%r 0 0); relaxation 1.0;
                 includeMovingBodies yes; }
""" % mig.J)
    (d / "0").mkdir(exist_ok=True)
    common.sh("postProcess -func writeCellCentres -time 0", d, "log.writeCellCentres")
    bU = "left { type cyclic; } right { type cyclic; } sym { type symmetryPlane; } outer { type noSlip; }"
    bp = "left { type cyclic; } right { type cyclic; } sym { type symmetryPlane; } outer { type zeroGradient; }"
    write_dict(d / "0" / "U", "volVectorField", "U",
               "dimensions [0 1 -1 0 0 0 0];\ninternalField %s;\nboundaryField { %s }\n" % (poiseuille(d), bU))
    write_dict(d / "0" / "p", "volScalarField", "p",
               "dimensions [0 2 -2 0 0 0 0];\ninternalField uniform 0;\nboundaryField { %s }\n" % bp)
    (d / "system" / "fvSchemes").write_text(
        (mig.CYL / "system" / "fvSchemes").read_text().replace("default steadyState;", "default Euler;"))
    write_dict(d / "system" / "fvSolution", "dictionary", "fvSolution", pimple_solution((-0.4*L, -0.21, 0.013),
               NOUTER_CAP if ptol else NOUTER, tight=False, ptol=ptol))
    dt = 0.4*hf/(2*mig.J)*dtfac                # Courant 0.4 on the centre-line velocity
    dt = float("%.3g" % dt)
    T = T or (3.0 if Re <= 8 else 15.0)
    n = int(round(T/dt))                       # writeInterval counts time steps
    geometry.control_dict(d, "ibmPimpleFoam", endTime=T, deltaT=dt, writeInterval=max(1, n//4))
    cells = common.last_match((d / "log.checkMesh").read_text(), r"cells:\s+(\d+)", cast=int)
    run_parallel(d, "ibmPimpleFoam", np_ or NP_CASE)
    log = (d / "log.ibmPimpleFoam").read_text()
    if coupling == "prescribed":
        m = re.findall(r"IBM force on sphere: pressure \(([^)]*)\) viscous \(([^)]*)\)", log)
        a = np.array([[float(x) for x in (p + " " + v).split()] for p, v in m])
        t = np.array([float(x) for x in re.findall(r"^Time = (\S+)", log, re.M)][:len(a)])
        np.save(d / "forces.npy", np.column_stack([t, a]))
        VB = np.pi*D**3/6
        Fy, dFy, sFy = window(t, a[:, 1] + a[:, 4])
        r = {"name": name, "Re": Re, "epsr": er, "e": e, "T": T, "dt": dt, "steps": len(a),
             "nouter": log.count("PIMPLE: iteration")/max(len(a), 1),
             "fRe": float(2*Fy/VB), "dfRe": float(2*dFy/VB), "sfRe": float(2*sFy/VB),
             "fRe_v6": v6_point(Re, er)["fRe"]}
        f.write_text(json.dumps(r))
        return r
    h = free_lines(log)
    g = gradients(log)
    n = min(len(g), len(h["t"]))
    t = h["t"][-n:]
    VB = np.pi*D**3/6
    V, dV, sV = window(t, h["U"][-n:, 0])
    Om, dOm, sOm = window(t, h["Om"][-n:, 2])
    Fy, _, sFy = window(t, h["F"][-n:, 1])
    G, dG, sG = window(t, g[-n:])
    aG = VB/(L*np.pi/4)
    r = {"name": name, "Re": Re, "epsr": er, "e": e, "coupling": coupling, "rho": rho, "T": T, "dt": dt,
         "steps": len(h["t"]), "cells": cells, "nouter": log.count("PIMPLE: iteration")/max(len(h["t"]), 1),
         "V": float(V), "dV": float(dV), "sV": float(sV), "Om": float(Om), "dOm": float(dOm), "sOm": float(sOm),
         "fRe": float(2*Fy/VB), "sfRe": float(2*sFy/VB), "G": float(G), "G0": pipe["G"],
         "beta": float((G/pipe["G"] - 1)/aG), "settled": bool(abs(dV) <= SETTLE*abs(V) and
                                                               abs(dOm) <= SETTLE*max(abs(Om), 1e-3)),
         "Fx_res": float(np.abs(h["F"][-n:, 0]).mean()*2/VB), "Mz_res": float(np.abs(h["M"][-n:, 2]).mean()*2/VB)}
    f.write_text(json.dumps(r))
    return r


def frame_case(ctx, Re, er, T=None):
    """diagnostic: the V6 configuration (sphere still, turning at the V6 Omega; pipe sliding at -V6;
    Ubar (J - V) Vf0/Vf) run transiently with ibmPimpleFoam on the V6G mesh and time step: the
    lift of V6 with the solver and mesh of V6G but a static geometry"""
    hf, h0 = v6.mesh_of(D, DH)
    e = er*float(ref.eps_star(D))
    p6 = v6_point(Re, er)
    V, Om = p6["V"], p6["Om"]
    name = "frame_d%g_e%.4f_Re%g_h%d" % (D, e, Re, round(1/hf))
    d = ctx.root / name
    f = d / "result.json"
    if f.exists():
        return json.loads(f.read_text())
    pipe = empty_pipe(ctx, Re, hf, e, None if abs(h0 - 2*hf) < 1e-12 else h0)
    d.mkdir(parents=True, exist_ok=True)
    src = ctx.root / ("pipe_e%.4f_Re%g" % (e, Re))
    if not (d / "constant" / "polyMesh").exists():
        shutil.copytree(src / "constant" / "polyMesh", d / "constant" / "polyMesh")
        shutil.copytree(src / "system", d / "system")
        shutil.copy(src / "log.checkMesh", d / "log.checkMesh")
    L = 3.0
    mig.setup(d, Re, hf, D, e, "rigid", V, Om, L)          # V6 bodies and fvOptions (Ubar J - V)
    Vf = mig.geometry_only(d)
    ub = float((mig.J - V)*pipe["Vf0"]/Vf)
    (d / "constant" / "fvOptions").write_text(re.sub(r"Ubar \([^)]*\)", "Ubar (%r 0 0)" % ub,
                                                    (d / "constant" / "fvOptions").read_text()))
    C, _ = read_field(str(d / "0" / "C"))
    r2 = C[:, 1]**2 + C[:, 2]**2
    u = np.where(r2 < mig.RP**2, 2*mig.J*(1 - 4*r2), 0.0) - V
    bU = "left { type cyclic; } right { type cyclic; } sym { type symmetryPlane; } outer { type noSlip; }"
    write_dict(d / "0" / "U", "volVectorField", "U", "dimensions [0 1 -1 0 0 0 0];\ninternalField "
               "nonuniform List<vector>\n%d\n(\n%s\n);\nboundaryField { %s }\n"
               % (len(u), "\n".join("(%.12g 0 0)" % x for x in u), bU))
    (d / "system" / "fvSchemes").write_text(
        (mig.CYL / "system" / "fvSchemes").read_text().replace("default steadyState;", "default Euler;"))
    write_dict(d / "system" / "fvSolution", "dictionary", "fvSolution", pimple_solution((-0.4*L, -0.21, 0.013), NOUTER, tight=False))
    dt = float("%.3g" % (0.4*hf/(2*mig.J)))
    T = T or 1.5
    n = int(round(T/dt))
    geometry.control_dict(d, "ibmPimpleFoam", endTime=T, deltaT=dt, writeInterval=max(1, n//4))
    run_parallel(d, "ibmPimpleFoam", 16)
    log = (d / "log.ibmPimpleFoam").read_text()
    m = re.findall(r"IBM force on sphere: pressure \(([^)]*)\) viscous \(([^)]*)\)", log)
    a = np.array([[float(x) for x in (p + " " + v).split()] for p, v in m])
    a = a[len(a)//2:]
    VB = np.pi*D**3/6
    r = {"name": name, "V": V, "Om": Om, "steps": len(m), "Fp": a[:, :3].mean(axis=0).tolist(),
         "Fv": a[:, 3:].mean(axis=0).tolist(), "fRe": float(2*(a[:, 1] + a[:, 4]).mean()/VB),
         "fRe_v6": p6["fRe"]}
    f.write_text(json.dumps(r))
    return r


def v6_point(Re, er):
    """the V6 baseline point at this (Re, eps/eps*) on the d/h 32 mesh"""
    b = common.load_baseline("V6")
    for p in b["data"]["points"]:
        if p["interface"] == "rigid" and p["d"] == D and p["Re"] == Re and p["dh"] == DH and abs(p["epsr"] - er) < 1e-6:
            return p
    return None


# --------------------------------------------------------------------------------------- run

def run(ctx, res):
    st = gate_static(ctx)
    gates = common.pmap(lambda j: gate_free(ctx, *j, st, GATE_DT), GATE_RUNS, max(1, ctx.jobs//NP_GATE))
    scan = common.pmap(lambda dt: gate_free(ctx, "implicit", 1.0, st, dt), GATE_DT_SCAN,
                       max(1, ctx.jobs//NP_GATE))
    jobs = [(Re, er, "implicit", 1.0) for Re, er in POINTS]
    jobs += [(POINTS[0][0], POINTS[0][1], c, rho) for c, rho in VARIANTS if (c, rho) != ("implicit", 1.0)]
    jobs += [(POINTS[0][0], POINTS[0][1], "implicit", 1.0, None, 0.5)]   # time-step error (G1: ~dt^2)
    cases = common.pmap(lambda j: free_case(ctx, *j), jobs, max(1, ctx.jobs//NP_CASE))
    res.data["gate_static"] = st
    res.data["gate"] = gates
    res.data["gate_dt_scan"] = scan
    res.data["cases"] = cases
    res.data["metrics"] = {}
    for g in gates:
        res.rel("G1 free sphere, Stokes box, %s rho_s/rho %g: terminal velocity (static drag on the same mesh)"
                % (g["coupling"], g["rho"]), g["Ut"], 1.0, GATE_TOL,
                "drift %.1e, std %.1e, force balance %.1e, %d steps" % (g["dUt"], g["sUt"], g["Fbal"], g["steps"]))
    for g in scan:
        res.report("G1 implicit rho_s/rho 1, dt %g: terminal velocity" % g["dt"], g["Ut"], 1.0,
                   "deviation %+.3g %%" % (100*(g["Ut"] - 1)))
    for c in cases:
        p = v6_point(c["Re"], c["epsr"])
        tag = "Re %g eps/eps* %g %s rho_s/rho %g dt %g" % (c["Re"], c["epsr"], c["coupling"], c["rho"], c["dt"])
        res.le("%s: settled (V, Omega over the halves of the last third)" % tag,
               max(abs(c["dV"])/abs(c["V"]), abs(c["dOm"])/max(abs(c["Om"]), 1e-3)), SETTLE)
        if p is None:
            res.report("%s: no V6 point on this mesh" % tag, c["V"], None)
            continue
        res.rel("%s: V vs V6" % tag, c["V"], p["V"], TOL["V"], "std %.1e" % c["sV"])
        res.rel("%s: Omega vs V6" % tag, c["Om"], p["Om"], TOL["Om"], "std %.1e" % c["sOm"])
        res.absolute("%s: f/Re vs V6" % tag, c["fRe"], p["fRe"], TOL["fRe"], "std %.1e" % c["sfRe"])
        res.rel("%s: beta vs V6" % tag, c["beta"], p["beta"], TOL["beta"], "G %.6g G0 %.6g" % (c["G"], c["G0"]))
