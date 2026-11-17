"""V6 MIG: inertial migration of a sphere in a pipe (Rivero-Rodriguez & Scheid 2018; specification
doc/validationSuite.md §8; driver mig.py, reference validation/reference_riveroScheid2018.py).

Gates first (V6_GATES_ONLY=1 in the environment: only these), then the cases.

MIG-G1 pipe: the empty pipe d_h = 1 (inverted IBM cylinder) sliding at -V = -0.7 along its axis,
Re = 2, flux J = 1 (fluid-volume average J - V): exact u = -V + 2 J (1 - 4 r^2), G = 32/Re. On the
plain pipe meshes of the MIG meshes' outer level (h = 1/40, 1/60, 1/80: the pipe wall of d = 0.4 at
d/h = 32, 48, 64): G within 0.2 % on the base mesh used (d/h 48), order of G and of the profile
error (alpha-weighted L2 over the fluid) >= 1.5.

MIG-G2 sphere, Stokes (Re = 0.01, U = d = 1): a sphere in a periodic box of side L (x, y
periodic, z in [0, L/2] with symmetry planes at both ends; the box is a cubic array), translating
(mean flow U through the fixed sphere: drag 6 pi mu R U), stress-free (4 pi mu R U,
Hadamard-Rybczynski) and rotating at Omega = 1 about z (torque 8 pi mu R^3 Omega). The array
correction (Hasimoto 1959: 1/(1 - 1.7601 s' + ...), s' = c^(1/3), c = (4 pi/3)(R/L)^3) is removed
by fitting F/F_Stokes = 1/(a0 - a1 s + a3 s^3) (s = R/L) through L/d = 4, 6, 8 (rigid, stress-free)
and T/T_Stokes = t0 + t3 s^3 (rotating). Drag of the array: the force on the sphere plus G V_B (the
driving gradient on the sphere volume), U the superficial velocity Ubar V_f/L^3. d/h = 16, 32, 48;
rigid and torque within 1 % at d/h 32 and decreasing to d/h 48; stress-free reported against 1 %.
"""
import json
import os
from pathlib import Path

import numpy as np

import common
import geometry
import meshes
import mig
import reference_riveroScheid2018 as ref
from common import write_dict
from foamio import latest_time, read_field

ID, KIND, TIER = "V6", "case", 2
TITLE = "MIG: inertial migration of a sphere in a pipe (gates G1, G2; cases R1-R4, S1, S2)"
STRESSES, RELIES = ["K1", "K3", "K4", "K7", "K8", "K10"], ["V10", "V11", "V1"]
REGRESSION = (1e-4,)

G1_RE, G1_V = 2.0, 0.7
G1_H = [1/40, 1/60, 1/80]
G2_RE = 0.01
G2_DHS = [16, 32, 48]
G2_LS = [4.0, 6.0, 8.0]
G2_KINDS = ["rigid", "stressFree", "rotating"]


# ---------------------------------------------------------------------------------------- G1
def g1(ctx, h0):
    d = ctx.root / ("g1_pipe_h%d" % round(1/h0))
    if not (d / "result.json").exists():
        d.mkdir(parents=True, exist_ok=True)
        mig.setup(d, G1_RE, h0/2, None, V=G1_V)      # hf = h0/2: the outer level of the MIG mesh
        # p is uniform up to round-off in the empty pipe: its normalised residual stays ~1e-4, the
        # solution is judged by stationarity (G steady to 1e-10 from ~500 iterations at h = 1/40)
        common.set_keys(d / "system" / "controlDict", endTime=3000, writeInterval=3000)
        mig.geometry_only(d)
        r = mig.solve(d)
        (d / "result.json").write_text(json.dumps(r))
    r = json.loads((d / "result.json").read_text())
    r.update(profile_error(d))
    r["h"] = h0
    return r


def profile_error(d):
    """alpha-weighted L2 and max error of u_x against -V + 2 J (1 - 4 r^2) at the value points,
    and |u_y, u_z|"""
    t = latest_time(str(d))
    U, _ = read_field(os.path.join(d, t, "U"))
    C, _ = read_field(os.path.join(d, "0", "C"))
    a, _ = read_field(os.path.join(d, "0", "ibmAlpha"))
    v, _ = read_field(os.path.join(d, "0", "V"))
    xv = common.value_points(d, t, C)          # where each cell's value lives (the fluid centroid of cut cells)
    r2 = xv[:, 1]**2 + xv[:, 2]**2
    ue = -G1_V + 2*mig.J*(1 - 4*r2)
    w = a*v
    m = a > 0.5
    return {"uL2": float(np.sqrt(np.sum(w*(U[:, 0] - ue)**2)/np.sum(w))),
            "uMax": float(np.abs(U[m, 0] - ue[m]).max()),
            "uCross": float(np.sqrt(np.sum(w*(U[:, 1]**2 + U[:, 2]**2))/np.sum(w)))}


# ---------------------------------------------------------------------------------------- G2
def g2_mesh(d, L, dh):
    """box [-L/2, L/2]^2 x [0, L/2] (d = 1), background h0 = 1/3 ... (L divisible), nested boxes
    halving to h = 1/dh; box k (cells h_k) reaches 4 cells of the level outside it (8 h_k) beyond the
    sphere"""
    hf = 1.0/dh
    n = int(round(np.log2(dh/4.0)))           # levels: h0 = 2^n hf ~ 1/4 ... 1/2
    h0 = hf*2**n
    N = int(round(L/h0))
    assert abs(N*h0 - L) < 1e-9, "L not a multiple of h0"
    boxes = []
    for k in range(1, n + 1):                  # box k refines to h0/2^k
        hk = h0/2**k
        w = 0.5 + 8*hk
        boxes.append(((-w, -w, -1.0), (w, w, w)))
    patches = {"left": ("cyclic", ["xmin"], "neighbourPatch right;"),
               "right": ("cyclic", ["xmax"], "neighbourPatch left;"),
               "front": ("cyclic", ["ymin"], "neighbourPatch back;"),
               "back": ("cyclic", ["ymax"], "neighbourPatch front;"),
               "sym": ("symmetryPlane", ["zmin"]), "top": ("symmetryPlane", ["zmax"])}
    return meshes.refined(d, (-L/2, -L/2, 0.0), (L/2, L/2, L/2), h0, boxes, patches, two_d=False)


def g2_setup(d, L, dh, kind, **ibm):
    d = Path(d)
    if not (d / "constant" / "polyMesh" / "owner").exists():
        g2_mesh(d, L, dh)
    s = {"type": "sphere", "centre": (0, 0, 0), "radius": 0.5,
         "wallType": "freeSlip" if kind == "stressFree" else "noSlip"}
    if kind == "rotating":
        s["motion"] = {"type": "none", "angularVelocity": (0, 0, 1.0), "invariant": True}
    geometry.ibm_properties(d, {"sphere": s}, **ibm)
    for x in ("transportProperties", "turbulenceProperties"):
        (d / "constant" / x).write_text((mig.CYL / "constant" / x).read_text())
    common.set_keys(d / "constant" / "transportProperties", nu=1.0/G2_RE)
    ub = 0.0 if kind == "rotating" else 1.0
    write_dict(d / "constant" / "fvOptions", "dictionary", "fvOptions", """
momentumSource { type ibmMeanVelocityForce; selectionMode all; fields (U); Ubar (%r 0 0); flowDir (1 0 0);
                 relaxation 1.0; }
""" % ub)   # rotating: Ubar 0 needs the direction; the zero mean along x holds by symmetry anyway
    bc = ("left { type cyclic; } right { type cyclic; } front { type cyclic; } back { type cyclic; } "
          "sym { type symmetryPlane; } top { type symmetryPlane; }")
    (d / "0").mkdir(exist_ok=True)
    write_dict(d / "0" / "U", "volVectorField", "U",
               "dimensions [0 1 -1 0 0 0 0];\ninternalField uniform (%r 0 0);\nboundaryField { %s }\n" % (ub, bc))
    write_dict(d / "0" / "p", "volScalarField", "p",
               "dimensions [0 2 -2 0 0 0 0];\ninternalField uniform 0;\nboundaryField { %s }\n" % bc)
    (d / "system" / "fvSchemes").write_text((mig.CYL / "system" / "fvSchemes").read_text())
    write_dict(d / "system" / "fvSolution", "dictionary", "fvSolution", """
solvers { p { solver PCG; preconditioner { preconditioner linkedGroup; base { preconditioner GAMG;
                smoother DICGaussSeidel; nCellsInCoarsestLevel 10; tolerance 0; relTol 0.01; } }
              tolerance 1e-13; relTol 0.01; maxIter 1000; }
          U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-13; relTol 0.1; maxIter 50; } }
SIMPLE { nNonOrthogonalCorrectors 0; consistent yes; pRefPoint (%r %r %r); pRefValue 0;
         residualControl { p 1e-9; U 1e-9; } }
relaxationFactors { equations { U 0.9; } fields { p 1.0; } }
""" % (0.45*L, 0.45*L, 0.2*L))
    geometry.control_dict(d, "ibmSimpleFoam", endTime=20000, deltaT=1, writeInterval=20000)


def g2(ctx, kind, L, dh):
    d = ctx.root / ("g2_%s_L%g_dh%d" % (kind, L, dh))
    if not (d / "result.json").exists():
        d.mkdir(parents=True, exist_ok=True)
        g2_setup(d, L, dh, kind)
        Vf = mig.geometry_only(d)
        r = mig.solve(d)
        r["Vf"] = 2*Vf
        (d / "result.json").write_text(json.dumps(r))
    r = json.loads((d / "result.json").read_text())
    mu, R = 1.0/G2_RE, 0.5
    VB = 4*np.pi/3*R**3
    if kind == "rotating":
        val, ref = 2*r["T"][2], -8*np.pi*mu*R**3        # torque on the sphere opposes Omega
    else:
        U = 1.0*r["Vf"]/L**3                              # superficial velocity
        val = 2*r["F"][0] + r["G"]*VB                     # force of the array on one sphere
        ref = (6 if kind == "rigid" else 4)*np.pi*mu*R*U
    r.update(kind=kind, L=L, dh=dh, value=val, stokes=ref, ratio=val/ref, s=R/L)
    return r


def g2_fit(rs, kind):
    """value at s = R/L -> 0 from the three boxes: 1/ratio = a0 - a1 s + a3 s^3 (translating) or
    ratio = t0 + t3 s^3 (rotating); returns the extrapolated ratio"""
    s = np.array([r["s"] for r in rs])
    q = np.array([r["ratio"] for r in rs])
    if kind == "rotating":
        A = np.column_stack([np.ones_like(s), s**3])
        return float(np.linalg.lstsq(A, q, rcond=None)[0][0])
    A = np.column_stack([np.ones_like(s), -s, s**3])
    return float(1.0/np.linalg.solve(A, 1.0/q)[0])


# ---------------------------------------------------------------------------------------- cases
DHS = [32, 48, 64]              # d/h of the sphere; judged on the finest
RES2 = [8.0, 32.0, 128.0]
EPS = [0.0, 0.2, 0.4, 0.6, 0.75, 0.9]
TOL = 3e-4                      # force balance |F_x|/V_B, |T_z|/V_B: 1e-3 of |f/Re| ~ 0.3
NPDH = {32: 8, 48: 16, 64: 32}  # ranks (mig.ranks of the d = 0.4 meshes)


def fig9():
    t = np.loadtxt(Path(ref.__file__).parent / "data" / "riveroScheid2018" / "fig9_rigid_d0.4.dat")
    return {k: t[:, i + 1] for i, k in enumerate(["lin", 8.0, 16.0, 32.0, 64.0, 128.0])}, t[:, 0]


def fig10(interface, Re, d):
    f = "fig10b_rigid_neutral.dat" if interface == "rigid" else "fig10d_stressFree_neutral.dat"
    t = np.loadtxt(Path(ref.__file__).parent / "data" / "riveroScheid2018" / f)
    m = (t[:, 0] == Re) & (t[:, 1] == 0)
    o = np.argsort(t[m, 2])
    return float(np.interp(d, t[m, 2][o], t[m, 3][o]))


def mesh_of(d, dh):
    hf = d/dh
    return hf, (2*hf if d >= 0.3 else 4*hf)     # d = 0.2: two levels, the pipe wall at h0 = d/(dh/4)


def x0_of(interface, d, e):
    V = float(ref.V0(interface, d, e))
    return [V, float(ref.Omega0(d, e))] if interface == "rigid" else [V]


def chain(ctx, interface, d, Re, dh, epss, tag="", **ibm):
    """balanced points at eps/eps_* in epss (one mesh family member), the Jacobian handed on"""
    hf, h0 = mesh_of(d, dh)
    root = ctx.root / "cases"
    out, Jac = [], None
    for s in epss:
        e = s*float(ref.eps_star(d))
        P = mig.Point(root, interface, d, e, Re, hf, h0=h0, tag=tag, **ibm)
        x0 = x0_of(interface, d, e) if not out else _guess(out, interface, d, e)
        r, Jac = P.balance(x0, TOL, dx=[0.02, 0.05][:len(x0)], Jac=Jac)
        r.update(dh=dh, tag=tag, epsr=s)
        out.append(r)
    return out


def _guess(prev, interface, d, e):
    """the fit shifted by the discrete offset of the previous point"""
    p = prev[-1]
    x = x0_of(interface, d, e)
    x[0] += p["V"] - float(ref.V0(interface, d, p["e"]))
    if len(x) > 1:
        x[1] += p["Om"] - float(ref.Omega0(d, p["e"]))
    return x


def neutral(ctx, interface, d, Re, dh, e_ref):
    """secant on eps for f = 0 from e_ref -/+ 0.05 eps_*; returns the points tried and eps_eq"""
    es = float(ref.eps_star(d))
    hf, h0 = mesh_of(d, dh)
    pts, Jac = [], None
    xs = [e_ref - 0.05*es, e_ref + 0.05*es]
    for k in range(6):
        e = xs[k] if k < 2 else xs[-1]
        P = mig.Point(ctx.root / "cases", interface, d, e, Re, hf, h0=h0)
        r, Jac = P.balance(x0_of(interface, d, e) if not pts else _guess(pts, interface, d, e), TOL,
                           dx=[0.02, 0.05][:1 + (interface == "rigid")], Jac=Jac)
        r.update(dh=dh, tag="", epsr=e/es)
        pts.append(r)
        if k >= 1:
            a, b = pts[-2], pts[-1]
            enew = b["e"] - b["fRe"]*(b["e"] - a["e"])/(b["fRe"] - a["fRe"])
            if abs(enew - b["e"]) < 0.002:
                return pts, float(enew)
            xs.append(enew)
    return pts, float(xs[-1])


def case_jobs(quick):
    """(label, function, args, ranks); the finest mesh for every point, all three on the
    convergence points, options on the middle mesh"""
    f = DHS[-1] if not quick else DHS[0]
    J = []
    if quick:
        return [("R1", chain, ("rigid", 0.4, 2.0, f, [0.4]), NPDH[f])]
    J.append(("R1", chain, ("rigid", 0.4, 2.0, f, EPS), NPDH[f]))
    for dh in DHS[:-1]:
        J.append(("R1conv", chain, ("rigid", 0.4, 2.0, dh, [0.4, 0.9]), NPDH[dh]))
    J.append(("R1Re1", chain, ("rigid", 0.4, 1.0, DHS[1], [0.4]), NPDH[DHS[1]]))
    for Re in RES2:
        J.append(("R2", chain, ("rigid", 0.4, Re, f, EPS[1:]), NPDH[f]))
    for dh in DHS[:-1]:
        J.append(("R2conv", chain, ("rigid", 0.4, 128.0, dh, [0.4]), NPDH[dh]))
    for dd in (0.2, 0.6):
        for Re in (2.0, 128.0):
            e = ref.neutral("rigid", dd)[0][0] if Re == 2.0 else fig10("rigid", 128, dd)
            J.append(("R3", neutral, ("rigid", dd, Re, f, e), NPDH[f]))
    for dd in (0.6, 0.8):
        J.append(("R4", chain, ("rigid", dd, 2.0, f, [0.1]), NPDH[f]))
    J.append(("S1", chain, ("stressFree", 0.4, 2.0, f, EPS), NPDH[f]))
    for Re in (2.0, 32.0):
        e = ref.neutral("stressFree", 0.2)[0][0] if Re == 2.0 else fig10("stressFree", 32, 0.2)
        J.append(("S2", neutral, ("stressFree", 0.2, Re, f, e), NPDH[f]))
    m = DHS[1]
    for tag, opt in (("_noconv", {"ibmConvection": False}), ("_merge0", {"mergeAlpha": 0.0})):
        J.append(("R1opt", chain, ("rigid", 0.4, 2.0, m, [0.4, 0.9], tag), NPDH[m], opt))
        J.append(("R2opt", chain, ("rigid", 0.4, 128.0, m, [0.4], tag), NPDH[m], opt))
    return J


def run_cases(ctx, res):
    import threading
    jobs = case_jobs(ctx.quick)
    free, cv = [max(ctx.jobs, max(j[3] for j in jobs))], threading.Condition()

    def one(j):
        with cv:
            cv.wait_for(lambda: free[0] >= j[3])
            free[0] -= j[3]
        try:
            return j[0], j[1](ctx, *j[2], **(j[4] if len(j) > 4 else {}))
        finally:
            with cv:
                free[0] += j[3]
                cv.notify_all()
    out = common.pmap(one, jobs, len(jobs))
    pts = [r for lab, o in out if lab not in ("R3", "S2") for r in o]
    neu = [(lab, j[2], o) for (lab, o), j in zip(out, jobs) if lab in ("R3", "S2")]
    pts += [r for _, _, o in neu for r in o[0]]
    res.data["points"] = pts
    res.data["neutral"] = [{"case": lab, "interface": a[0], "d": a[1], "Re": a[2], "dh": a[3], "e_ref": a[4],
                            "e_eq": o[1]} for lab, a, o in neu]
    evaluate(res, out, jobs)


def evaluate(res, out, jobs):
    M = res.data.setdefault("metrics", {})
    pts = res.data["points"]
    for r in pts:
        k = "%s_d%g_Re%g_dh%d_s%.2f%s" % (r["interface"], r["d"], r["Re"], r["dh"], r["epsr"], r["tag"])
        M[k + "_fRe"], M[k + "_V"] = r["fRe"], r["V"]
    res.le("cases: points not balanced or not converged", sum(not (r["balanced"] and r["converged"]) for r in pts),
           0, " ".join(r["name"] for r in pts if not (r["balanced"] and r["converged"])))
    fmax = lambda itf, d: float(np.max(ref.f1(itf, d, np.linspace(0, 0.95, 96)*ref.eps_star(d))))
    F9, s9 = fig9()
    F4 = ref.fig4()
    sel = lambda itf, d, Re, dh, tag="": sorted([r for r in pts if r["interface"] == itf and r["d"] == d and
                                                 r["Re"] == Re and r["dh"] == dh and r["tag"] == tag],
                                                key=lambda r: r["epsr"])
    fin = max(r["dh"] for r in pts)
    # R1: f/Re against figure 9 (linear), V, Omega, beta against figure 4 (fit in the notes)
    for r in sel("rigid", 0.4, 2.0, fin):
        s = r["epsr"]
        lab = "R1 d 0.4 eps/eps* %.2f Re 2, d/h %d" % (s, fin)
        fl, ff = float(np.interp(s, s9, F9["lin"])), float(ref.f1("rigid", 0.4, r["e"]))
        res.absolute(lab + ": f/Re (figure 9, linear)", r["fRe"], fl, 0.03*fmax("rigid", 0.4), "fit %.4f" % ff)
        V4, b4, O4 = (float(np.interp(s, F4[:, 0], F4[:, i])) for i in (1, 2, 3))
        res.rel(lab + ": V (figure 4b)", r["V"], V4, 0.005, "fit %.5f" % float(ref.V0("rigid", 0.4, r["e"])))
        if s > 0:
            res.rel(lab + ": Omega (figure 4d)", r["Om"], O4, 0.01, "fit %.5f" % float(ref.Omega0(0.4, r["e"])))
        aG = r["beta_aG"]/r["beta"] if r["beta"] else 0
        chk = res.rel if abs(b4*aG) > 1e-4 else res.report
        chk(lab + ": beta (figure 4c)", r["beta"], b4, *((0.05,) if chk == res.rel else ()),
            note="fit %.4f; beta alpha_G %.2e" % (float(ref.beta0("rigid", 0.4, r["e"])), r["beta_aG"]))
    # convergence (Richardson over the three meshes) and the options
    for itf, Re, ss in (("rigid", 2.0, (0.4, 0.9)), ("rigid", 128.0, (0.4,))):
        for s in ss:
            rs = sorted([r for r in pts if r["interface"] == itf and r["d"] == 0.4 and r["Re"] == Re and
                         not r["tag"] and abs(r["epsr"] - s) < 1e-9], key=lambda r: r["dh"])
            if len(rs) < 3:
                continue
            for q in ("fRe", "V", "Om", "beta"):
                Ce, p = common.richardson([1/r["dh"] for r in rs], [r[q] for r in rs])
                res.report("Re %g eps/eps* %.2f: %s on d/h %s, extrapolated (order)" % (
                    Re, s, q, "/".join(str(r["dh"]) for r in rs)), Ce,
                    note="%s; order %.2f" % (" ".join("%.5g" % r[q] for r in rs), p))
            b = [r for r in rs if r["dh"] == DHS[1]][0]
            for tag in ("_noconv", "_merge0"):
                o = [r for r in pts if r["tag"] == tag and r["Re"] == Re and abs(r["epsr"] - s) < 1e-9]
                if o:
                    res.report("Re %g eps/eps* %.2f, d/h %d: f/Re with %s against the default" % (Re, s, DHS[1], tag),
                               o[0]["fRe"] - b["fRe"], note="%.5g against %.5g; V %.6g against %.6g"
                               % (o[0]["fRe"], b["fRe"], o[0]["V"], b["V"]))
    r1 = [r for r in pts if r["Re"] == 1.0]
    if r1:
        b = [r for r in pts if r["Re"] == 2.0 and r["dh"] == r1[0]["dh"] and not r["tag"] and
             r["interface"] == "rigid" and r["d"] == 0.4 and abs(r["epsr"] - 0.4) < 1e-9]
        if b:
            res.report("R1 f/Re at Re 1 against Re 2 (linear regime), d/h %d" % r1[0]["dh"], r1[0]["fRe"] - b[0]["fRe"],
                       note="%.5g against %.5g" % (r1[0]["fRe"], b[0]["fRe"]))
    # R2: f/Re against figure 9
    for Re in RES2:
        for r in sel("rigid", 0.4, Re, fin):
            res.absolute("R2 d 0.4 Re %g eps/eps* %.2f, d/h %d: f/Re (figure 9)" % (Re, r["epsr"], fin), r["fRe"],
                     float(np.interp(r["epsr"], s9, F9[Re])), 0.03*fmax("rigid", 0.4), "V %.5f" % r["V"])
    # R3, S2: neutral positions
    for n in res.data["neutral"]:
        e_ref = n["e_ref"]
        src = "fit" if n["Re"] == 2.0 else "figure 10%s" % ("b" if n["interface"] == "rigid" else "d")
        res.absolute("%s %s d %g Re %g, d/h %d: neutral position eps_eq (%s)" % (
            "R3" if n["interface"] == "rigid" else "S2", n["interface"], n["d"], n["Re"], n["dh"], src),
            n["e_eq"], e_ref, 0.01)
    # R4: sign and value of f at eps/eps* 0.1
    for d in (0.6, 0.8):
        for r in sel("rigid", d, 2.0, fin):
            fr = float(ref.f1("rigid", d, r["e"]))
            res.ge("R4 d %g eps/eps* 0.1, d/h %d: sign of f/Re as the fit (%+.4f)" % (d, fin, fr),
                   float(np.sign(r["fRe"])*np.sign(fr)), 1.0, "f/Re %+.5f" % r["fRe"])
            res.report("R4 d %g eps/eps* 0.1, d/h %d: f/Re (fit)" % (d, fin), r["fRe"], fr)
    # S1: stress-free d 0.4 against the fit; beta reported
    for r in sel("stressFree", 0.4, 2.0, fin):
        lab = "S1 d 0.4 eps/eps* %.2f Re 2, d/h %d" % (r["epsr"], fin)
        res.absolute(lab + ": f/Re (fit)", r["fRe"], float(ref.f1("stressFree", 0.4, r["e"])),
                 0.03*fmax("stressFree", 0.4))
        res.rel(lab + ": V (fit)", r["V"], float(ref.V0("stressFree", 0.4, r["e"])), 0.005)
        res.report(lab + ": beta (fit, table 4)", r["beta"], float(ref.beta0("stressFree", 0.4, r["e"])),
                   "beta alpha_G %.2e" % r["beta_aG"])


# ---------------------------------------------------------------------------------------- run
def run(ctx, res):
    ctx.gates_only = getattr(ctx, "gates_only", None) or bool(os.environ.get("V6_GATES_ONLY"))
    hs = G1_H[:2] if ctx.quick else G1_H
    dhs = G2_DHS[:1] if ctx.quick else G2_DHS
    jobs = [("g1", h, None) for h in hs] + [("g2", (k, L), dh) for k in G2_KINDS for L in G2_LS for dh in dhs]

    def one(j):
        return g1(ctx, j[1]) if j[0] == "g1" else g2(ctx, j[1][0], j[1][1], j[2])
    out = common.pmap(one, jobs, max(1, ctx.jobs // 8))
    G1 = sorted([r for j, r in zip(jobs, out) if j[0] == "g1"], key=lambda r: -r["h"])
    G2 = [r for j, r in zip(jobs, out) if j[0] == "g2"]
    res.data["g1"], res.data["g2"] = G1, G2
    M = res.data["metrics"] = {}

    # G1
    Gex = 32.0/G1_RE
    for r in G1:
        M["g1_h%d_G" % round(1/r["h"])] = r["G"]
    note = "G %s against %g; u L2 %s; cross flow %s" % (" ".join("%.6g" % r["G"] for r in G1), Gex,
                                                        " ".join("%.3g" % r["uL2"] for r in G1),
                                                        " ".join("%.2g" % r["uCross"] for r in G1))
    base = [r for r in G1 if abs(r["h"] - 1/60) < 1e-9] or G1[-1:]
    res.rel("G1 pipe pressure gradient, h = 1/%d (pipe wall of d/h 48)" % round(1/base[0]["h"]),
            base[0]["G"], Gex, 0.002, note)
    res.le("G1 runs not converged", sum(not r["converged"] for r in G1), 0)
    if len(G1) >= 2:
        a, b = G1[-2], G1[-1]
        for q, lab in (("G", "G"), ("uL2", "profile error (L2)")):
            ea = abs(a[q] - Gex) if q == "G" else a[q]
            eb = abs(b[q] - Gex) if q == "G" else b[q]
            res.ge("G1 observed order of the %s, h = 1/%d to 1/%d" % (lab, round(1/a["h"]), round(1/b["h"])),
                   np.log(ea/eb)/np.log(a["h"]/b["h"]), 1.5)

    # G2
    res.le("G2 runs not converged", sum(not r["converged"] for r in G2), 0,
           " ".join("%s/L%g/%d" % (r["kind"], r["L"], r["dh"]) for r in G2 if not r["converged"]))
    for kind, lab in (("rigid", "drag of the translating rigid sphere / 6 pi mu R U"),
                      ("stressFree", "drag of the stress-free sphere / 4 pi mu R U"),
                      ("rotating", "torque of the rotating sphere / 8 pi mu R^3 Omega")):
        ext = {}
        for dh in dhs:
            rs = sorted([r for r in G2 if r["kind"] == kind and r["dh"] == dh], key=lambda r: r["L"])
            if len(rs) < 3:
                continue
            ext[dh] = g2_fit(rs, kind)
            M["g2_%s_dh%d" % (kind, dh)] = ext[dh]
            res.data.setdefault("g2_ext", {})["%s_dh%d" % (kind, dh)] = ext[dh]
        if not ext:
            continue
        note = "; ".join("d/h %d: %.5f (boxes %s)" % (dh, ext[dh], " ".join(
            "%.4f" % r["ratio"] for r in sorted([r for r in G2 if r["kind"] == kind and r["dh"] == dh],
                                                key=lambda r: r["L"]))) for dh in ext)
        judge = 32 if 32 in ext else max(ext)
        if kind == "stressFree":
            res.report("G2 %s, d/h %d, array correction removed" % (lab, judge), ext[judge], 1.0, note)
        else:
            res.rel("G2 %s, d/h %d, array correction removed" % (lab, judge), ext[judge], 1.0, 0.01, note)
        if 48 in ext and 32 in ext:
            res.le("G2 %s: error at d/h 48 against d/h 32" % kind, abs(ext[48] - 1), abs(ext[32] - 1))
    if ctx.gates_only:
        return
    run_cases(ctx, res)
