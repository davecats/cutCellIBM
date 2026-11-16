"""V10 SPH: rigid sphere D = 1 in a uniform stream U = 1, Re = 20 ... 300 (Johnson & Patel 1999;
Magnaudet, Rivero & Fabre 1995; references in validation/reference_sphere.py, specification
doc/validationSuite.md §10).

Half domain [-10, 20] x [-10, 10] x [0, 10], sphere centre at the origin, symmetry plane z = 0
through it (the plane of the Re 250/300 wake: the lateral force is along y), lateral boundaries
slip, inflow U e_x, outflow zero gradient / p = 0. Mesh family D/h = 24, 32, 48, 64: background
h0 = 16/(D/h) (cubic cells), four nested boxes fixed in space, one level each, the finest the box
[-1, 3] x [-1.2, 1.2] x [0, 1.2] of the spec; each box ~7 cells of its coarser level inside the
next (D/h 24), so the family is geometrically similar.

Steady ibmSimpleFoam (SIMPLEC as V5 2D-1) for Re 20, 50, 100 on D/h 24, 32, 48 and Re 150, 200
on 32, 48, 64; Re 250 (steady, not axisymmetric) on 32, 48, 64 from the Re 200 solution with a
small lateral velocity; Re 300 ibmPimpleFoam from the Re 250 solution, dt = 0.01, statistics
over the last 15 shedding periods. Sampling: alpha-weighted quadratic least squares over the
value points, the cells mirrored at z = 0 (bubble length, vortex centre); the separation angle
from the wall shear the solver applies (separation()). Re_c of the non-axisymmetric wake from
|C_L|^2 at Re 230 and 250 on the finest mesh (reported).

The D/h 64 chain ran on HoreKa in
2026-11 (build of 2190505 against its v2412 module, n16 V14 3D-3Z bit-identical, G1, G2 PASS)."""
import os
import shutil

import numpy as np

import common
import geometry
import meshes
import reference_sphere as ref
from common import write_dict
from foamio import latest_time, read_field

ID, KIND, TIER = "V10", "case", 2
TITLE = "SPH: rigid sphere in a uniform stream, Re 20-300"
STRESSES, RELIES = ["K1", "K3", "K7", "K8"], ["G1", "G2", "G3", "V2"]
REGRESSION = (1e-4,)

R = 0.5
CYL = common.TUTORIALS / "cylinderRe1"
DOM = ((-10.0, -10.0, 0.0), (20.0, 10.0, 10.0))
BOXES = [((-5.2, -5.4, 0.0), (11.0, 5.4, 5.4)),     # level 1
         ((-2.8, -3.0, 0.0), (6.5, 3.0, 3.0)),      # level 2
         ((-1.6, -1.8, 0.0), (4.2, 1.8, 1.8)),      # level 3
         ((-1.0, -1.2, 0.0), (3.0, 1.2, 1.2))]      # level 4, D/h at the sphere
NP = {16: 8, 24: 16, 32: 32, 48: 64, 64: 128}          # ranks (0.27, 0.64, 2.2, 5.1 M cells)
LOW, HIGH = [24, 32, 48], [32, 48, 64]
RES_LOW, RES_HIGH = [20, 50, 100], [150, 200]
DT, END300, WINDOW = 0.01, 300.0, 15/0.137       # Re 300: statistics over the last 15 periods
# the C_D references: Magnaudet et al. (table II of Blanco & Magnaudet) where given, else J&P fig. 5
CD_REF = {Re: ("Magnaudet et al. 1995", ref.CD_TABLE["Magnaudet et al. 1995"][Re]) for Re in (20, 50, 100)}


def setup(d, dh, Re, unsteady=False, lateral=None, full=False, **ibm):
    """case d at D/h = dh; lateral: the half width of the domain (default 10, the outer-boundary
    check 15); full: the whole domain z in [-10, 10] (no symmetry plane, the boxes mirrored);
    ibm: extra top-level ibmProperties entries (ibmConvection, mergeAlpha)"""
    lo, hi = np.array(DOM[0]), np.array(DOM[1])
    if lateral:
        lo[1], hi[1], hi[2] = -lateral, lateral, lateral
    h0 = 16.0/dh
    patches = {"inlet": ("patch", ["xmin"]), "outlet": ("patch", ["xmax"]),
               "lateral": ("patch", ["ymin", "ymax", "zmax"]), "sym": ("symmetryPlane", ["zmin"])}
    boxes = BOXES
    if full:
        lo[2] = -hi[2]
        patches = {"inlet": ("patch", ["xmin"]), "outlet": ("patch", ["xmax"]),
                   "lateral": ("patch", ["ymin", "ymax", "zmin", "zmax"])}
        boxes = [((a[0], a[1], -b[2]), b) for a, b in BOXES]
    ncells = meshes.refined(d, lo, hi, h0, boxes, patches, two_d=False)
    meshes.face_centres(d)
    geometry.ibm_properties(d, {"sphere": {"type": "sphere", "centre": (0, 0, 0), "radius": R,
                                           "wallType": "noSlip"}}, **ibm)
    for x in ("transportProperties", "turbulenceProperties"):
        (d / "constant" / x).write_text((CYL / "constant" / x).read_text())
    common.set_keys(d / "constant" / "transportProperties", nu=1.0/Re)
    meshes.write_fields(d, {
        "U": ("volVectorField", "[0 1 -1 0 0 0 0]", "uniform (1 0 0)",
              {"inlet": "type fixedValue; value uniform (1 0 0);", "lateral": "type slip;",
               "outlet": "type zeroGradient;", **({} if full else {"sym": "type symmetryPlane;"})}),
        "p": ("volScalarField", "[0 2 -2 0 0 0 0]", "uniform 0",
              {"inlet": "type zeroGradient;", "lateral": "type zeroGradient;",
               "outlet": "type fixedValue; value uniform 0;", **({} if full else {"sym": "type symmetryPlane;"})})},
        two_d=False)
    t = (CYL / "system" / "fvSchemes").read_text()
    if unsteady:
        (d / "system" / "fvSchemes").write_text(t.replace("default steadyState;", "default backward;"))
        write_dict(d / "system" / "fvSolution", "dictionary", "fvSolution", """
solvers { p { solver PCG; preconditioner { preconditioner linkedGroup; base { preconditioner GAMG;
                smoother DICGaussSeidel; nCellsInCoarsestLevel 10; tolerance 0; relTol 0.01; } }
              tolerance 1e-9; relTol 0.01; maxIter 1000; }
          pFinal { $p; relTol 0; }
          U { solver PBiCGStab; preconditioner { preconditioner linkedGroup; base DILU; } tolerance 1e-10; relTol 0; }
          UFinal { $U; } }
PIMPLE { momentumPredictor yes; nOuterCorrectors 3; nCorrectors 2; nNonOrthogonalCorrectors 0;
         consistent yes; ddtCorr yes; }
""")
    else:
        (d / "system" / "fvSchemes").write_text(t)
        write_dict(d / "system" / "fvSolution", "dictionary", "fvSolution", """
solvers { p { solver PCG; preconditioner { preconditioner linkedGroup; base { preconditioner GAMG;
                smoother DICGaussSeidel; nCellsInCoarsestLevel 10; tolerance 0; relTol 0.01; } }
              tolerance 1e-12; relTol 0.01; maxIter 1000; }
          U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-12; relTol 0.1; } }
SIMPLE { nNonOrthogonalCorrectors 0; consistent yes; residualControl { p 1e-9; U 1e-9; } }
relaxationFactors { equations { U 0.9; } fields { p 1.0; } }
""")
        geometry.control_dict(d, "ibmSimpleFoam", endTime=20000, deltaT=1, writeInterval=20000)
    return ncells


def decomposed(d, app, np_, extra=""):
    write_dict(d / "system" / "decomposeParDict", "dictionary", "decomposeParDict",
               "numberOfSubdomains %d;\nmethod scotch;\n%s" % (np_, common.LINK_CONSTRAINT))
    common.link_zone(d)
    common.sh("decomposePar -force" + extra, d, "log.decomposePar")
    common.sh(app, d, "log." + app, np_=np_)
    common.sh("reconstructPar -latestTime", d, "log.reconstructPar")
    common.sh("ibmSetGeometry", d, "log.ibmSetGeometry")


def start_from(d, src, perturb=0.0):
    """initial U, p of case d from the latest time of src (the same mesh); perturb: a lateral
    velocity perturb*e_y added in the near wake (0.5 < x < 3, r < 1), to leave the axisymmetric
    branch"""
    t = latest_time(str(src))
    for f in ("U", "p"):
        shutil.copy(src / t / f, d / "0" / f)
    if perturb:
        C, _ = read_field(os.path.join(d, "0", "C"))
        U, _ = read_field(os.path.join(d, "0", "U"))
        r = np.hypot(C[:, 1], C[:, 2])
        m = (C[:, 0] > 0.5) & (C[:, 0] < 3) & (r < 1)
        U[m, 1] += perturb
        txt = (d / "0" / "U").read_text()
        i0 = txt.index("internalField")
        i1 = txt.index("boundaryField")
        txt = txt[:i0] + "internalField %s;\n" % meshes.vector_list(U) + txt[i1:]
        (d / "0" / "U").write_text(txt)


class Field:
    """value points, alpha, U, p of the fluid near the sphere, mirrored at z = 0, with a k-d tree"""

    def __init__(self, d, box=((-1.5, -2.0, 0.0), (4.0, 2.0, 2.0))):
        from scipy.spatial import cKDTree
        t = latest_time(str(d))
        xc, a = meshes.fluid_cells(str(d), t)
        U, _ = read_field(os.path.join(d, t, "U"))
        p, _ = read_field(os.path.join(d, t, "p"))
        m = (a > 0) & np.all(xc >= box[0], axis=1) & np.all(xc <= box[1], axis=1)
        xc, a, U, p = xc[m], a[m], U[m], p[m]
        xm, Um = xc*(1, 1, -1), U*(1, 1, -1)
        self.x, self.a = np.vstack([xc, xm]), np.concatenate([a, a])
        self.f = np.column_stack([np.vstack([U, Um]), np.concatenate([p, p])])
        self.tree = cKDTree(self.x)

    def fit(self, x0, h, radius=4.0):
        """values (u, v, w, p) at x0 and their gradients (4 x 3) by a weighted quadratic fit"""
        idx = self.tree.query_ball_point(x0, radius*h)
        X = (self.x[idx] - x0)/h
        cols = [np.ones(len(X))] + [X[:, i] for i in range(3)] + [X[:, i]*X[:, j] for i in range(3) for j in range(i, 3)]
        w = np.sqrt(self.a[idx])[:, None]
        c, *_ = np.linalg.lstsq(np.column_stack(cols)*w, self.f[idx]*w, rcond=None)
        return c[0], c[1:4].T/h


def separation(d, t, bands=(2.0, 1.0, 3.0)):
    """separation angle from the front stagnation point (deg): the sign change of the meridional
    wall shear the solver applies, nu c V u_P / |S_w| per wall cell (the viscous force of the
    C_D), area-averaged in bands of theta (all meridians: the flow is axisymmetric), linearly
    interpolated between the bands. The wall cells cluster on the grid planes x = const (a band
    spacing of ~h/(R sin theta)), so the first band width is the result, the others its spread.
    Least-squares fits of u_theta at fixed wall distances (h, 2h) were biased by -1 to -1.5 deg
    and did not converge (Re 50, 100 at D/h 24-48); 180 if the shear does not change sign."""
    if not (d / "0" / "V").exists():
        common.sh("postProcess -func writeCellVolumes -time 0", d, "log.writeCellVolumes")
    rf = lambda *p: read_field(os.path.join(d, *p))[0]
    C, V, c, Sw, U = rf("0", "C"), rf("0", "V"), rf("0", "ibmNoSlipCoeff"), rf("0", "ibmSw"), rf(t, "U")
    A = np.linalg.norm(Sw, axis=1)
    w = (A > 0) & (c > 0)
    x, A = C[w], A[w]
    r = np.linalg.norm(x, axis=1)
    th = np.arccos(np.clip(-x[:, 0]/r, -1, 1))
    rho = np.maximum(np.hypot(x[:, 1], x[:, 2]), 1e-12)
    e = np.column_stack([np.sin(th), np.cos(th)*x[:, 1]/rho, np.cos(th)*x[:, 2]/rho])
    f = c[w]*V[w]*np.sum(U[w]*e, axis=1)
    out = []
    for bw in bands:
        k = np.digitize(th, np.radians(np.arange(90.0, 180.01, bw)))
        ks = [i for i in np.unique(k) if A[k == i].sum() > 0]
        T = np.array([th[k == i] @ A[k == i]/A[k == i].sum() for i in ks])
        S = np.array([f[k == i].sum()/A[k == i].sum() for i in ks])
        z = np.where((S[:-1] > 0) & (S[1:] < 0))[0]
        out.append(180.0 if len(z) == 0 else
                   float(np.degrees(T[z[-1]] + (T[z[-1]+1] - T[z[-1]])*S[z[-1]]/(S[z[-1]] - S[z[-1]+1]))))
    return out


def bubble(F, h):
    """x_s: length of the recirculation on the axis behind the sphere (from x = R), the last
    sign change of u_x from - to + on (x, 0, 0)"""
    x = np.arange(R + h, R + 3.0, h/2)
    u = np.array([F.fit((xi, 0, 0), h)[0][0] for xi in x])
    k = np.where((u[:-1] < 0) & (u[1:] >= 0))[0]
    if len(k) == 0:
        return 0.0
    k = k[-1]
    return float(x[k] + (x[k+1] - x[k])*u[k]/(u[k] - u[k+1]) - R)


def vortex(F, h, x0):
    """centre (x_c, y_c) of the ring vortex in the plane z = 0: Newton on (u_x, u_y) = 0"""
    x = np.array([x0[0], x0[1], 0.0])
    for _ in range(20):
        v, g = F.fit(x, h)
        J = g[:2, :2]
        dx = np.linalg.solve(J, -v[:2])
        x[:2] += dx
        if np.linalg.norm(dx) < 1e-6:
            break
    return float(x[0]), float(x[1])


def forces(d, full=False):
    """C_D, C_L (lateral, y) of the whole sphere: the (half-domain: twice the) force over pi D^2/8"""
    log = (d / "log.ibmSimpleFoam").read_text() if (d / "log.ibmSimpleFoam").exists() else None
    t, Fh = meshes.force_history(d, 1.0)
    c = (1 if full else 2)*8/np.pi
    return t, c*Fh[:, 0], c*Fh[:, 1], log


def steady(ctx, dh, Re, tag="", src=None, perturb=0.0, **opts):
    d = ctx.root / ("re%d_dh%d%s" % (Re, dh, tag))
    ncells = setup(d, dh, Re, **opts)
    if src is not None:
        start_from(d, ctx.root / src, perturb)
    decomposed(d, "ibmSimpleFoam", NP[dh]*(2 if opts.get("full") else 1))
    t, cd, cl, log = forces(d, opts.get("full", False))
    h = 1.0/dh
    out = {"Re": Re, "Dh": dh, "tag": tag, "cells": ncells, "C_D": float(cd[-1]), "C_L": float(cl[-1]),
           "iterations": int(t[-1]), "converged": "SIMPLE solution converged" in log,
           "C_Dchange": float(np.ptp(cd[-200:])/abs(cd[-1]))}
    if Re <= 200 and not opts.get("full"):
        F = Field(d)
        th = separation(d, latest_time(str(d)))
        out.update({"theta_s": th[0], "theta_s_spread": float(np.ptp(th)), "x_s": bubble(F, h)})
        s = ref.steady()
        if Re in s["Re"]:
            i = list(s["Re"]).index(Re)
            if np.isfinite(s["x_c"][i]):
                out["x_c"], out["y_c"] = vortex(F, h, (s["x_c"][i], s["y_c"][i]))
    return out


def unsteady(ctx, dh, Re=300, dtfactor=1.0, src=None):
    d = ctx.root / ("re%d_dh%d_dt%g" % (Re, dh, dtfactor))
    ncells = setup(d, dh, Re, unsteady=True)
    dt = DT*dtfactor
    steps = int(round(END300/dt))
    geometry.control_dict(d, "ibmPimpleFoam", endTime=END300, deltaT=repr(dt), writeInterval=steps)
    if src is not None:
        start_from(d, ctx.root / src)
    decomposed(d, "ibmPimpleFoam", NP[dh])
    t, cd, cl, _ = forces(d)
    log = (d / "log.ibmPimpleFoam").read_text()
    import re
    co = max(float(x) for x in re.findall(r"Courant Number mean: \S+ max: (\S+)", log))
    st = meshes.periodic_statistics(t, cl, cd, WINDOW)
    return {"Re": Re, "Dh": dh, "dtfactor": dtfactor, "cells": ncells, "Co": co, "stats": st,
            "series": [t[::20].tolist(), cd[::20].tolist(), cl[::20].tolist()]}


def chain(ctx, dh, quick=False):
    """Re 200 -> 250 (perturbed) -> 300 on one mesh; Re 230 (onset) on the finest"""
    out = [("steady", steady(ctx, dh, 200))]
    out.append(("steady", steady(ctx, dh, 250, src="re200_dh%d" % dh, perturb=0.01)))
    if dh == HIGH[-1]:
        out.append(("steady", steady(ctx, dh, 230, src="re200_dh%d" % dh, perturb=0.01)))
    if not quick:
        out.append(("unsteady", unsteady(ctx, dh, src="re250_dh%d" % dh)))
        if dh == HIGH[0]:
            out.append(("unsteady", unsteady(ctx, dh, dtfactor=0.5, src="re250_dh%d" % dh)))
    return out


def run(ctx, res):
    low, high = (LOW[:1], HIGH[:1]) if ctx.quick else (LOW, HIGH)
    mid = LOW[1]
    jobs = [(chain, dh, (dh, ctx.quick)) for dh in high]
    jobs += [(steady, dh, (dh, Re)) for dh in low for Re in RES_LOW]
    jobs += [(steady, dh, (dh, Re)) for dh in high for Re in RES_HIGH if Re != 200]
    if not ctx.quick:   # once per case: half against full domain, outer boundary, the IBM options
        jobs += [(steady, mid, (mid, 50), {"tag": "_full", "full": True}),
                 (steady, mid, (mid, 50), {"tag": "_lat15", "lateral": 15.0}),
                 (steady, mid, (mid, 100), {"tag": "_noconv", "ibmConvection": False}),
                 (steady, mid, (mid, 100), {"tag": "_merge0", "mergeAlpha": 0.0})]
    jobs = [j + ({},) if len(j) == 3 else j for j in jobs]
    import threading
    need = lambda j: NP[j[1]]*(2 if j[3].get("full") else 1)
    free, cv = [max(ctx.jobs, max(need(j) for j in jobs))], threading.Condition()

    def one(j):
        with cv:
            cv.wait_for(lambda: free[0] >= need(j))
            free[0] -= need(j)
        try:
            r = j[0](ctx, *j[2], **j[3])
            return r if j[0] is chain else [("steady", r)]
        finally:
            with cv:
                free[0] += need(j)
                cv.notify_all()
    jobs.sort(key=lambda j: (j[0] is chain, need(j)))   # steady cases smallest first, the long chains last
    out = [x for r in common.pmap(one, jobs, len(jobs)) for x in r]
    st = sorted([r for k, r in out if k == "steady"], key=lambda r: (r["tag"], r["Re"], r["Dh"]))
    un = sorted([r for k, r in out if k == "unsteady"], key=lambda r: (r["Dh"], -r["dtfactor"]))
    res.data["steady"], res.data["unsteady"] = st, un
    base = [r for r in st if not r["tag"]]
    M = res.data["metrics"] = {}
    for r in st:
        for q in ("C_D", "C_L", "theta_s", "x_s", "x_c", "y_c"):
            if q in r:
                M["re%d_dh%d%s_%s" % (r["Re"], r["Dh"], r["tag"], q)] = r[q]
    for r in un:
        for q in ("C_D", "C_L", "C_Lrms", "C_Drms", "f"):
            M["re300_dh%d_dt%g_%s" % (r["Dh"], r["dtfactor"], q)] = r["stats"][q]
    res.le("steady runs not converged", sum(not r["converged"] for r in st), 0,
           " ".join("%d/%d:%.1e" % (r["Re"], r["Dh"], r["C_Dchange"]) for r in st if not r["converged"]))
    s = ref.steady()
    jp = {int(Re): i for i, Re in enumerate(s["Re"])}
    for Re in RES_LOW + RES_HIGH:
        rs = [r for r in base if r["Re"] == Re]
        if not rs:
            continue
        f, seq = rs[-1], lambda q: " ".join("%.5g" % r[q] for r in rs)
        name, cref = CD_REF.get(Re, ("Johnson & Patel 1999", s["C_D"][jp[Re]] if Re in jp else np.nan))
        res.rel("Re %d C_D, D/h %d (%s)" % (Re, f["Dh"], name), f["C_D"], cref, 0.015, seq("C_D"))
        if len(rs) >= 3:
            Ce, p = common.richardson([1/r["Dh"] for r in rs], [r["C_D"] for r in rs])
            res.report("Re %d observed order of C_D" % Re, p,
                       note="extrapolated %.5g (%s %.5g, %+.2f %%)" % (Ce, name, cref, 100*(Ce/cref - 1)))
        if Re in jp and Re >= 50:
            i = jp[Re]
            res.absolute("Re %d separation angle theta_s (deg), D/h %d" % (Re, f["Dh"]), f["theta_s"],
                      s["theta_s"][i], 1.5, seq("theta_s") + "; spread over band widths 1-3 deg %.2f" % f["theta_s_spread"])
            res.absolute("Re %d bubble length x_s / D, D/h %d" % (Re, f["Dh"]), f["x_s"], s["x_s"][i], 0.03, seq("x_s"))
            res.absolute("Re %d vortex centre x_c / D, D/h %d" % (Re, f["Dh"]), f["x_c"], s["x_c"][i], 0.02, seq("x_c"))
            res.absolute("Re %d vortex centre y_c / D, D/h %d" % (Re, f["Dh"]), f["y_c"], s["y_c"][i], 0.02, seq("y_c"))
        elif Re == 20:
            res.report("Re 20 separation angle theta_s (deg; J&P: onset at Re 20), D/h %d" % f["Dh"], f["theta_s"],
                       note=seq("theta_s") + "; x_s " + seq("x_s"))
    ReL, CL = ref.lateral()
    r250 = [r for r in base if r["Re"] == 250]
    if r250:
        f = r250[-1]
        cref = float(CL[list(ReL).index(250)])
        res.rel("Re 250 |C_L| (steady, not axisymmetric), D/h %d" % f["Dh"], abs(f["C_L"]), cref, 0.05,
                " ".join("%.5g" % abs(r["C_L"]) for r in r250))
        res.report("Re 250 C_D, D/h %d" % f["Dh"], f["C_D"], note=" ".join("%.5g" % r["C_D"] for r in r250))
        r230 = [r for r in base if r["Re"] == 230]
        if r230:   # Landau: C_L^2 linear in Re - Re_c
            a, b = abs(r230[-1]["C_L"]), abs(f["C_L"])
            rec = 230 - a*a*(250 - 230)/(b*b - a*a)
            res.report("onset of the non-axisymmetric wake Re_c (from |C_L|^2 at Re 230, 250; J&P 211, "
                       "Natarajan & Acrivos 210)", rec, note="|C_L| %.4g (J&P %.4g) at Re 230" % (
                           a, abs(float(CL[list(ReL).index(230)]))))
    u1 = [r for r in un if r["dtfactor"] == 1.0]
    if u1:
        g, B = u1[-1], ref.RE300
        sq = lambda q: " ".join("%.5g" % r["stats"][q] for r in u1)
        res.report("Re 300 max Courant number", max(r["Co"] for r in un), note=" ".join("%.3f" % r["Co"] for r in un))
        res.rel("Re 300 mean C_D, D/h %d" % g["Dh"], g["stats"]["C_D"], B["C_D mean"], 0.02, sq("C_D"))
        res.rel("Re 300 St, D/h %d" % g["Dh"], g["stats"]["f"], B["St"], 0.02, sq("f"))
        res.rel("Re 300 mean |C_L|, D/h %d" % g["Dh"], abs(g["stats"]["C_L"]), abs(B["C_L mean"]), 0.10, sq("C_L"))
        res.report("Re 300 C_L amplitude (sqrt 2 rms), D/h %d" % g["Dh"], np.sqrt(2)*g["stats"]["C_Lrms"],
                   B["C_L amplitude"], " ".join("%.4g" % (np.sqrt(2)*r["stats"]["C_Lrms"]) for r in u1))
        h = [r for r in un if r["Dh"] == HIGH[0]]
        if len(h) == 2:
            a, b = h
            res.report("Re 300 dt/2 against dt at D/h %d, mean C_D, St (rel.)" % a["Dh"],
                       abs(b["stats"]["C_D"]/a["stats"]["C_D"] - 1),
                       note="C_D %+.2e, St %+.2e" % (b["stats"]["C_D"]/a["stats"]["C_D"] - 1,
                                                     b["stats"]["f"]/a["stats"]["f"] - 1))
    ex = {r["tag"]: r for r in st if r["tag"]}
    for tag, Re, what in (("_full", 50, "full domain against the half domain (symmetry plane through the cut cells)"),
                          ("_lat15", 50, "lateral boundaries at 15 D against 10 D"),
                          ("_noconv", 100, "ibmConvection off against on"),
                          ("_merge0", 100, "mergeAlpha 0 against 0.5")):
        if tag in ex:
            r0 = [r for r in base if r["Re"] == Re and r["Dh"] == mid][0]
            res.report("Re %d C_D, D/h %d: %s (rel.)" % (Re, mid, what), ex[tag]["C_D"]/r0["C_D"] - 1,
                       note="%.6g against %.6g" % (ex[tag]["C_D"], r0["C_D"]))
