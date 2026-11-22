"""V15 TMR 2DZP: turbulent flat plate, Spalart-Allmaras, Re = 5e6 per unit length (NASA TMR;
reference_tmr.py; doc/validationSuite.md §11). Body-fitted (simpleFoam + SpalartAllmaras on the
TMR grid) against immersed (the same grid lines with one cell row added below the plate and the
grid shifted by f dy1: the plate, a box body, cuts that row with fluid fraction f). The body-fitted
run is the reference; the TMR values check it (REPORT: OpenFOAM's SA is not the TMR SA).
Tier 1: grid 273 x 193, wall-resolved, f = 0.5.
Tier 2 (--tier2, doc/validationSuite.md §11): the TMR grids 137, 273, 545, f = 0.25, 0.5, 0.75,
wall law SA and none; options on 273 (mergeAlpha 0 at f 0.25, ibmConvection off,
ibmConvectionWallLaw on, a freeSlip box ahead of the plate); wall functions on the 273 grid
coarsened at the wall to a first cell of h+ 30 and 90 at x = 0.97 (as V7), body-fitted with
Spalding's law at the same first-cell height. C_f(x) of the immersed plate: 2 ibmWallShear, the
traction the wall cells apply (per cell the term whose sum is the plate force)."""
import gzip
import os

import numpy as np

import common
import reference_tmr as tmr
from common import write_dict
from foamio import latest_time, read_field

ID, KIND, TIER = "V15", "case", 1
TITLE = "TMR flat plate 2DZP, SA, immersed vs body-fitted (wall-resolved)"
STRESSES, RELIES = ["K12", "K13"], ["V7"]
REGRESSION = (1e-4,)

NU = 1/5e6
GRIDS = common.RUNS / "tmr_grids" / "zip"   # flatplategrids.zip of the TMR, extracted
GRID = {137: "flatplate_clust2_2levelsdown_137x97.p2dfmt", 273: "flatplate_clust2_1leveldown_273x193.p2dfmt",
        545: "flatplate_clust2.p2dfmt"}
DZ = 0.01
XCF = tmr.FLATPLATE["x_cf"]
XPROF = (0.97008, 1.90334)                    # TMR profile stations
UTAU097 = (0.002729/2)**0.5                   # TMR incompressible C_f(0.97): the h+ of the wall-function meshes
HPLUS = (30, 90)
GRIDS2 = (137, 273, 545)
FS = (0.25, 0.5, 0.75)
XWIN = (0.1, 1.9)                             # C_f(x) window: the leading edge excluded
RTOL_WR, RTOL_CFX, RTOL_NUT, RTOL_SPREAD, RTOL_WF = 0.01, 0.02, 0.02, 0.01, 0.03
END = {545: 80000}                            # SIMPLE iteration cap (residualControl ends the runs):
                                              # the 545 grid needs ~42000 (Uy 1.5e-8 at 40000, HK2 2026-11-23)


def read_grid(n):
    p = GRIDS / GRID[n]
    # read the .gz directly: cases set up in parallel threads must not race on unpacking it
    tok = (p.read_text() if p.exists() else gzip.open(str(p) + ".gz", "rt").read()).split()
    ni, nj = int(tok[1]), int(tok[2])
    v = np.array(tok[3:3 + 2*ni*nj], float)
    x, y = v[:ni*nj].reshape(nj, ni), v[ni*nj:].reshape(nj, ni)
    return x[0], y[:, 0]                      # tensor product: x_i, y_j


def grid_lines(g):
    """x_i, y_j of a TMR grid (int) or of a wall-function mesh "h30", "h90": the 273 grid with its
    lines below the first one spaced >= h replaced by a uniform layer ending on that line, cells of
    height ~h (within 3 %), h+ = h u_tau/nu at x = 0.97"""
    if not isinstance(g, str):
        return read_grid(g)
    xs, ys = read_grid(273)
    h = int(g[1:])*NU/UTAU097
    j = int(np.argmax(np.diff(ys) >= h))      # first spacing >= h
    k = max(int(round(ys[j]/h)), 1)
    return xs, np.concatenate([np.linspace(0.0, ys[j], k + 1), ys[j + 1:]])


def write_p2d(path, xs, ys):
    X, Y = np.meshgrid(xs, ys)
    with open(path, "w") as fh:
        fh.write("1\n%d %d\n" % (len(xs), len(ys)))
        for a in (X.ravel(), Y.ravel()):
            for i in range(0, len(a), 4):
                fh.write(" ".join("%.15e" % v for v in a[i:i + 4]) + "\n")


def mesh(d, n, ibm, f):
    xs, ys = grid_lines(n)
    if ibm:
        dy1 = ys[1] - ys[0]
        ys = np.concatenate([[ys[0] - dy1], ys]) + f*dy1
    (d / "system").mkdir(parents=True, exist_ok=True)
    write_p2d(d / "grid.p2dfmt", xs, ys)
    import geometry
    geometry.control_dict(d, "simpleFoam")
    y0 = float(ys[0])
    write_dict(d / "system" / "topoSetDict", "dictionary", "topoSetDict", """
actions
(
    { name bnd;    type faceSet; action new; source patchToFace; patch defaultFaces; }
    { name fb;     type faceSet; action new; source normalToFace; normal (0 0 1); cos 0.9; }
    { name fb;     type faceSet; action add; source normalToFace; normal (0 0 -1); cos 0.9; }
    { name fb;     type faceSet; action subset; source faceToFace; set bnd; }
    { name inlet;  type faceSet; action new; source boxToFace; box (-0.4 -1 -1) (-0.3333 2 1); }
    { name inlet;  type faceSet; action subset; source faceToFace; set bnd; }
    { name outlet; type faceSet; action new; source boxToFace; box (1.9999 -1 -1) (2.1 2 1); }
    { name outlet; type faceSet; action subset; source faceToFace; set bnd; }
    { name top;    type faceSet; action new; source boxToFace; box (-1 0.9999 -1) (3 2 1); }
    { name top;    type faceSet; action subset; source faceToFace; set bnd; }
    { name sym;    type faceSet; action new; source boxToFace; box (-1 %r -1) (%r %r 1); }
    { name sym;    type faceSet; action subset; source faceToFace; set bnd; }
    { name sym;    type faceSet; action subtract; source faceToFace; set fb; }
    { name plate;  type faceSet; action new; source boxToFace; box (0 %r -1) (2.1 %r 1); }
    { name plate;  type faceSet; action subset; source faceToFace; set bnd; }
    { name plate;  type faceSet; action subtract; source faceToFace; set fb; }
    { name plate;  type faceSet; action subtract; source faceToFace; set sym; }
);
""" % (y0 - 1e-9, 0.0 if not ibm else 3.0, y0 + 1e-9, y0 - 1e-9, y0 + 1e-9))
    pl = '    { name plate; patchInfo { type wall; }; constructFrom set; set plate; }\n' if not ibm else ""
    write_dict(d / "system" / "createPatchDict", "dictionary", "createPatchDict", """
pointSync false;
patches
(
    { name inlet;  patchInfo { type patch; }; constructFrom set; set inlet; }
    { name outlet; patchInfo { type patch; }; constructFrom set; set outlet; }
    { name top;    patchInfo { type patch; }; constructFrom set; set top; }
    { name sym;    patchInfo { type symmetryPlane; }; constructFrom set; set sym; }
%s    { name frontAndBack; patchInfo { type empty; }; constructFrom set; set fb; }
);
""" % pl)
    return xs, ys


def build_mesh(d):
    """after all dictionaries: createPatch reads system/fvSchemes from v2412 on"""
    common.sh("plot3dToFoam -2D %r -noBlank grid.p2dfmt" % DZ, d, "log.plot3dToFoam")
    common.sh("topoSet", d, "log.topoSet")
    common.sh("createPatch -overwrite", d, "log.createPatch")


def case(d, n, ibm, f, ma=0.5, law="spalartAllmaras", opts=None, lead=False):
    end = END.get(n, 40000)
    xs, ys = mesh(d, n, ibm, f)
    nut_inf = 3*NU*27/(27 + 7.1**3)          # nu_t = nuTilda fv1, chi = 3
    bc = {"inlet": {"U": "fixedValue; value uniform (1 0 0);", "p": "zeroGradient;",
                    "nuTilda": "fixedValue; value uniform %r;" % (3*NU), "nut": "calculated; value uniform %r;" % nut_inf},
          "outlet": {"U": "zeroGradient;", "p": "fixedValue; value uniform 0;", "nuTilda": "zeroGradient;",
                     "nut": "calculated; value uniform %r;" % nut_inf},
          "top": {"U": "slip;", "p": "slip;", "nuTilda": "slip;", "nut": "calculated; value uniform %r;" % nut_inf},
          "sym": {k: "symmetryPlane;" for k in ("U", "p", "nuTilda", "nut")},
          "frontAndBack": {k: "empty;" for k in ("U", "p", "nuTilda", "nut")}}
    if not ibm:
        bc["plate"] = {"U": "noSlip;", "p": "zeroGradient;", "nuTilda": "fixedValue; value uniform 0;",
                       "nut": "nutUSpaldingWallFunction; value uniform 0;"}
    for name, cls, dim, val in (("U", "volVectorField", "[0 1 -1 0 0 0 0]", "(1 0 0)"),
                                ("p", "volScalarField", "[0 2 -2 0 0 0 0]", "0"),
                                ("nuTilda", "volScalarField", "[0 2 -1 0 0 0 0]", repr(3*NU)),
                                ("nut", "volScalarField", "[0 2 -1 0 0 0 0]", repr(nut_inf))):
        write_dict(d / "0" / name, cls, name, "dimensions %s;\ninternalField uniform %s;\nboundaryField\n{\n%s}\n"
                   % (dim, val, "".join("    %s { type %s }\n" % (p, b[name]) for p, b in bc.items())))
    write_dict(d / "constant" / "transportProperties", "dictionary", "transportProperties",
               "transportModel Newtonian;\nnu %r;\n" % NU)
    write_dict(d / "constant" / "turbulenceProperties", "dictionary", "turbulenceProperties",
               "simulationType RAS;\nRAS { RASModel %s; turbulence on; printCoeffs on; }\n"
               % ("ibmSpalartAllmaras" if ibm else "SpalartAllmaras"))
    if ibm:
        dxle = float(xs[np.argmin(np.abs(xs)) + 1] - xs[np.argmin(np.abs(xs))])
        import geometry
        # leading edge half a cell downstream of the mesh line x = 0 (a body face on a mesh line
        # is degenerate for the level set): Re_x at x = 0.97 changes by 0.05 %
        bodies = {"plate": {"type": "box", "min": (0.5*dxle, -1.0, -1.0), "max": (3.0, 0.0, 1.0),
                            "wallType": "noSlip", "wallFunction": law}}
        if lead:
            # variant: a freeSlip box ahead of the plate at its height, touching it (two bodies;
            # their common face is no wall, G1)
            bodies["lead"] = {"type": "box", "min": (-1.0, -1.0, -1.0), "max": (0.5*dxle, 0.0, 1.0),
                              "wallType": "freeSlip"}
        geometry.ibm_properties(d, bodies, mergeAlpha=ma, nSubSamples=32, **(opts or {}))
    write_dict(d / "system" / "fvSchemes", "dictionary", "fvSchemes", """
ddtSchemes      { default steadyState; }
gradSchemes     { default Gauss linear; }
divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,nuTilda) bounded Gauss linearUpwind grad(nuTilda);
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear %s; }
interpolationSchemes { default linear; }
snGradSchemes   { default %s; }
wallDist        { method %s; }
""" % (("cutCellCorrected", "cutCellCorrected", "ibmLevelSet") if ibm else ("corrected", "corrected", "meshWave")))
    fv_solution(d, ibm)
    write_dict(d / "system" / "controlDict", "dictionary", "controlDict", """
application %s;
%sstartFrom startTime; startTime 0; stopAt endTime; endTime %d; deltaT 1;
writeControl timeStep; writeInterval %d; purgeWrite 0; writeFormat ascii; writePrecision 12;
timeFormat general; timePrecision 6; runTimeModifiable false;
""" % ("ibmSimpleFoam" if ibm else "simpleFoam", 'libs ("libcutCellIBMTurbulence.so");\n' if ibm else "", end, end))
    if not ibm:
        write_dict(d / "system" / "controlDict", "dictionary", "controlDict",
                   (d / "system" / "controlDict").read_text().split("\n", 1)[1] + """
functions { wss { type wallShearStress; libs (fieldFunctionObjects); patches (plate);
            writeControl writeTime; } }
""")
    build_mesh(d)


def fv_solution(d, ibm):
    """immersed: CG preconditioned by GAMG with the incomplete-Cholesky smoother for the pressure,
    which resolves the stiff two-cell systems of linked cells (a point Gauss-Seidel smoother
    cannot: plain GAMG stalls) and is as fast without links; U and nuTilda need no special solver
    (linkCells redistributes their group residuals). linkedGroup: the same for linked groups
    split over processors (README section 4); relTol 0.01 in the GAMG dictionary is the
    tolerance of its coarsest level (to round-off costs 40 % of the run time, same iterations).
    Body-fitted: the standard solvers."""
    solvers = ('    p { solver PCG; preconditioner { preconditioner linkedGroup; base { preconditioner GAMG; '
               'smoother DICGaussSeidel; nCellsInCoarsestLevel 10; tolerance 0; relTol 0.01; } } '
               'tolerance 1e-10; relTol 0.05; }\n'
               if ibm else
               '    p { solver GAMG; smoother GaussSeidel; tolerance 1e-10; relTol 0.05; }\n') + (
               '    "(U|nuTilda)" { solver smoothSolver; smoother GaussSeidel; nSweeps 2; tolerance 1e-12; relTol 0.1; }\n')
    write_dict(d / "system" / "fvSolution", "dictionary", "fvSolution", """
solvers
{
%s}
SIMPLE
{
    nNonOrthogonalCorrectors 0; consistent yes;
    residualControl { p 1e-7; U 1e-8; nuTilda 1e-8; }
}
relaxationFactors { fields { p 1; } equations { U 0.7; nuTilda 0.7; } }
""" % solvers)


def cf_bodyfitted(d, t):
    """x, C_f = 2 |tau_w| (U = 1) on the plate faces from the wallShearStress function object"""
    _, pt = read_field(os.path.join(d, t, "wallShearStress"))
    _, pc = read_field(os.path.join(d, "0", "C"))
    return pc["plate"][:, 0], 2*np.abs(pt["plate"][:, 0])


OPTIONS = {"noConv": dict(ibmConvection=False), "convLaw": dict(ibmConvectionWallLaw=True)}


def name_of(n, ibm, f=0.5, ma=0.5, law="spalartAllmaras", opt=None):
    if not ibm:
        return "bf_%s" % n
    return "ibm_%s_f%g_ma%g%s%s" % (n, f, ma, "_noLaw" if law == "none" else "", "_" + opt if opt else "")


def column(x, xc):
    """mask of the cells of the mesh column nearest to x"""
    return np.abs(x - x[np.argmin(np.abs(x - xc))]) < 1e-9


def profiles(C, U, nut, y0, cfx, cf):
    """u+(y+) at the TMR stations (y from the plate, u_tau from the local C_f) and max nu_t/nu at
    x = 0.97"""
    out = {}
    for xs in XPROF:
        m = column(C[:, 0], xs)
        ut = (0.5*float(np.interp(xs, cfx, cf)))**0.5
        o = np.argsort(C[m, 1])
        y = C[m, 1][o] - y0
        out["prof_%g" % xs] = {"yplus": (y*ut/NU).tolist(), "uplus": (U[m, 0][o]/ut).tolist()}
    out["nut097"] = float(np.max(nut[column(C[:, 0], XCF)])/NU)
    return out


def one(ctx, n, ibm, f=0.5, ma=0.5, law="spalartAllmaras", opt=None):
    name = name_of(n, ibm, f, ma, law, opt)
    d = ctx.root / name
    try:
        return _one(ctx, d, name, n, ibm, f, ma, law, opt)
    except common.RunError as e:      # tier 2: a diverged run is a result, not a driver crash
        if ctx.tier2 and name not in TIER1_NAMES:
            return {"name": name, "n": n, "ibm": ibm, "f": f, "ma": ma, "law": law, "opt": opt,
                    "C_D": float("nan"), "converged": False, "error": str(e)[:200]}
        raise


def _one(ctx, d, name, n, ibm, f, ma, law, opt):
    if not (common.REUSE and (d / "system").exists()):
        case(d, n, ibm, f, ma, law, OPTIONS.get(opt), lead=(opt == "lead"))
    (d / "0").mkdir(exist_ok=True)
    common.sh("checkMesh -constant", d, "log.checkMesh")
    common.sh("postProcess -func writeCellCentres -time 0", d, "log.postProcess")
    if ibm:
        common.sh("ibmSetGeometry", d, "log.ibmSetGeometry")
    app = "ibmSimpleFoam" if ibm else "simpleFoam"
    common.sh(app, d, "log." + app)
    log = (d / ("log." + app)).read_text()
    t = latest_time(str(d))
    if ibm:
        # C_D = F_x per unit depth/(rho U^2 L/2), L = 2
        F = common.ibm_force(log, "plate")
        CD, CDp, CDv = (F[k][0]/DZ for k in ("total", "pressure", "viscous"))
        x, cf = cf_ibm(d, t)
    else:
        x, cf = cf_bodyfitted(d, t)
        o = np.argsort(x)
        x, cf = x[o], cf[o]
        dx = np.gradient(x)
        CD = float(np.sum(cf*dx)/2.0)
        CDp, CDv = 0.0, CD
    out = {"name": name, "n": n, "ibm": ibm, "f": f, "ma": ma, "law": law, "opt": opt, "C_D": float(CD),
           "C_Dp": float(CDp), "C_Dv": float(CDv), "converged": "SIMPLE solution converged" in log,
           "iterations": int(t) if t.isdigit() else None}
    if x is not None:
        out["cf097"] = float(np.interp(XCF, x, cf))
        if ctx.tier2:
            C, _ = read_field(os.path.join(d, "0", "C"))
            U, _ = read_field(os.path.join(d, t, "U"))
            nut, _ = read_field(os.path.join(d, t, "nut"))
            if ibm:
                a, _ = read_field(os.path.join(d, "0", "ibmAlpha"))
                off = os.path.join(d, t, "ibmCentroidOffset")
                if os.path.exists(off):
                    C = C + read_field(off)[0]
                m = a > 0
                C, U, nut = C[m], U[m], nut[m]
            out.update(profiles(C, U, nut, 0.0, x, cf))
            out["cfx"], out["cf"] = x.tolist(), cf.tolist()
    return out


def cf_ibm(d, t):
    """x, C_f = 2 ibmWallShear of the plate's wall cells (one per mesh column)"""
    tw = os.path.join(d, t, "ibmWallShear")
    if not os.path.exists(tw):
        return None, None
    C, _ = read_field(os.path.join(d, "0", "C"))
    Sw, _ = read_field(os.path.join(d, "0", "ibmSw"))
    tau, _ = read_field(tw)
    m = (np.sum(Sw**2, axis=1) > 0) & (C[:, 0] > 0) & (tau > 0)
    o = np.argsort(C[m, 0])
    return C[m, 0][o], 2*tau[m][o]


TIER1 = [dict(n=273, ibm=False), dict(n=273, ibm=True, f=0.5)]
TIER1_NAMES = {name_of(**j) for j in TIER1}


def tier2_jobs():
    bf = [dict(n=g, ibm=False) for g in GRIDS2 + tuple("h%d" % h for h in HPLUS)]
    jobs = [dict(n=g, ibm=True, f=f, law=law) for g in GRIDS2 for f in FS for law in ("spalartAllmaras", "none")]
    jobs += [dict(n=273, ibm=True, f=0.25, ma=0.0), dict(n=273, ibm=True, f=0.5, opt="noConv"),
             dict(n=273, ibm=True, f=0.5, opt="convLaw"), dict(n=273, ibm=True, f=0.5, opt="lead")]
    for h in HPLUS:
        jobs += [dict(n="h%d" % h, ibm=True, f=f) for f in FS]
        jobs += [dict(n="h%d" % h, ibm=True, f=0.5, opt=o) for o in ("noConv", "convLaw")]
    names = set()
    out = []
    for j in bf + jobs:                       # tier-1 runs are not repeated
        nm = name_of(**j)
        if nm not in TIER1_NAMES and nm not in names:
            names.add(nm)
            out.append(j)
    return out


def run(ctx, res):
    jobs = list(TIER1)
    if ctx.tier2:
        # longest first: the 545 grid
        jobs += sorted(tier2_jobs(), key=lambda j: j["n"] != 545)
    rows = common.pmap(lambda j: one(ctx, **j), jobs, ctx.jobs)
    res.data["rows"] = rows
    res.data["metrics"] = {r["name"] + "_C_D": r["C_D"] for r in rows if np.isfinite(r["C_D"])}
    if ctx.tier2:
        for r in rows:
            for k in ("cf097", "nut097"):
                if r.get(k) is not None:
                    res.data["metrics"]["%s_%s" % (r["name"], k)] = r[k]
    t1 = [r for r in rows if r["name"] in TIER1_NAMES]
    res.le("runs not converged (residualControl, 40000 iterations)", sum(not r["converged"] for r in t1), 0)
    bf = [r for r in t1 if not r["ibm"]][0]
    ib = [r for r in t1 if r["ibm"]][0]
    res.rel("plate C_D, IBM f 0.5 vs body-fitted, grid 273 x 193", ib["C_D"], bf["C_D"], 0.01)
    inc = tmr.flatplate_cf()["FUN3D incompressible SA"]
    res.report("body-fitted C_f(0.97) vs TMR incompressible FUN3D (finest grid)", bf["cf097"], float(inc[1][0]),
               "deviation %+.2f %%" % (100*(bf["cf097"]/inc[1][0] - 1)))
    cd = tmr.flatplate_drag()["CFL3D"]
    res.report("body-fitted C_D vs TMR CFL3D (compressible M 0.2, finest)", bf["C_D"], float(cd[1][0]),
               "deviation %+.2f %%" % (100*(bf["C_D"]/cd[1][0] - 1)))
    if ctx.tier2:
        tier2_report(res, rows)


def cfx_dev(r, b):
    """max |C_f,IBM/C_f,BF - 1| over the body-fitted face centres in XWIN"""
    if r.get("cfx") is None or b.get("cfx") is None:
        return float("nan")
    bx, bc = np.array(b["cfx"]), np.array(b["cf"])
    m = (bx >= XWIN[0]) & (bx <= XWIN[1])
    return float(np.max(np.abs(np.interp(bx[m], r["cfx"], r["cf"])/bc[m] - 1)))


def tier2_report(res, rows):
    """IBM against body-fitted on the same grid (doc/validationSuite.md §11 V15): wall-resolved
    C_f(0.97), C_D 1 %, C_f(x) 2 % on 0.1 <= x <= 1.9, max nu_t/nu at 0.97 2 %, spread over f 1 %;
    wall functions (user, 2026-11-22, as V7): C_f(0.97), C_D within 3 % of the wall-resolved
    body-fitted 273 run (same x lines) and closer to it than the body-fitted wall function on the
    same mesh (Spalding's law is 5 % off at h+ 90); off-default options there REPORT"""
    bf = {r["n"]: r for r in rows if not r["ibm"]}
    t2 = [r for r in rows if r["name"] not in TIER1_NAMES or r["ibm"]]
    bad = [r["name"] for r in rows if not r["converged"]]
    res.le("tier 2: runs not converged (residualControl) or failed", len(bad), 0, ", ".join(bad) or "-")
    inc = float(tmr.flatplate_cf()["FUN3D incompressible SA"][1][0])
    for g in GRIDS2:
        if g in bf and np.isfinite(bf[g].get("cf097", np.nan)):
            res.report("tier 2 body-fitted %s: C_f(0.97) vs TMR incompressible" % g, bf[g]["cf097"], inc,
                       "deviation %+.2f %%" % (100*(bf[g]["cf097"]/inc - 1)))
    for r in t2:
        if not r["ibm"] or r["n"] not in bf:
            continue
        b = bf[r["n"]]
        wf = isinstance(r["n"], str)
        tag = "tier 2 %s" % r["name"]
        if not np.isfinite(r["C_D"]) or r.get("cf097") is None:
            res.report(tag, float("nan"), None, r.get("error", "no result"))
            continue
        if wf:
            wr = bf[273]
            for k, lab in (("cf097", "C_f(0.97)"), ("C_D", "C_D")):
                e, ebf = r[k]/wr[k] - 1, b[k]/wr[k] - 1
                note = "IBM %+.2f %%, body-fitted wall function %+.2f %%" % (100*e, 100*ebf)
                if r.get("opt"):
                    res.report(tag + " %s vs wall-resolved body-fitted (option)" % lab, r[k], wr[k], note)
                else:
                    res.add(tag + " %s vs wall-resolved body-fitted" % lab,
                            "PASS" if abs(e) <= RTOL_WF and abs(e) < abs(ebf) else "FAIL", r[k], wr[k],
                            "%g %% and < body-fitted wall function" % (100*RTOL_WF), note)
            res.report(tag + " max C_f(x) deviation from wall-resolved body-fitted, 0.1..1.9", cfx_dev(r, wr), None,
                       "body-fitted wall function %.2f %%" % (100*cfx_dev(b, wr)))
        else:
            res.rel(tag + " C_f(0.97) vs body-fitted", r["cf097"], b["cf097"], RTOL_WR)
            res.rel(tag + " C_D vs body-fitted", r["C_D"], b["C_D"], RTOL_WR)
            res.le(tag + " max C_f(x) deviation, 0.1..1.9", cfx_dev(r, b), RTOL_CFX)
            res.rel(tag + " max nu_t/nu at x 0.97", r["nut097"], b["nut097"], RTOL_NUT)
    for g in GRIDS2:
        for law in ("spalartAllmaras", "none"):
            v = [r["cf097"] for r in rows if r["ibm"] and r["n"] == g and r["law"] == law and not r.get("opt")
                 and r["ma"] == 0.5 and r.get("cf097") is not None]
            if len(v) == len(FS) and g in bf:
                res.le("tier 2 spread of C_f(0.97) over f, grid %s, law %s" % (g, law),
                       (max(v) - min(v))/bf[g]["cf097"], RTOL_SPREAD)
