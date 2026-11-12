"""V15 TMR 2DZP: turbulent flat plate, Spalart-Allmaras, Re = 5e6 per unit length (NASA TMR;
reference_tmr.py; doc/validationSuite.md §11). Body-fitted (simpleFoam + SpalartAllmaras on the
TMR grid) against immersed (the same grid lines with one cell row added below the plate and the
grid shifted by f dy1: the plate, a box body, cuts that row with fluid fraction f). The body-fitted
run is the reference; the TMR values check it (REPORT: OpenFOAM's SA is not the TMR SA).
Tier 1: grid 273 x 193, wall-resolved, f = 0.5."""
import os
import subprocess

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
GRID = {273: "flatplate_clust2_1leveldown_273x193.p2dfmt", 545: "flatplate_clust2.p2dfmt"}
DZ = 0.01
XCF = tmr.FLATPLATE["x_cf"]


def read_grid(n):
    p = GRIDS / GRID[n]
    if not p.exists():
        subprocess.run(["gunzip", "-k", str(p) + ".gz"], check=True)
    tok = p.read_text().split()
    ni, nj = int(tok[1]), int(tok[2])
    v = np.array(tok[3:3 + 2*ni*nj], float)
    x, y = v[:ni*nj].reshape(nj, ni), v[ni*nj:].reshape(nj, ni)
    return x[0], y[:, 0]                      # tensor product: x_i, y_j


def write_p2d(path, xs, ys):
    X, Y = np.meshgrid(xs, ys)
    with open(path, "w") as fh:
        fh.write("1\n%d %d\n" % (len(xs), len(ys)))
        for a in (X.ravel(), Y.ravel()):
            for i in range(0, len(a), 4):
                fh.write(" ".join("%.15e" % v for v in a[i:i + 4]) + "\n")


def mesh(d, n, ibm, f):
    xs, ys = read_grid(n)
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


def case(d, n, ibm, f, ma=0.5, law="spalartAllmaras"):
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
        geometry.ibm_properties(d, {"plate": {"type": "box", "min": (0.5*dxle, -1.0, -1.0), "max": (3.0, 0.0, 1.0),
                                              "wallType": "noSlip", "wallFunction": law}}, mergeAlpha=ma,
                                nSubSamples=32)
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
%sstartFrom startTime; startTime 0; stopAt endTime; endTime 40000; deltaT 1;
writeControl timeStep; writeInterval 40000; purgeWrite 0; writeFormat ascii; writePrecision 12;
timeFormat general; timePrecision 6; runTimeModifiable false;
""" % ("ibmSimpleFoam" if ibm else "simpleFoam", 'libs ("libcutCellIBMTurbulence.so");\n' if ibm else ""))
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
    split over processors (README section 4). Body-fitted: the standard solvers."""
    solvers = ('    p { solver PCG; preconditioner { preconditioner linkedGroup; base { preconditioner GAMG; '
               'smoother DICGaussSeidel; nCellsInCoarsestLevel 10; tolerance 0; relTol 0; } } '
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


def one(ctx, n, ibm, f=0.5, ma=0.5):
    name = "%s_%d%s" % ("ibm" if ibm else "bf", n, "_f%g_ma%g" % (f, ma) if ibm else "")
    d = ctx.root / name
    if not (common.REUSE and (d / "system").exists()):
        case(d, n, ibm, f, ma)
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
        # C_D = F_x per unit depth/(rho U^2 L/2), L = 2; no wall shear field yet for C_f(x)
        CD = common.ibm_force(log, "plate")["total"][0]/DZ
        x, cf = None, None
    else:
        x, cf = cf_bodyfitted(d, t)
        o = np.argsort(x)
        x, cf = x[o], cf[o]
        dx = np.gradient(x)
        CD = float(np.sum(cf*dx)/2.0)
    out = {"name": name, "n": n, "ibm": ibm, "f": f, "ma": ma, "C_D": float(CD),
           "converged": "SIMPLE solution converged" in log}
    if x is not None:
        out["cf097"] = float(np.interp(XCF, x, cf))
    return out


def run(ctx, res):
    jobs = [(273, False), (273, True, 0.5)]
    rows = common.pmap(lambda j: one(ctx, *j), jobs, ctx.jobs)
    res.data["rows"] = rows
    res.data["metrics"] = {r["name"] + "_C_D": r["C_D"] for r in rows}
    res.le("runs not converged (residualControl, 40000 iterations)", sum(not r["converged"] for r in rows), 0)
    bf = [r for r in rows if not r["ibm"]][0]
    ib = [r for r in rows if r["ibm"]][0]
    res.rel("plate C_D, IBM f 0.5 vs body-fitted, grid 273 x 193", ib["C_D"], bf["C_D"], 0.01)
    inc = tmr.flatplate_cf()["FUN3D incompressible SA"]
    res.report("body-fitted C_f(0.97) vs TMR incompressible FUN3D (finest grid)", bf["cf097"], float(inc[1][0]),
               "deviation %+.2f %%" % (100*(bf["cf097"]/inc[1][0] - 1)))
    cd = tmr.flatplate_drag()["CFL3D"]
    res.report("body-fitted C_D vs TMR CFL3D (compressible M 0.2, finest)", bf["C_D"], float(cd[1][0]),
               "deviation %+.2f %%" % (100*(bf["C_D"]/cd[1][0] - 1)))
