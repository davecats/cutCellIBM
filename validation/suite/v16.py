"""V16 TMR 2DBC: turbulent bump in a channel, Spalart-Allmaras, Re = 3e6 per unit length (NASA TMR;
reference_tmr.py; doc/validationSuite.md §11). The immersed wall needs a wall function: y+ <= 1 on a
curved wall is out of reach for a Cartesian mesh. Immersed meshes: tensor-product grids with a uniform
band dx = 2e-3, dy = 2, 4, 8e-4 (h+ ~ 33, 67, 133 at x = 0.75) for -0.1 <= x <= 1.6, y <= 0.13,
stretched (ratio 1.1) to the TMR domain [-25, 26.5] x [0, 5]; the wall 0 <= x <= 1.5 (bump
y = 0.05 sin^4(pi x/0.9 - pi/3) on 0.3..1.2) an STL body, its ends mid-cell, the flat parts cutting the
bottom row with fluid fraction f = 0.5; the mesh bottom a symmetry plane (1 - f) dy below. Body-fitted
references: the TMR 705 x 321 grid coarsened at the wall to the same first-cell height (Spalding's
law, as V8), and the TMR 705 x 321 grid itself, wall-resolved, against the TMR values (REPORT: the
OpenFOAM SA is not the TMR SA). Quantities: C_p(x), C_f(x), C_f at x = 0.632, 0.75, 0.868, C_L,
C_Dp, C_Dv, u and nu_t/nu profiles at x = 0.75, 1.20148. Tier 2, parallel runs of NP ranks."""
import gzip
import os
import re

import numpy as np

import common
import reference_tmr as tmr
import v15
from common import write_dict
from foamio import latest_time, read_field

ID, KIND, TIER = "V16", "case", 2
TITLE = "TMR bump-in-channel 2DBC, SA, immersed vs body-fitted wall functions"
STRESSES, RELIES = ["K13"], ["V15"]
REGRESSION = (1e-4,)

NU = 1/3e6
DZ = 0.01
NP = 8
END = 120000                    # SIMPLE cap; runs past 40000 continue from their latest time
END_MESH = {"h33": 240000}      # h+ 33 band: immersed C_L still converging geometrically at 120000 (user, 2026-11-24)
STAT = (0.1, 1e-4)              # stationary: C_L, C_D within 1e-4 relative over the last 10 % of the iterations
GRIDS = common.RUNS / "tmr_grids" / "bump"   # bumpgrids-grids.zip of the TMR, the p2dfmt files
GRID705 = "bump_1levelsdown_705x321.p2dfmt"
DX = 2e-3
X0 = -0.101                     # band lines X0 + k DX: x = 0 and 1.5 are mid-cell
NXB = 851                       # band cells: -0.101 .. 1.601
YBAND = 0.13                    # bump top 0.05 + 0.08
F = 0.5
MESHES = {"h33": 2e-4, "h67": 4e-4, "h133": 8e-4}
DEFAULT = "h67"
XCF = tmr.BUMP["x_cf"]
XPROF = (0.75, 1.20148)
LREF = tmr.BUMP["ref_length"]
CP_WIN = (0.1, 1.4)
TOL_CP, RTOL_CF, RTOL_CL = 0.005, 0.03, 0.02
OPTIONS = {"noConv": dict(ibmConvection=False), "convLaw": dict(ibmConvectionWallLaw=True)}


# ---------------------------------------------------------------------------------------- meshes

def stretch(d0, length, r=1.1):
    """offsets of the lines of a geometric stretch from spacing d0 r, scaled to end at length"""
    d = []
    while sum(d) < length:
        d.append(d0*r**(len(d) + 1))
    d = np.array(d)*length/sum(d)
    return np.cumsum(d)


def ibm_lines(dy, f=F):
    xb = X0 + DX*np.arange(NXB + 1)
    xs = np.concatenate([(xb[0] - stretch(DX, xb[0] + 25.0))[::-1], xb, xb[-1] + stretch(DX, 26.5 - xb[-1])])
    y0 = -(1 - f)*dy
    yb = y0 + dy*np.arange(int(np.ceil((YBAND - y0)/dy)) + 1)
    ys = np.concatenate([yb, yb[-1] + stretch(dy, 5.0 - yb[-1])])
    return xs, ys


def read_tmr(name):
    p = GRIDS / name
    tok = (p.read_text() if p.exists() else gzip.open(str(p) + ".gz", "rt").read()).split()
    ni, nj = int(tok[1]), int(tok[2])
    v = np.array(tok[3:3 + 2*ni*nj], float)
    return v[:ni*nj].reshape(nj, ni), v[ni*nj:].reshape(nj, ni)


def coarsened(h):
    """the 705 x 321 grid with the lines below the first one spaced >= h (column x = 0.75) replaced by
    k uniform cells along every i-line: first-cell height h within 1 % along the wall"""
    X, Y = read_tmr(GRID705)
    s = np.vstack([np.zeros((1, X.shape[1])), np.cumsum(np.hypot(np.diff(X, axis=0), np.diff(Y, axis=0)), axis=0)])
    i0 = int(np.argmin(np.abs(X[0] - 0.75)))
    j = int(np.argmax(np.diff(s[:, i0]) >= h))
    k = max(int(round(s[j, i0]/h)), 1)
    Xn, Yn = np.empty((k + X.shape[0] - j, X.shape[1])), np.empty((k + X.shape[0] - j, X.shape[1]))
    for i in range(X.shape[1]):
        sn = np.linspace(0.0, s[j, i], k + 1)
        Xn[:k + 1, i] = np.interp(sn, s[:j + 1, i], X[:j + 1, i])
        Yn[:k + 1, i] = np.interp(sn, s[:j + 1, i], Y[:j + 1, i])
    Xn[k + 1:], Yn[k + 1:] = X[j + 1:], Y[j + 1:]
    return Xn, Yn


def write_p2d_xy(path, X, Y):
    with open(path, "w") as fh:
        fh.write("1\n%d %d\n" % (X.shape[1], X.shape[0]))
        for a in (X.ravel(), Y.ravel()):
            for i in range(0, len(a), 4):
                fh.write(" ".join("%.15e" % v for v in a[i:i + 4]) + "\n")


def write_stl(path, x0=0.0, x1=1.5, ds=2.5e-4):
    """closed binary STL of the wall body: the wall y = bump_y(x) on x0..x1 over y = -1, z = -1..1,
    triangles oriented outward"""
    x = np.linspace(x0, x1, int(round((x1 - x0)/ds)) + 1)
    y = tmr.bump_y(x)
    tris = []

    def tri(a, b, c, out):
        a, b, c = np.array(a, float), np.array(b, float), np.array(c, float)
        n = np.cross(b - a, c - a)
        if np.dot(n, out) < 0:
            b, c, n = c, b, -n
        tris.append((n/np.linalg.norm(n), a, b, c))

    lo, hi, yb = -1.0, 1.0, -1.0
    for k in range(len(x) - 1):
        xa, xc, ya, yc = x[k], x[k + 1], y[k], y[k + 1]
        up = (-(yc - ya), xc - xa, 0.0)
        tri((xa, ya, lo), (xc, yc, lo), (xc, yc, hi), up)
        tri((xa, ya, lo), (xc, yc, hi), (xa, ya, hi), up)
        # bottom split like the wall: the edges of the caps match (closed, no T-junctions)
        tri((xa, yb, lo), (xc, yb, lo), (xc, yb, hi), (0, -1, 0))
        tri((xa, yb, lo), (xc, yb, hi), (xa, yb, hi), (0, -1, 0))
        for z, out in ((hi, (0, 0, 1)), (lo, (0, 0, -1))):
            tri((xa, yb, z), (xc, yb, z), (xc, yc, z), out)
            tri((xa, yb, z), (xc, yc, z), (xa, ya, z), out)
    for xe, ye, out in ((x[0], y[0], (-1, 0, 0)), (x[-1], y[-1], (1, 0, 0))):
        tri((xe, yb, lo), (xe, ye, lo), (xe, ye, hi), out)
        tri((xe, yb, lo), (xe, ye, hi), (xe, yb, hi), out)
    rec = np.zeros(len(tris), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    for i, (n, a, b, c) in enumerate(tris):
        rec["n"][i], rec["v"][i] = n, (a, b, c)
    with open(path, "wb") as fh:
        fh.write(b"V16 bump wall".ljust(80, b" "))
        fh.write(np.uint32(len(tris)).tobytes())
        fh.write(rec.tobytes())


# ------------------------------------------------------------------------------------------ cases

def name_of(mesh, ibm, opt=None):
    return ("ibm_%s%s" % (mesh, "_" + opt if opt else "")) if ibm else "bf_%s" % mesh


def case(d, mesh, ibm, opt=None):
    """mesh: "h33", "h67", "h133" (immersed or body-fitted wall function) or "wr" (body-fitted 705)"""
    (d / "system").mkdir(parents=True, exist_ok=True)
    if ibm:
        xs, ys = ibm_lines(MESHES[mesh])
        v15.write_p2d(d / "grid.p2dfmt", xs, ys)
    else:
        X, Y = read_tmr(GRID705) if mesh == "wr" else coarsened(MESHES[mesh])
        write_p2d_xy(d / "grid.p2dfmt", X, Y)
    import geometry
    geometry.control_dict(d, "simpleFoam")
    wall = '    { name bump;   type faceSet; action new; source boxToFace; box (-1e-6 -1 -1) (1.500001 0.06 1); }\n' \
           '    { name bump;   type faceSet; action subset; source faceToFace; set bottom; }\n' if not ibm else ""
    write_dict(d / "system" / "topoSetDict", "dictionary", "topoSetDict", """
actions
(
    { name bnd;    type faceSet; action new; source patchToFace; patch defaultFaces; }
    { name fb;     type faceSet; action new; source normalToFace; normal (0 0 1); cos 0.9; }
    { name fb;     type faceSet; action add; source normalToFace; normal (0 0 -1); cos 0.9; }
    { name fb;     type faceSet; action subset; source faceToFace; set bnd; }
    { name inlet;  type faceSet; action new; source boxToFace; box (-26 -2 -1) (-24.9999 6 1); }
    { name inlet;  type faceSet; action subset; source faceToFace; set bnd; }
    { name outlet; type faceSet; action new; source boxToFace; box (26.4999 -2 -1) (27 6 1); }
    { name outlet; type faceSet; action subset; source faceToFace; set bnd; }
    { name top;    type faceSet; action new; source boxToFace; box (-26 4.9999 -1) (27 6 1); }
    { name top;    type faceSet; action subset; source faceToFace; set bnd; }
    { name bottom; type faceSet; action new; source faceToFace; set bnd; }
    { name bottom; type faceSet; action subtract; source faceToFace; set fb; }
    { name bottom; type faceSet; action subtract; source faceToFace; set inlet; }
    { name bottom; type faceSet; action subtract; source faceToFace; set outlet; }
    { name bottom; type faceSet; action subtract; source faceToFace; set top; }
%s    { name sym;    type faceSet; action new; source faceToFace; set bottom; }
%s);
""" % (wall, "    { name sym;    type faceSet; action subtract; source faceToFace; set bump; }\n" if not ibm else ""))
    pl = '    { name bump; patchInfo { type wall; }; constructFrom set; set bump; }\n' if not ibm else ""
    write_dict(d / "system" / "createPatchDict", "dictionary", "createPatchDict", """
pointSync false;
patches
(
    { name inlet;  patchInfo { type patch; }; constructFrom set; set inlet; }
    { name outlet; patchInfo { type patch; }; constructFrom set; set outlet; }
    { name top;    patchInfo { type symmetryPlane; }; constructFrom set; set top; }
    { name sym;    patchInfo { type symmetryPlane; }; constructFrom set; set sym; }
%s    { name frontAndBack; patchInfo { type empty; }; constructFrom set; set fb; }
);
""" % pl)
    nut_inf = 3*NU*27/(27 + 7.1**3)
    bc = {"inlet": {"U": "fixedValue; value uniform (1 0 0);", "p": "zeroGradient;",
                    "nuTilda": "fixedValue; value uniform %r;" % (3*NU), "nut": "calculated; value uniform %r;" % nut_inf},
          "outlet": {"U": "zeroGradient;", "p": "fixedValue; value uniform 0;", "nuTilda": "zeroGradient;",
                     "nut": "calculated; value uniform %r;" % nut_inf},
          "top": {k: "symmetryPlane;" for k in ("U", "p", "nuTilda", "nut")},
          "sym": {k: "symmetryPlane;" for k in ("U", "p", "nuTilda", "nut")},
          "frontAndBack": {k: "empty;" for k in ("U", "p", "nuTilda", "nut")}}
    if not ibm:
        bc["bump"] = {"U": "noSlip;", "p": "zeroGradient;", "nuTilda": "fixedValue; value uniform 0;",
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
        (d / "constant" / "triSurface").mkdir(parents=True, exist_ok=True)
        write_stl(d / "constant" / "triSurface" / "bump.stl")
        geometry.ibm_properties(d, {"bump": {"type": "triSurfaceMesh", "file": '"bump.stl"', "wallType": "noSlip",
                                             "wallFunction": "spalartAllmaras"}},
                                mergeAlpha=0.5, nSubSamples=32, **(OPTIONS.get(opt) or {}))
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
    v15.fv_solution(d, ibm)
    fo = "" if ibm else """
functions
{
    wss { type wallShearStress; libs (fieldFunctionObjects); patches (bump); writeControl writeTime; }
    forces { type forces; libs (forces); patches (bump); rho rhoInf; rhoInf 1; CofR (0 0 0);
             writeControl timeStep; writeInterval 1; }
}
"""
    write_dict(d / "system" / "controlDict", "dictionary", "controlDict", """
application %s;
%sstartFrom latestTime; startTime 0; stopAt endTime; endTime %d; deltaT 1;
writeControl timeStep; writeInterval %d; purgeWrite 0; writeFormat ascii; writePrecision 12;
timeFormat general; timePrecision 6; runTimeModifiable false;
%s""" % ("ibmSimpleFoam" if ibm else "simpleFoam", 'libs ("libcutCellIBMTurbulence.so");\n' if ibm else "",
         END_MESH.get(mesh, END), END_MESH.get(mesh, END), fo))
    write_dict(d / "system" / "decomposeParDict", "dictionary", "decomposeParDict",
               "numberOfSubdomains %d;\nmethod scotch;\n%s" % (NP, common.LINK_CONSTRAINT if ibm else ""))
    common.sh("plot3dToFoam -2D %r -noBlank grid.p2dfmt" % DZ, d, "log.plot3dToFoam")
    common.sh("topoSet", d, "log.topoSet")
    common.sh("createPatch -overwrite", d, "log.createPatch")


def one(ctx, mesh, ibm, opt=None):
    name = name_of(mesh, ibm, opt)
    d = ctx.root / name
    try:
        return _one(d, name, mesh, ibm, opt)
    except common.RunError as e:          # a diverged run is a result, not a driver crash
        return {"name": name, "mesh": mesh, "ibm": ibm, "opt": opt, "converged": False, "error": str(e)[:300]}


def _one(d, name, mesh, ibm, opt):
    if not (common.REUSE and (d / "system").exists()):
        case(d, mesh, ibm, opt)
    (d / "0").mkdir(exist_ok=True)
    common.sh("checkMesh -constant", d, "log.checkMesh")
    common.sh("postProcess -func writeCellCentres -time 0", d, "log.postProcess")
    app = "ibmSimpleFoam" if ibm else "simpleFoam"
    if ibm:
        common.link_zone(d)
    common.sh("decomposePar -force", d, "log.decomposePar")
    common.sh(app, d, "log." + app, np_=NP)
    common.sh("reconstructPar -latestTime", d, "log.reconstructPar")
    if ibm:
        common.sh("ibmSetGeometry", d, "log.ibmSetGeometry")   # geometry fields at 0
    log = (d / ("log." + app)).read_text()
    t = latest_time(str(d))
    out = {"name": name, "mesh": mesh, "ibm": ibm, "opt": opt, "converged": "SIMPLE solution converged" in log,
           "iterations": int(t) if t.isdigit() else None}
    out.update(evaluate_ibm(d, t, log) if ibm else evaluate_bf(d, t))
    out["stat_dev"] = stationarity(force_history(d, ibm, log))
    out["stationary"] = out["stat_dev"] <= STAT[1]
    return out


def force_history(d, ibm, log):
    """C_L, C_D per iteration of the latest solver log (immersed) or of all forces/<t>/force.dat
    (body-fitted, written every step since the runs were continued past 40000)"""
    if ibm:
        pat = r"IBM force on bump: pressure \(([^)]*)\) viscous \(([^)]*)\)"
        F = np.array([[float(v) for v in (a + " " + b).split()] for a, b in re.findall(pat, log)])
        if not len(F):
            return np.zeros((0, 2))
        Fx, Fy = F[:, 0] + F[:, 3], F[:, 1] + F[:, 4]
    else:
        rows = []
        for f in sorted((d / "postProcessing" / "forces").glob("*/force.dat"), key=lambda q: float(q.parent.name)):
            for line in f.read_text().splitlines():
                if line.strip() and not line.startswith("#"):
                    v = [float(x) for x in line.replace("(", " ").replace(")", " ").split()]
                    rows.append((v[0], v[1], v[2]))
        if not rows:
            return np.zeros((0, 2))
        rows = np.array(sorted(dict((r[0], r) for r in rows).values()))
        Fx, Fy = rows[:, 1], rows[:, 2]
    s = 2/LREF/DZ
    return np.column_stack([s*Fy, s*Fx])


def stationarity(h):
    """max relative change of C_L, C_D over the last STAT[0] of the history from its last value
    (inf for fewer than 1000 entries)"""
    if len(h) < 1000:
        return float("inf")
    w = h[int((1 - STAT[0])*len(h)):]
    return float(np.max(np.abs(w - h[-1])/np.abs(h[-1])))


# ------------------------------------------------------------------------------------ evaluation

def coeffs(Fp, Fv):
    """C_L, C_Dp, C_Dv of the force per unit depth (U = 1, reference length 1.5)"""
    s = 2/LREF
    return {"C_L": s*float(Fp[1] + Fv[1]), "C_Dp": s*float(Fp[0]), "C_Dv": s*float(Fv[0]),
            "C_D": s*float(Fp[0] + Fv[0])}


def stations(x, cf):
    return {"cf_%g" % xs: float(np.interp(xs, x, cf)) for xs in XCF}


def evaluate_ibm(d, t, log):
    f = common.ibm_force(log, "bump")
    out = coeffs(f["pressure"]/DZ, f["viscous"]/DZ)
    C, _ = read_field(os.path.join(d, "0", "C"))
    Sw, _ = read_field(os.path.join(d, "0", "ibmSw"))
    a, _ = read_field(os.path.join(d, "0", "ibmAlpha"))
    p, _ = read_field(os.path.join(d, t, "p"))
    tau, _ = read_field(os.path.join(d, t, "ibmWallShear"))
    A = np.sqrt(np.sum(Sw**2, axis=1))
    w = A > 0
    # per mesh column: area-weighted wall pressure and wall shear of its wall cells
    col = np.round((C[w, 0] - X0)/DX - 0.5).astype(int)
    cols = np.unique(col)
    Aw, pw, tw, xw = A[w], p[w], tau[w], C[w, 0]
    xc = np.array([xw[col == c].mean() for c in cols])
    cp = np.array([np.sum(Aw[col == c]*2*pw[col == c])/np.sum(Aw[col == c]) for c in cols])
    cf = np.array([np.sum(Aw[col == c]*2*tw[col == c])/np.sum(Aw[col == c]) for c in cols])
    m = (xc > 0) & (xc < 1.5)
    out.update({"x": xc[m].tolist(), "cp": cp[m].tolist(), "cf": cf[m].tolist()})
    out.update(stations(xc[m], cf[m]))
    off = os.path.join(d, t, "ibmCentroidOffset")
    Cc = C + read_field(off)[0] if os.path.exists(off) else C
    U, _ = read_field(os.path.join(d, t, "U"))
    nut, _ = read_field(os.path.join(d, t, "nut"))
    fl = a > 0
    out.update(profiles_ibm(C[fl], Cc[fl], U[fl], nut[fl]))
    return out


def profiles_ibm(C, Cc, U, nut):
    """the mesh column (cell centres C) at each station; y from the wall at the fluid centroids Cc"""
    out = {}
    for xs in XPROF + (0.75,):
        m = np.abs(C[:, 0] - C[np.argmin(np.abs(C[:, 0] - xs)), 0]) < 1e-9
        o = np.argsort(Cc[m, 1])
        y = Cc[m, 1][o] - float(tmr.bump_y(xs))
        out["prof_%g" % xs] = {"y": y.tolist(), "u": U[m, 0][o].tolist(), "nut_nu": (nut[m][o]/NU).tolist()}
    out["nut075"] = float(np.max(out["prof_0.75"]["nut_nu"]))
    return out


def owner_of_patch(d, patch):
    pm = os.path.join(d, "constant", "polyMesh")
    m = re.search(patch + r"\s*\{[^}]*nFaces\s+(\d+);\s*startFace\s+(\d+);", open(os.path.join(pm, "boundary")).read())
    n, s = int(m.group(1)), int(m.group(2))
    txt = open(os.path.join(pm, "owner")).read()
    txt = txt[txt.index("(", txt.index("}")) + 1:]
    return np.array(txt[:txt.index(")")].split(), int)[s:s + n]


def evaluate_bf(d, t):
    fo = sorted((d / "postProcessing" / "forces").glob("*/force.dat"), key=lambda q: float(q.parent.name))
    vals = [float(v) for v in fo[-1].read_text().strip().splitlines()[-1].replace("(", " ").replace(")", " ").split()]
    out = coeffs(np.array(vals[4:7])/DZ, np.array(vals[7:10])/DZ)
    _, pc = read_field(os.path.join(d, "0", "C"))
    p, _ = read_field(os.path.join(d, t, "p"))
    _, pt = read_field(os.path.join(d, t, "wallShearStress"))
    x = pc["bump"][:, 0]
    o = np.argsort(x)
    # p is zeroGradient on the wall: the pressure of the wall cells
    x, cp, cf = x[o], 2*p[owner_of_patch(d, "bump")][o], 2*np.sqrt(np.sum(pt["bump"]**2, axis=1))[o]
    out.update({"x": x.tolist(), "cp": cp.tolist(), "cf": cf.tolist()})
    out.update(stations(x, cf))
    C, _ = read_field(os.path.join(d, "0", "C"))
    U, _ = read_field(os.path.join(d, t, "U"))
    nut, _ = read_field(os.path.join(d, t, "nut"))
    X, _ = read_tmr(GRID705)
    out.update(profiles_bf(C, U, nut, X.shape[1] - 1))   # plot3dToFoam: cell (i, j) = j ncol + i
    return out


def profiles_bf(C, U, nut, ncol):
    out = {}
    idx = np.arange(len(C))
    for xs in XPROF:
        i = int(np.argmin(np.abs(C[:ncol, 0] - xs)))
        m = idx % ncol == i
        o = np.argsort(C[m, 1])
        y = C[m, 1][o] - float(tmr.bump_y(C[i, 0]))
        out["prof_%g" % xs] = {"y": y.tolist(), "u": U[m, 0][o].tolist(), "nut_nu": (nut[m][o]/NU).tolist()}
    i = int(np.argmin(np.abs(C[:ncol, 0] - 0.75)))
    if not np.allclose(C[i::ncol, 0], C[i, 0], atol=0.01):
        raise ValueError("cell order of the body-fitted mesh is not (i, j) = j ncol + i")
    out["nut075"] = float(np.max(nut[idx % ncol == i])/NU)
    return out


# ---------------------------------------------------------------------------------------- driver

def jobs():
    j = [dict(mesh="wr", ibm=False)]
    j += [dict(mesh=m, ibm=False) for m in MESHES] + [dict(mesh=m, ibm=True) for m in MESHES]
    j += [dict(mesh=DEFAULT, ibm=True, opt=o) for o in OPTIONS]
    return j


def run(ctx, res):
    rows = common.pmap(lambda j: one(ctx, **j), jobs(), max(1, ctx.jobs//NP))
    res.data["rows"] = rows
    res.data["metrics"] = {"%s_%s" % (r["name"], k): r[k] for r in rows
                           for k in ("C_L", "C_Dp", "C_Dv") + tuple("cf_%g" % x for x in XCF) if r.get(k) is not None}
    report(res, rows)


def cp_dev(r, b):
    """max |C_p,IBM - C_p,BF| at the immersed columns in CP_WIN"""
    x = np.array(r["x"])
    m = (x >= CP_WIN[0]) & (x <= CP_WIN[1])
    return float(np.max(np.abs(np.array(r["cp"])[m] - np.interp(x[m], b["x"], b["cp"]))))


def report(res, rows):
    """immersed against the body-fitted wall function on the same first-cell height (§11 V16): C_p
    within 0.005 on 0.1 <= x <= 1.4, C_f at the three stations 3 %, C_L 2 %; C_Dp REPORT; against
    the body-fitted wall-resolved 705 run REPORT; off-default options REPORT"""
    by = {r["name"]: r for r in rows}
    bad = [r["name"] for r in rows if not (r["converged"] or r.get("stationary"))]
    res.le("runs neither converged (residualControl, %s iterations) nor stationary (C_L, C_D within %g over the "
           "last %d %%) or failed" % ("%d, h33 %d" % (END, END_MESH["h33"]), STAT[1], 100*STAT[0]), len(bad), 0,
           ", ".join("%s (%.1e)" % (r["name"], r.get("stat_dev", float("nan"))) for r in rows
                     if r["name"] in bad) or "-")
    wr = by.get("bf_wr")
    f = tmr.bump_forces()["CFL3D"]
    cfr = tmr.bump_cf()
    if wr and "C_L" in wr:
        for k in ("C_L", "C_Dp", "C_Dv"):
            ref = float(f[k][np.argmin(f["h"])])
            res.report("body-fitted 705 (wall-resolved) %s vs TMR CFL3D finest" % k, wr[k], ref,
                       "deviation %+.2f %%" % (100*(wr[k]/ref - 1)))
        for x in XCF:
            h, q = cfr[("CFL3D", x)]
            ref = float(q[np.argmin(h)])
            res.report("body-fitted 705 (wall-resolved) C_f(%g) vs TMR CFL3D finest" % x, wr["cf_%g" % x], ref,
                       "deviation %+.2f %%" % (100*(wr["cf_%g" % x]/ref - 1)))
    for r in rows:
        if not r["ibm"]:
            continue
        b = by.get(name_of(r["mesh"], False))
        tag = r["name"]
        if "C_L" not in r or not b or "C_L" not in b:
            res.report(tag, float("nan"), None, r.get("error", "no result"))
            continue
        vs_wr = ("; vs wall-resolved: IBM %+.2f %%, body-fitted wall function %+.2f %%"
                 if wr and "C_L" in wr else "")
        if r.get("opt"):
            res.report(tag + " C_L vs body-fitted wall function (option)", r["C_L"], b["C_L"],
                       "%+.2f %%" % (100*(r["C_L"]/b["C_L"] - 1)))
            for x in XCF:
                k = "cf_%g" % x
                res.report(tag + " C_f(%g) vs body-fitted wall function (option)" % x, r[k], b[k],
                           "%+.2f %%" % (100*(r[k]/b[k] - 1)))
            continue
        res.le(tag + " max |C_p - C_p,BF|, 0.1..1.4", cp_dev(r, b), TOL_CP)
        res.rel(tag + " C_L vs body-fitted wall function", r["C_L"], b["C_L"], RTOL_CL,
                (vs_wr % (100*(r["C_L"]/wr["C_L"] - 1), 100*(b["C_L"]/wr["C_L"] - 1))) if vs_wr else "")
        for x in XCF:
            k = "cf_%g" % x
            res.rel(tag + " C_f(%g) vs body-fitted wall function" % x, r[k], b[k], RTOL_CF,
                    (vs_wr % (100*(r[k]/wr[k] - 1), 100*(b[k]/wr[k] - 1))) if vs_wr else "")
        res.report(tag + " C_Dp vs body-fitted wall function", r["C_Dp"], b["C_Dp"],
                   "%+.2f %%" % (100*(r["C_Dp"]/b["C_Dp"] - 1)))
        res.report(tag + " C_Dv vs body-fitted wall function", r["C_Dv"], b["C_Dv"],
                   "%+.2f %%" % (100*(r["C_Dv"]/b["C_Dv"] - 1)))
