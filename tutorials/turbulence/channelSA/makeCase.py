#!/usr/bin/env python3
"""Spalart-Allmaras plane channel, body-fitted or with immersed walls.

    makeCase.py <caseDir> --kind bf|ibm --mesh wr|wf [--offset f] [--wallFunction spalartAllmaras|none]
                [--nh n --h1 h] [--mergeAlpha a] [--secondOrder] [--hPlus lo hi] [--wallSample s]
                [--incline p q]
                [--init bfCase]

Channel of half-height 1 between y = -1 and 1, periodic in x with one cell
(a one-dimensional problem), bulk velocity 1, nu = 1e-4 (Re_b = 2e4,
Re_tau about 550), driven by a mean-velocity force.

  --mesh wr   wall-resolved: 100 cells across, graded, first cell 1.8e-3 (y+ about 1)
  --mesh wf   wall functions: 12 uniform cells, first cell centre at y+ about 45

bf:  walls are patches (nuTilda 0, nut nutUSpaldingWallFunction), simpleFoam.
ibm: the same distribution of cells, extended by four wall-sized cells into
     the solid on either side and shifted upwards by f times the wall cell,
     with the walls as two immersed planes at y = -1 and 1: the bottom wall
     cuts its cell at the fluid fraction f, the top wall at 1 - f.
     ibmSimpleFoam with ibmSpalartAllmaras and wallDist ibmLevelSet.

--incline p q (ibm): walls inclined to the grid at tan(alpha) = p/q, the
     channel still of half-height 1 normal to the walls and driven along
     them. Uniform square cells of size h1 (--h1, or the first cell of
     --mesh), written directly as a polyMesh: q columns of a staircase band
     covering the channel and NS cells beyond each wall, periodic with the
     cyclic left/right shifted by (q h, p h) (the wall rises p cells over q
     columns), the rest of the band boundary a slip patch "outer" in the
     solid. --offset f shifts the lattice upwards by f h (alpha = 0: the
     same cut as without --incline). The q columns see q different cuts.
--init bfCase (with --incline): initial U and nuTilda from the converged
     profiles of a body-fitted run (its cell centres 0/C and latest time),
     interpolated in the wall-normal coordinate (U along the channel; a
     uniform fine mesh takes very long to develop them from rest).
"""
import argparse
import math
import os
import re
import numpy as np

NU = 1e-4
LX = 0.2
NS = 4          # solid cells beyond each immersed wall


def grading(n, h1):
    """Expansion ratio last/first of n geometric cells of first size h1 over length 1."""
    if abs(h1 * n - 1) < 1e-12:
        return 1.0, 1.0
    # total length h1 (r^n - 1)/(r - 1) increases with r; r < 1 when h1 > 1/n
    lo, hi = (1.0 + 1e-12, 2.0) if h1 * n < 1 else (0.05, 1.0 - 1e-12)
    for _ in range(200):
        r = 0.5 * (lo + hi)
        if h1 * (r**n - 1) / (r - 1) > 1:
            hi = r
        else:
            lo = r
    return r, r**(n - 1)


def header(cls, obj, loc=None):
    return ("FoamFile { version 2.0; format ascii; class %s; %sobject %s; }\n"
            % (cls, ('location "%s"; ' % loc) if loc else "", obj))


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "w").write(text)


def inclined_mesh(d, p, q, h, f):
    """polyMesh of the inclined channel (see --incline); returns the patch
    names and the inclination."""
    alpha = math.atan2(p, q)
    ca, sa = math.cos(alpha), math.sin(alpha)
    y0 = -1 + f * h                      # lattice lines y = y0 + k h

    # rows of column i: cells touching the band [lower wall - NS h, upper
    # wall + NS h] (normal distances), exact in integers so that column q
    # is column 0 shifted by p rows
    def rows(i):
        lo = (-1 - NS * h) / ca          # lower edge of the band at x = 0
        hi = (1 + NS * h) / ca
        # y of the band edges over the column: x from i h to (i + 1) h
        ylo = lo + i * h * sa / ca
        yhi = hi + (i + 1) * h * sa / ca
        jlo = math.floor((ylo - y0) / h)
        jhi = math.ceil((yhi - y0) / h)
        return jlo, jhi                  # cells jlo .. jhi - 1
    R = [rows(i) for i in range(q)]
    cells = {}
    for i in range(q):
        for j in range(*R[i]):
            cells[(i, j)] = len(cells)
    pts = {}

    def pt(i, k, z):
        key = (i, k, z)
        if key not in pts:
            pts[key] = len(pts)
        return pts[key]

    def xface(i, j):                     # face x = i h, normal +x
        return [pt(i, j, 0), pt(i, j + 1, 0), pt(i, j + 1, 1), pt(i, j, 1)]

    def yface(i, k):                     # face y = y0 + k h, normal +y
        return [pt(i, k, 0), pt(i, k, 1), pt(i + 1, k, 1), pt(i + 1, k, 0)]

    def zface(i, j, z):                  # normal +z
        return [pt(i, j, z), pt(i + 1, j, z), pt(i + 1, j + 1, z), pt(i, j + 1, z)]

    internal, outer, left, right, empty = [], [], [], [], []
    for (i, j), c in cells.items():
        # +y neighbour
        n = cells.get((i, j + 1))
        if n is not None:
            internal.append((c, n, yface(i, j + 1)))
        else:
            outer.append((c, yface(i, j + 1)))
        if (i, j - 1) not in cells:
            outer.append((c, yface(i, j)[::-1]))
        # +x neighbour
        if i + 1 < q:
            n = cells.get((i + 1, j))
            if n is not None:
                internal.append((c, n, xface(i + 1, j)))
            else:
                outer.append((c, xface(i + 1, j)))
        else:
            if (0, j - p) in cells:
                right.append((j, c, xface(q, j)))
            else:
                outer.append((c, xface(q, j)))
        if i > 0:
            if (i - 1, j) not in cells:
                outer.append((c, xface(i, j)[::-1]))
        else:
            if (q - 1, j + p) in cells:
                # vertex 0 opposite vertex 0 of the right face, the rest reversed
                left.append((j + p, c, [pt(0, j, 0), pt(0, j, 1), pt(0, j + 1, 1), pt(0, j + 1, 0)]))
            else:
                outer.append((c, xface(0, j)[::-1]))
        empty.append((c, zface(i, j, 1)))
        empty.append((c, zface(i, j, 0)[::-1]))
    internal.sort(key=lambda t: (t[0], t[1]))
    left.sort(); right.sort()
    assert [t[0] for t in left] == [t[0] for t in right], "cyclic faces do not match"

    faces, owner, neigh = [], [], []
    for o, n, fc in internal:
        faces.append(fc); owner.append(o); neigh.append(n)
    patches = []
    for name, lst, extra in (("left", [(c, fc) for _, c, fc in left], "type cyclic; neighbourPatch right; inGroups 1(cyclic);"),
                             ("right", [(c, fc) for _, c, fc in right], "type cyclic; neighbourPatch left; inGroups 1(cyclic);"),
                             ("outer", outer, "type patch;"),
                             ("frontAndBack", empty, "type empty; inGroups 1(empty);")):
        patches.append((name, len(lst), len(faces), extra))
        for c, fc in lst:
            faces.append(fc); owner.append(c)
    P = [None] * len(pts)
    for (i, k, z), n in pts.items():
        P[n] = (i * h, y0 + k * h, z * h)
    pm = os.path.join(d, "constant/polyMesh")
    note = 'note "nPoints:%d nCells:%d nFaces:%d nInternalFaces:%d";' % (len(P), len(cells), len(faces), len(internal))
    write(os.path.join(pm, "points"), header("vectorField", "points", "constant/polyMesh")
          + "%d\n(\n%s\n)\n" % (len(P), "\n".join("(%.15g %.15g %.15g)" % x for x in P)))
    write(os.path.join(pm, "faces"), header("faceList", "faces", "constant/polyMesh")
          + "%d\n(\n%s\n)\n" % (len(faces), "\n".join("4(%d %d %d %d)" % tuple(fc) for fc in faces)))
    write(os.path.join(pm, "owner"), header("labelList", "owner", "constant/polyMesh").replace("object", note + " object")
          + "%d\n(\n%s\n)\n" % (len(owner), "\n".join(map(str, owner))))
    write(os.path.join(pm, "neighbour"), header("labelList", "neighbour", "constant/polyMesh").replace("object", note + " object")
          + "%d\n(\n%s\n)\n" % (len(neigh), "\n".join(map(str, neigh))))
    write(os.path.join(pm, "boundary"), header("polyBoundaryMesh", "boundary", "constant/polyMesh")
          + "%d\n(\n%s\n)\n" % (len(patches), "\n".join(
              "    %s { %s nFaces %d; startFace %d; }" % (nm, ex, n, st) for nm, n, st, ex in patches)))
    xc = np.zeros((len(cells), 3))
    for (i, j), c in cells.items():
        xc[c] = ((i + 0.5) * h, y0 + (j + 0.5) * h, 0.5 * h)
    return len(cells), alpha, xc


def read_internal(path):
    """Internal field of an ascii OpenFOAM field file (nonuniform list)."""
    s = open(path).read()
    i = s.index("internalField")
    m = re.search(r"List<(\w+)>\s*(\d+)", s[i:])
    n, nc = int(m.group(2)), 3 if m.group(1) == "vector" else 1
    body = s[i + m.end():s.index("boundaryField", i)]
    v = np.array(body.replace("(", " ").replace(")", " ").split()[:n * nc], float)
    return v.reshape(n, nc)


def latest(case):
    ts = [t for t in os.listdir(case) if re.fullmatch(r"[0-9.eE+-]+", t) and t != "0"]
    return max(ts, key=float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--kind", choices=("bf", "ibm"), required=True)
    ap.add_argument("--mesh", choices=("wr", "wf"), required=True)
    ap.add_argument("--offset", type=float, default=0.5)
    ap.add_argument("--wallFunction", choices=("spalartAllmaras", "none"), default="spalartAllmaras")
    ap.add_argument("--nh", type=int, help="cells per half-channel (overrides --mesh)")
    ap.add_argument("--h1", type=float, help="first cell size (overrides --mesh)")
    ap.add_argument("--mergeAlpha", type=float, default=0.0,
                    help="ibm: link cut cells with alpha below this to a neighbour")
    ap.add_argument("--secondOrder", action="store_true",
                    help="ibm: second-order wall treatment (cutCellCorrected schemes)")
    ap.add_argument("--hPlus", type=float, nargs=2, metavar=("RESOLVED", "MODELLED"),
                    help="ibm: h+ range over which the wall law takes over (default 1 10)")
    ap.add_argument("--wallSample", type=float, default=1.5,
                    help="ibm: distance of the wall-law sample in wall-cell sizes (0: the wall cell)")
    ap.add_argument("--incline", type=int, nargs=2, metavar=("P", "Q"),
                    help="ibm: walls inclined at tan(alpha) = P/Q to the grid, uniform cells")
    ap.add_argument("--init", help="with --incline: initial profiles from this body-fitted run")
    a = ap.parse_args()
    d = a.case
    if a.incline and a.kind != "ibm":
        raise SystemExit("--incline needs --kind ibm")

    # cells across the channel: half-channel distribution, mirrored
    if a.mesh == "wr":
        nh, h1 = 50, 1.8e-3
    else:
        nh, h1 = 6, 1.0 / 6
    nh = a.nh or nh
    h1 = a.h1 or h1
    r, G = grading(nh, h1)
    grade = "((0.5 0.5 %.10g) (0.5 0.5 %.10g))" % (G, 1 / G)

    shift = a.offset * h1 if a.kind == "ibm" else 0.0
    if a.kind == "bf":
        ys = [(-1.0, 1.0, 2 * nh, grade)]
    else:
        s = NS * h1
        ys = [(-1.0 - s, -1.0, NS, "1"), (-1.0, 1.0, 2 * nh, grade), (1.0, 1.0 + s, NS, "1")]

    # blocks stacked in y, one cell in x and z
    verts, blocks = [], []
    yb = [ys[0][0]] + [b[1] for b in ys]
    for y in yb:
        y += shift
        verts += ["(0 %.12g 0)" % y, "(%g %.12g 0)" % (LX, y),
                  "(%g %.12g %g)" % (LX, y, LX), "(0 %.12g %g)" % (y, LX)]
    for k, (_, _, n, g) in enumerate(ys):
        b0, b1 = 4 * k, 4 * (k + 1)
        # hex ordered x, y, z: (x0y0z0 x1y0z0 x1y1z0 x0y1z0  x0y0z1 x1y0z1 x1y1z1 x0y1z1)
        blocks.append("hex (%d %d %d %d %d %d %d %d) (1 %d 1) simpleGrading (1 %s 1)"
                      % (b0, b0 + 1, b1 + 1, b1, b0 + 3, b0 + 2, b1 + 2, b1 + 3, n, g))
    nb = len(ys)
    left = " ".join("(%d %d %d %d)" % (4*k, 4*k + 3, 4*(k+1) + 3, 4*(k+1)) for k in range(nb))
    right = " ".join("(%d %d %d %d)" % (4*k + 1, 4*(k+1) + 1, 4*(k+1) + 2, 4*k + 2) for k in range(nb))
    front = " ".join("(%d %d %d %d)" % (4*k, 4*(k+1), 4*(k+1) + 1, 4*k + 1) for k in range(nb))
    back = " ".join("(%d %d %d %d)" % (4*k + 3, 4*k + 2, 4*(k+1) + 2, 4*(k+1) + 3) for k in range(nb))
    top = 4 * nb
    wtype = "wall" if a.kind == "bf" else "patch"
    alpha = 0.0
    if a.incline:
        ncells, alpha, xc = inclined_mesh(d, a.incline[0], a.incline[1], h1, a.offset)
    else:
      write(os.path.join(d, "system/blockMeshDict"), header("dictionary", "blockMeshDict") + """
// %s, %s mesh: first cell %.4g, geometric ratio %.6f%s
scale 1;
vertices
(
    %s
);
blocks
(
    %s
);
boundary
(
    left   { type cyclic; neighbourPatch right; faces (%s); }
    right  { type cyclic; neighbourPatch left;  faces (%s); }
    bottom { type %s; faces ((0 1 2 3)); }
    top    { type %s; faces ((%d %d %d %d)); }
    frontAndBack { type empty; faces (%s %s); }
);
""" % (a.kind, a.mesh, h1, r, (", shifted by %.4g (offset %.3g of the wall cell)" % (shift, a.offset)) if a.kind == "ibm" else "",
       "\n    ".join(verts), "\n    ".join(blocks), left, right, wtype, wtype,
       top, top + 3, top + 2, top + 1, front, back))

    # boundary conditions
    side = "    left { type cyclic; }\n    right { type cyclic; }\n    frontAndBack { type empty; }\n"
    walls = "outer" if a.incline else "(bottom|top)"
    ca, sa = math.cos(alpha), math.sin(alpha)
    Ubar = "(%.15g %.15g 0)" % (ca, sa)
    if a.kind == "bf":
        bU, bp = "noSlip;", "zeroGradient;"
        bNuT, bNut = "fixedValue; value uniform 0;", "nutUSpaldingWallFunction; value uniform 0;"
    else:
        bU, bp = "slip;", "zeroGradient;"
        bNuT, bNut = "zeroGradient;", "calculated; value uniform 0;"
    init = {}
    if a.init and a.incline:
        t = latest(a.init)
        ry = read_internal(os.path.join(a.init, "0/C"))[:, 1]
        o = np.argsort(ry)
        ru = read_internal(os.path.join(a.init, t, "U"))[o, 0]
        rn = read_internal(os.path.join(a.init, t, "nuTilda"))[o, 0]
        eta = xc @ np.array([-sa, ca, 0])
        inside = np.abs(eta) < 1
        u = np.where(inside, np.interp(eta, ry[o], ru), 0)
        nt = np.where(inside, np.interp(eta, ry[o], rn), 0)
        init["U"] = "nonuniform List<vector> %d\n(\n%s\n)" % (len(u), "\n".join(
            "(%.10g %.10g 0)" % (ui * ca, ui * sa) for ui in u))
        init["nuTilda"] = "nonuniform List<scalar> %d\n(\n%s\n)" % (len(nt), "\n".join("%.10g" % x for x in nt))
    for name, cls, dim, val, bc in (("U", "volVectorField", "[0 1 -1 0 0 0 0]", Ubar, bU),
                                    ("p", "volScalarField", "[0 2 -2 0 0 0 0]", "0", bp),
                                    ("nuTilda", "volScalarField", "[0 2 -1 0 0 0 0]", "%g" % (3 * NU), bNuT),
                                    ("nut", "volScalarField", "[0 2 -1 0 0 0 0]", "0", bNut)):
        write(os.path.join(d, "0", name), header(cls, name, "0") + """dimensions %s;
internalField %s;
boundaryField
{
%s    "%s" { type %s }
}
""" % (dim, init.get(name, "uniform " + val), side, walls, bc))

    write(os.path.join(d, "constant/transportProperties"), header("dictionary", "transportProperties")
          + "transportModel Newtonian;\nnu %g;\n" % NU)
    write(os.path.join(d, "constant/turbulenceProperties"), header("dictionary", "turbulenceProperties")
          + "simulationType RAS;\nRAS { RASModel %s; turbulence on; printCoeffs on; }\n"
          % ("SpalartAllmaras" if a.kind == "bf" else "ibmSpalartAllmaras"))
    write(os.path.join(d, "constant/fvOptions"), header("dictionary", "fvOptions") + """
momentumSource
{
    type            %s;
    selectionMode   all;
    fields          (U);
    Ubar            %s;
    relaxation      1.0;
}
""" % ("meanVelocityForce" if a.kind == "bf" else "ibmMeanVelocityForce", Ubar))
    if a.kind == "ibm":
        write(os.path.join(d, "constant/ibmProperties"), header("dictionary", "ibmProperties") + """
alphaMin      0.02;
mergeAlpha    %g;
secondOrder   %s;
%snSubSamples   32;
nBisect       40;
writeFields   true;
writeForces   true;
bodies
{
    lowerWall { type plane; point (%s); normal (%s); wallType noSlip; wallFunction %s; }
    upperWall { type plane; point (%s); normal (%s); wallType noSlip; wallFunction %s; }
}
""" % (a.mergeAlpha, "true" if a.secondOrder else "false",
       (("hPlusResolved %g;\nhPlusModelled %g;\n" % tuple(a.hPlus)) if a.hPlus else "")
       + ("wallSampleDistance %g;\n" % a.wallSample if a.wallSample != 1.5 else ""),
       "%.15g %.15g 0" % (sa, -ca), "%.15g %.15g 0" % (-sa, ca), a.wallFunction,
       "%.15g %.15g 0" % (-sa, ca), "%.15g %.15g 0" % (sa, -ca), a.wallFunction))

    write(os.path.join(d, "system/controlDict"), header("dictionary", "controlDict") + """
application     %s;
%sstartFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         20000;
deltaT          1;
writeControl    timeStep;
writeInterval   20000;
purgeWrite      0;
writeFormat     ascii;
writePrecision  12;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
""" % ("simpleFoam" if a.kind == "bf" else "ibmSimpleFoam",
       'libs ("libcutCellIBMTurbulence.so");\n' if a.kind == "ibm" else ""))
    write(os.path.join(d, "system/fvSchemes"), header("dictionary", "fvSchemes") + """
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
""" % ((("cutCellCorrected",) * 2 if a.secondOrder else ("corrected",) * 2)
       + ("meshWave" if a.kind == "bf" else "ibmLevelSet",)))
    write(os.path.join(d, "system/fvSolution"), header("dictionary", "fvSolution") + """
solvers
{
    p { solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0.01; }
    "(U|nuTilda)" { %s tolerance 1e-14; relTol 0.1; }
}
SIMPLE
{
    nNonOrthogonalCorrectors 0;
    consistent      yes;
    pRefPoint       %s;
    pRefValue       0;
    residualControl { p 1e-10; U 1e-10; nuTilda 1e-10; }
}
relaxationFactors
{
    equations { U 0.9; nuTilda 0.9; }
    fields    { p 1.0; }
}
""" % ("solver GAMG; smoother symGaussSeidel;" if a.incline else "solver smoothSolver; smoother symGaussSeidel;",
       "(%.15g %.15g %.15g)" % (0.5 * h1, 0.5 * h1 * sa / ca + 0.01, 0.5 * h1) if a.incline else "(0.1 0.01 0.1)"))
    print("%s: %s %s, %d cells across the channel, first cell %.4g (ratio %.5f), offset %s"
          % (d, a.kind, a.mesh, 2 * nh, h1, r, a.offset if a.kind == "ibm" else "-"))


if __name__ == "__main__":
    main()
