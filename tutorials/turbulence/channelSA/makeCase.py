#!/usr/bin/env python3
"""Spalart-Allmaras plane channel, body-fitted or with immersed walls.

    makeCase.py <caseDir> --kind bf|ibm --mesh wr|wf [--offset f] [--wallFunction spalding|none]

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
"""
import argparse
import os
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--kind", choices=("bf", "ibm"), required=True)
    ap.add_argument("--mesh", choices=("wr", "wf"), required=True)
    ap.add_argument("--offset", type=float, default=0.5)
    ap.add_argument("--wallFunction", choices=("spalding", "none"), default="spalding")
    ap.add_argument("--nh", type=int, help="cells per half-channel (overrides --mesh)")
    ap.add_argument("--h1", type=float, help="first cell size (overrides --mesh)")
    ap.add_argument("--mergeAlpha", type=float, default=0.0,
                    help="ibm: link cut cells with alpha below this to a neighbour")
    ap.add_argument("--secondOrder", action="store_true",
                    help="ibm: second-order wall treatment (cutCellCorrected schemes)")
    a = ap.parse_args()
    d = a.case

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
    if a.kind == "bf":
        bU, bp = "noSlip;", "zeroGradient;"
        bNuT, bNut = "fixedValue; value uniform 0;", "nutUSpaldingWallFunction; value uniform 0;"
    else:
        bU, bp = "slip;", "zeroGradient;"
        bNuT, bNut = "zeroGradient;", "calculated; value uniform 0;"
    for name, cls, dim, val, bc in (("U", "volVectorField", "[0 1 -1 0 0 0 0]", "(1 0 0)", bU),
                                    ("p", "volScalarField", "[0 2 -2 0 0 0 0]", "0", bp),
                                    ("nuTilda", "volScalarField", "[0 2 -1 0 0 0 0]", "%g" % (3 * NU), bNuT),
                                    ("nut", "volScalarField", "[0 2 -1 0 0 0 0]", "0", bNut)):
        write(os.path.join(d, "0", name), header(cls, name, "0") + """dimensions %s;
internalField uniform %s;
boundaryField
{
%s    "(bottom|top)" { type %s }
}
""" % (dim, val, side, bc))

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
    Ubar            (1 0 0);
    relaxation      1.0;
}
""" % ("meanVelocityForce" if a.kind == "bf" else "ibmMeanVelocityForce"))
    if a.kind == "ibm":
        write(os.path.join(d, "constant/ibmProperties"), header("dictionary", "ibmProperties") + """
alphaMin      0.02;
mergeAlpha    %g;
secondOrder   %s;
nSubSamples   32;
nBisect       40;
writeFields   true;
writeForces   true;
bodies
{
    lowerWall { type plane; point (0 -1 0); normal (0 1 0); wallType noSlip; wallFunction %s; }
    upperWall { type plane; point (0 1 0); normal (0 -1 0); wallType noSlip; wallFunction %s; }
}
""" % (a.mergeAlpha, "true" if a.secondOrder else "false", a.wallFunction, a.wallFunction))

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
    "(U|nuTilda)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-14; relTol 0.1; }
}
SIMPLE
{
    nNonOrthogonalCorrectors 0;
    consistent      yes;
    pRefPoint       (0.1 0.01 0.1);
    pRefValue       0;
    residualControl { p 1e-10; U 1e-10; nuTilda 1e-10; }
}
relaxationFactors
{
    equations { U 0.9; nuTilda 0.9; }
    fields    { p 1.0; }
}
""")
    print("%s: %s %s, %d cells across the channel, first cell %.4g (ratio %.5f), offset %s"
          % (d, a.kind, a.mesh, 2 * nh, h1, r, a.offset if a.kind == "ibm" else "-"))


if __name__ == "__main__":
    main()
