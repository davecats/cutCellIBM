#!/usr/bin/env python3
"""NACA 0012 at Re = 1e6, Spalart-Allmaras, body-fitted or immersed.

    makeCase.py <caseDir> --kind bf|ibm [--alpha deg] [--mesh wf|wr] [--h h]
                [--mergeAlpha a] [--secondOrder] [--Nn n] [--grid g.p2dfmt] [--R r]
                [--endTime n] [--levels n --nFine n --nBand n] [--np n]
                [--wallFunction spalding|spalartAllmaras|none] [--wallSample s]
                [--wallLawDistance centroid|groupCentroid|cellAverage|halfCell] [--noNuTildaWallLaw]

Chord 1, U = 1, nu = 1e-6; the angle of attack is set by the inflow direction,
so both meshes keep the airfoil at zero incidence. NACA 0012 with the closed
trailing edge of the NASA Turbulence Modeling Resource,
    y = 0.594689181 (0.298222773 sqrt(x) - 0.127125232 x - 0.357907906 x^2
                     + 0.291984971 x^3 - 0.105174606 x^4).
Far field --R chords away (default 50; 500 to match the TMR grids),
freestream conditions, nuTilda_inf = 3 nu.

bf:  6-block C-mesh (blockMesh): surface and wake blocks above and below the
     airfoil and its wake; first cell 2e-3 (--mesh wf, y+ about 50, Spalding
     wall function) or 2.5e-5 (--mesh wr, y+ about 1). simpleFoam.
     --grid: a NASA TMR C-grid instead (e.g. n0012_449-129.p2dfmt, far field
     about 500 chords; https://tmbwg.github.io/turbmodels/naca0012_grids.html),
     converted with tmrGrid.py and plot3dToFoam; --mesh wf removes the grid
     lines next to the wall up to a first cell of 2e-3 at mid-chord. The
     reference runs use these grids (the blockMesh C-mesh with wall functions
     separates at alpha 10, as does the TMR grid with wall functions).
ibm: Cartesian mesh with a uniform core of cell size h (default 2e-3,
     y+ about 50 at the first cell centroid) over x in [-0.05, 1.05],
     y in [-0.1, 0.1], graded geometrically to the far field; the airfoil as
     an extruded STL (constant/triSurface/naca0012.stl). ibmSimpleFoam with
     ibmSpalartAllmaras, wallDist ibmLevelSet, cell linking mergeAlpha.
     --levels n: n levels of refinement (refineMesh, in the plane) in
     nested bands around the airfoil, down to h/2^n at the wall: the finest
     level within nFine of its cells from the wall (default 10), every
     coarser one nBand of its cells further (default 6), each also one cell
     of the next coarser level into the solid (refineSet.py). The wall law
     (--wallFunction, default spalding; --wallSample, --wallLawDistance,
     --noNuTildaWallLaw) as in channelSA/makeCase.py. --np n: run on n
     processors (scotch).
Forces: postProcessing/forces (bf) or postProcessing/ibmForces (ibm);
validation/naca0012_compare.py turns both into lift and drag coefficients.
"""
import argparse
import math
import os
import shutil
import subprocess
import sys
import numpy as np

NU = 1e-6
NUT_INF = 3 * NU
DZ = 0.01
R = 50.0          # far field distance


def yt(x):
    x = np.asarray(x, float)
    return 0.594689181 * (0.298222773 * np.sqrt(np.maximum(x, 0)) - 0.127125232 * x
                          - 0.357907906 * x**2 + 0.291984971 * x**3 - 0.105174606 * x**4)


def header(cls, obj, loc=None):
    return ("FoamFile { version 2.0; format ascii; class %s; %sobject %s; }\n"
            % (cls, ('location "%s"; ' % loc) if loc else "", obj))


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "w").write(text)


def ratio_first(L, n, first):
    """Geometric ratio r (cell i+1 / cell i) of n cells over L with the given first cell."""
    if abs(first * n - L) < 1e-14 * L:
        return 1.0
    lo, hi = (1 + 1e-12, 3.0) if first * n < L else (0.05, 1 - 1e-12)
    for _ in range(300):
        r = 0.5 * (lo + hi)
        if first * (r**n - 1) / (r - 1) > L:
            hi = r
        else:
            lo = r
    return r


def n_for(L, first, r):
    """Number of geometric cells with ratio r and first cell 'first' to fill about L."""
    return max(1, int(round(math.log(1 + L * (r - 1) / first) / math.log(r))))


def bf_mesh(d, first, Nn):
    xs = 0.3                                   # split point of the surface blocks
    ys = float(yt(xs))
    # surface distributions: LE block clustered at the LE, aft block towards the TE
    N1, N2, N3 = 60, 80, 80
    dLE, dTE = 1.5e-3, 2.5e-3
    s1 = np.linspace(0, 1, 2001); x1 = xs * (1 - np.cos(0.5 * np.pi * s1))   # LE..xs
    L1 = np.sum(np.hypot(np.diff(x1), np.diff(yt(x1))))
    r1 = ratio_first(L1, N1, dLE); last1 = dLE * r1**(N1 - 1)
    L2 = np.sum(np.hypot(np.diff(np.linspace(xs, 1, 2001)), np.diff(yt(np.linspace(xs, 1, 2001)))))
    r2 = ratio_first(L2, N2, last1); last2 = last1 * r2**(N2 - 1)
    if last2 > dTE:      # keep the trailing edge resolved
        r2 = ratio_first(L2, N2, dTE); r2 = 1 / r2   # distribution from the TE backwards
        last2 = dTE
    L3 = R - 1.0
    r3 = ratio_first(L3, N3, last2)
    # wall-normal distribution
    rn = ratio_first(R, Nn, first)
    gN, g1, g2, g3 = rn**(Nn - 1), r1**(N1 - 1), r2**(N2 - 1), r3**(N3 - 1)
    # far-field edges: their own, near-uniform spacing (the LE block spans the
    # front arc from 180 to 135 degrees, the aft block the arc from 135 to 90
    # degrees and the top line to x = 1, the wake block the top line from x = 1)
    t1 = math.radians(135.0)
    LF3 = R - 1.0
    gF1, gF2 = 1.0, 1.0
    gF3 = ratio_first(LF3, N3, 0.25 * LF3 / N3)**(N3 - 1)

    pts = {"A0": (0, 0), "A1": (xs, ys), "A2": (1, 0), "B1": (xs, -ys),
           "F0": (xs - R, 0), "F1": (xs + R * math.cos(t1), R * math.sin(t1)), "F2": (1, R),
           "F3": (R, R), "W": (R, 0),
           "G1": (xs + R * math.cos(t1), -R * math.sin(t1)), "G2": (1, -R), "G3": (R, -R)}
    names = list(pts)
    idx = {}
    verts = []
    for z in (0.0, DZ):
        for n in names:
            idx[(n, z)] = len(verts)
            verts.append("(%.12g %.12g %g)" % (pts[n][0], pts[n][1], z))

    def hexb(a, b, c, e, nx, ny, gxs, gys):
        """gxs, gys: gradings of the two x edges (p0-p1, p3-p2) and y edges (p0-p3, p1-p2)."""
        v = [idx[(p, 0.0)] for p in (a, b, c, e)] + [idx[(p, DZ)] for p in (a, b, c, e)]
        g = [gxs[0], gxs[1], gxs[1], gxs[0], gys[0], gys[1], gys[1], gys[0], 1, 1, 1, 1]
        return "hex (%s) (%d %d 1) edgeGrading (%s)" % (" ".join(map(str, v)), nx, ny,
                                                         " ".join("%.10g" % x for x in g))

    # upper blocks: x along the surface (surface edge, far edge), y wall-normal;
    # lower blocks: x wall-normal, y along the surface
    blocks = [
        hexb("A0", "A1", "F1", "F0", N1, Nn, (g1, gF1), (gN, gN)),      # upper, LE block
        hexb("A1", "A2", "F2", "F1", N2, Nn, (g2, gF2), (gN, gN)),      # upper, aft block
        hexb("A2", "W", "F3", "F2", N3, Nn, (g3, gF3), (gN, gN)),       # upper wake
        hexb("A0", "F0", "G1", "B1", Nn, N1, (gN, gN), (g1, gF1)),      # lower, LE block
        hexb("B1", "G1", "G2", "A2", Nn, N2, (gN, gN), (g2, gF2)),      # lower, aft block
        hexb("A2", "G2", "G3", "W", Nn, N3, (gN, gN), (g3, gF3)),       # lower wake
    ]

    def poly(a, b, xa, xb, sign, n=200):
        s = np.linspace(0, 1, n + 2)[1:-1]
        if xa < 0.3 + 1e-9 and xb <= 0.3 + 1e-9:
            x = xa + (xb - xa) * (1 - np.cos(0.5 * np.pi * s))
        else:
            x = xa + (xb - xa) * s
        return "polyLine %d %d (%s)" % (a, b, " ".join("(%.12g %.12g %%g)" % (xi, sign * yt(xi)) for xi in x))

    edges = []
    for z in (0.0, DZ):
        e = [poly(idx[("A0", z)], idx[("A1", z)], 0, xs, 1),
             poly(idx[("A1", z)], idx[("A2", z)], xs, 1, 1),
             poly(idx[("A0", z)], idx[("B1", z)], 0, xs, -1),
             poly(idx[("B1", z)], idx[("A2", z)], xs, 1, -1)]
        edges += [t.replace("%g", "%g" % z) for t in e]
        tm = math.radians(157.5)
        for sgn, a1 in ((1, "F1"), (-1, "G1")):
            edges.append("arc %d %d (%.12g %.12g %g)" % (idx[("F0", z)], idx[(a1, z)],
                                                        xs + R * math.cos(tm), sgn * R * math.sin(tm), z))
        # aft far edge: arc from 135 to 90 degrees, then the top line to x = 1
        th = np.linspace(t1, 0.5 * np.pi, 40)[1:]
        for sgn, a, b in ((1, "F1", "F2"), (-1, "G1", "G2")):
            edges.append("polyLine %d %d (%s)" % (idx[(a, z)], idx[(b, z)], " ".join(
                "(%.12g %.12g %g)" % (xs + R * math.cos(t), sgn * R * math.sin(t), z) for t in th)))

    def face(a, b):
        return "(%d %d %d %d)" % (idx[(a, 0.0)], idx[(b, 0.0)], idx[(b, DZ)], idx[(a, DZ)])

    write(os.path.join(d, "system/blockMeshDict"), header("dictionary", "blockMeshDict") + """
// NACA 0012 C-mesh: first cell %.3g, wall-normal ratio %.5f over %d cells,
// surface %d + %d cells per side (LE cell %.2g), wake %d cells
scale 1;
vertices
(
    %s
);
blocks
(
    %s
);
edges
(
    %s
);
boundary
(
    airfoil { type wall; faces (%s %s %s %s); }
    farfield { type patch; faces (%s %s %s %s %s %s %s %s); }
    frontAndBack { type empty; faces (%s); }
);
""" % (first, rn, Nn, N1, N2, dLE, N3, "\n    ".join(verts), "\n    ".join(blocks), "\n    ".join(edges),
       face("A0", "A1"), face("A1", "A2"), face("A0", "B1"), face("B1", "A2"),
       face("F0", "F1"), face("F1", "F2"), face("F2", "F3"), face("F3", "W"),
       face("F0", "G1"), face("G1", "G2"), face("G2", "G3"), face("G3", "W"),
       " ".join("(%d %d %d %d) (%d %d %d %d)" % tuple(
           [idx[(p, 0.0)] for p in q][::-1] + [idx[(p, DZ)] for p in q])
           for q in (("A0", "A1", "F1", "F0"), ("A1", "A2", "F2", "F1"), ("A2", "W", "F3", "F2"),
                     ("A0", "F0", "G1", "B1"), ("B1", "G1", "G2", "A2"), ("A2", "G2", "G3", "W")))))
    return 2 * (N1 + N2 + N3) * Nn


def tmr_mesh(d, grid, firstCell):
    """Body-fitted mesh from a NASA TMR C-grid (tmrGrid.py), patches split by topoSet/createPatch."""
    here = os.path.dirname(os.path.abspath(__file__))
    os.makedirs(d, exist_ok=True)
    cmd = [sys.executable, os.path.join(here, "tmrGrid.py"), grid, os.path.join(d, "grid.p2dfmt")]
    if firstCell > 0:
        cmd += ["--firstCell", str(firstCell)]
    print(subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.strip())
    write(os.path.join(d, "system/topoSetDict"), header("dictionary", "topoSetDict") + """
actions
(
    { name bnd;     type faceSet; action new;      source patchToFace; patch defaultFaces; }
    { name fb;      type faceSet; action new;      source normalToFace; normal (0 0 1); cos 0.9; }
    { name fb;      type faceSet; action add;      source normalToFace; normal (0 0 -1); cos 0.9; }
    { name fb;      type faceSet; action subset;   source faceToFace; set bnd; }
    { name airfoil; type faceSet; action new;      source boxToFace; box (-0.001 -0.1 -1) (1.001 0.1 1); }
    { name airfoil; type faceSet; action subset;   source faceToFace; set bnd; }
    { name airfoil; type faceSet; action subtract; source faceToFace; set fb; }
    { name far;     type faceSet; action new;      source faceToFace; set bnd; }
    { name far;     type faceSet; action subtract; source faceToFace; set fb; }
    { name far;     type faceSet; action subtract; source faceToFace; set airfoil; }
);
""")
    write(os.path.join(d, "system/createPatchDict"), header("dictionary", "createPatchDict") + """
pointSync false;
patches
(
    { name airfoil;      patchInfo { type wall; };  constructFrom set; set airfoil; }
    { name farfield;     patchInfo { type patch; }; constructFrom set; set far; }
    { name frontAndBack; patchInfo { type empty; }; constructFrom set; set fb; }
);
""")
    return None


def ibm_mesh(d, h):
    x0, x1, y0 = -0.05, 1.05, 0.1
    r = 1.15
    nc_x, nc_y = int(round((x1 - x0) / h)), int(round(2 * y0 / h))
    xs = [x0 - (R - 0.3) - 0.0, x0, x1, R]
    ys = [-R, -y0, y0, R]
    nL = n_for(xs[1] - xs[0], h, r); nR = n_for(xs[3] - xs[2], h, r); nB = n_for(R - y0, h, r)
    gL = 1 / ratio_first(xs[1] - xs[0], nL, h)**(nL - 1)
    gR = ratio_first(xs[3] - xs[2], nR, h)**(nR - 1)
    gT = ratio_first(R - y0, nB, h)**(nB - 1)
    nxs, gxs = [nL, nc_x, nR], [gL, 1, gR]
    nys, gys = [nB, nc_y, nB], [1 / gT, 1, gT]
    verts = []
    for z in (0.0, DZ):
        for j in range(4):
            for i in range(4):
                verts.append("(%.12g %.12g %g)" % (xs[i], ys[j], z))
    vid = lambda i, j, k: k * 16 + j * 4 + i
    blocks = []
    for j in range(3):
        for i in range(3):
            v = [vid(i, j, 0), vid(i + 1, j, 0), vid(i + 1, j + 1, 0), vid(i, j + 1, 0),
                 vid(i, j, 1), vid(i + 1, j, 1), vid(i + 1, j + 1, 1), vid(i, j + 1, 1)]
            blocks.append("hex (%s) (%d %d 1) simpleGrading (%.10g %.10g 1)"
                          % (" ".join(map(str, v)), nxs[i], nys[j], gxs[i], gys[j]))
    far = []
    for i in range(3):
        far.append("(%d %d %d %d)" % (vid(i, 0, 0), vid(i + 1, 0, 0), vid(i + 1, 0, 1), vid(i, 0, 1)))
        far.append("(%d %d %d %d)" % (vid(i, 3, 0), vid(i, 3, 1), vid(i + 1, 3, 1), vid(i + 1, 3, 0)))
    for j in range(3):
        far.append("(%d %d %d %d)" % (vid(0, j, 0), vid(0, j, 1), vid(0, j + 1, 1), vid(0, j + 1, 0)))
        far.append("(%d %d %d %d)" % (vid(3, j, 0), vid(3, j + 1, 0), vid(3, j + 1, 1), vid(3, j, 1)))
    fb = []
    for j in range(3):
        for i in range(3):
            fb.append("(%d %d %d %d)" % (vid(i, j, 0), vid(i, j + 1, 0), vid(i + 1, j + 1, 0), vid(i + 1, j, 0)))
            fb.append("(%d %d %d %d)" % (vid(i, j, 1), vid(i + 1, j, 1), vid(i + 1, j + 1, 1), vid(i, j + 1, 1)))
    write(os.path.join(d, "system/blockMeshDict"), header("dictionary", "blockMeshDict") + """
// Cartesian mesh: uniform core h = %g over x [%g, %g], y [-%g, %g], graded by %g to the far field
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
    farfield { type patch; faces (%s); }
    frontAndBack { type empty; faces (%s); }
);
""" % (h, x0, x1, y0, y0, r, "\n    ".join(verts), "\n    ".join(blocks), " ".join(far), " ".join(fb)))

    # extruded airfoil surface, beyond the mesh in z, outward normals
    n = 400
    x = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n + 1)))          # LE..TE
    up = [(xi, yt(xi)) for xi in x[::-1]]                           # TE..LE upper
    lo = [(xi, -yt(xi)) for xi in x[1:-1]]                          # LE..TE lower, without ends
    poly = up + lo                                                  # counterclockwise
    z0, z1 = -0.5, DZ + 0.5
    tris = []
    m = len(poly)
    for k in range(m):
        a, b = poly[k], poly[(k + 1) % m]
        tris.append(((a[0], a[1], z0), (b[0], b[1], z0), (b[0], b[1], z1)))
        tris.append(((a[0], a[1], z0), (b[0], b[1], z1), (a[0], a[1], z1)))
    c = (0.5, 0.0)
    for k in range(m):
        a, b = poly[k], poly[(k + 1) % m]
        tris.append(((c[0], c[1], z1), (a[0], a[1], z1), (b[0], b[1], z1)))
        tris.append(((c[0], c[1], z0), (b[0], b[1], z0), (a[0], a[1], z0)))
    out = ["solid naca0012"]
    for t in tris:
        p = np.array(t)
        nrm = np.cross(p[1] - p[0], p[2] - p[0]); nrm /= np.linalg.norm(nrm)
        out.append("facet normal %.9g %.9g %.9g\n outer loop" % tuple(nrm))
        out += [" vertex %.12g %.12g %.12g" % tuple(q) for q in p]
        out.append(" endloop\nendfacet")
    out.append("endsolid naca0012")
    write(os.path.join(d, "constant/triSurface/naca0012.stl"), "\n".join(out) + "\n")
    return (nc_x + nL + nR) * (nc_y + 2 * nB)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--kind", choices=("bf", "ibm"), required=True)
    ap.add_argument("--alpha", type=float, default=0.0)
    ap.add_argument("--mesh", choices=("wf", "wr"), default="wf")
    ap.add_argument("--Nn", type=int, help="bf: wall-normal cells (default 60 wf, 120 wr)")
    ap.add_argument("--grid", help="bf: NASA TMR C-grid (.p2dfmt) instead of the blockMesh C-mesh")
    ap.add_argument("--R", type=float, default=50.0, help="ibm: far-field distance (500 with --grid runs)")
    ap.add_argument("--h", type=float, default=2e-3, help="ibm: core cell size")
    ap.add_argument("--mergeAlpha", type=float, default=0.5)
    ap.add_argument("--secondOrder", action="store_true")
    ap.add_argument("--endTime", type=int, default=8000)
    ap.add_argument("--levels", type=int, default=0, help="ibm: refinement levels around the airfoil")
    ap.add_argument("--nFine", type=int, default=10, help="ibm: finest-level cells from the wall")
    ap.add_argument("--nBand", type=int, default=6, help="ibm: cells of each coarser level beyond the next finer")
    ap.add_argument("--np", type=int, default=1, help="processors")
    ap.add_argument("--wallFunction", choices=("spalding", "spalartAllmaras", "none"), default="spalding")
    ap.add_argument("--wallSample", type=float, default=0.0)
    ap.add_argument("--wallLawDistance", choices=("centroid", "groupCentroid", "cellAverage", "halfCell"), default="centroid")
    ap.add_argument("--noNuTildaWallLaw", action="store_true")
    a = ap.parse_args()
    d = a.case
    ca, sa = math.cos(math.radians(a.alpha)), math.sin(math.radians(a.alpha))
    Uinf = "(%.12g %.12g 0)" % (ca, sa)

    global R
    R = a.R
    if a.kind == "bf" and a.grid:
        nc = tmr_mesh(d, a.grid, 2e-3 if a.mesh == "wf" else 0.0) or 0
        wall = "airfoil"
    elif a.kind == "bf":
        first = 2e-3 if a.mesh == "wf" else 2.5e-5
        nc = bf_mesh(d, first, a.Nn or (60 if a.mesh == "wf" else 120))
        wall = "airfoil"
    else:
        nc = ibm_mesh(d, a.h)
        wall = None

    far = {"U": "freestreamVelocity; freestreamValue uniform %s;" % Uinf,
           "p": "freestreamPressure; freestreamValue uniform 0;",
           "nuTilda": "freestream; freestreamValue uniform %g;" % NUT_INF,
           "nut": "freestream; freestreamValue uniform %g;" % NUT_INF}
    wallbc = {"U": "noSlip;", "p": "zeroGradient;", "nuTilda": "fixedValue; value uniform 0;",
              "nut": "nutUSpaldingWallFunction; value uniform 0;"}
    for name, cls, dim, val in (("U", "volVectorField", "[0 1 -1 0 0 0 0]", Uinf),
                                ("p", "volScalarField", "[0 2 -2 0 0 0 0]", "0"),
                                ("nuTilda", "volScalarField", "[0 2 -1 0 0 0 0]", "%g" % NUT_INF),
                                ("nut", "volScalarField", "[0 2 -1 0 0 0 0]", "%g" % NUT_INF)):
        bc = "    farfield { type %s }\n    frontAndBack { type empty; }\n" % far[name]
        if wall:
            bc += "    %s { type %s }\n" % (wall, wallbc[name])
        write(os.path.join(d, "0.orig" if a.levels else "0", name), header(cls, name, "0")
              + "dimensions %s;\ninternalField uniform %s;\nboundaryField\n{\n%s}\n" % (dim, val, bc))

    write(os.path.join(d, "constant/transportProperties"), header("dictionary", "transportProperties")
          + "transportModel Newtonian;\nnu %g;\n" % NU)
    write(os.path.join(d, "constant/turbulenceProperties"), header("dictionary", "turbulenceProperties")
          + "simulationType RAS;\nRAS { RASModel %s; turbulence on; printCoeffs on; }\n"
          % ("SpalartAllmaras" if a.kind == "bf" else "ibmSpalartAllmaras"))
    if a.kind == "ibm":
        write(os.path.join(d, "constant/ibmProperties"), header("dictionary", "ibmProperties") + """
alphaMin      0.02;
mergeAlpha    %g;
secondOrder   %s;
nSubSamples   32;
nBisect       40;
writeFields   true;
writeForces   true;
%sbodies
{
    airfoil { type triSurfaceMesh; file "naca0012.stl"; wallType noSlip; wallFunction %s; }
}
""" % (a.mergeAlpha, "true" if a.secondOrder else "false",
       ("wallSampleDistance %g;\n" % a.wallSample if a.wallSample > 0 else "")
       + ("wallLawDistance %s;\n" % a.wallLawDistance if a.wallLawDistance != "centroid" else "")
       + ("nuTildaWallLaw false;\n" if a.noNuTildaWallLaw else ""),
       a.wallFunction))

    so = a.kind == "ibm" and a.secondOrder
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
""" % ("cutCellCorrected" if so else "corrected", "cutCellCorrected" if so else "corrected",
       "meshWave" if a.kind == "bf" else "ibmLevelSet"))
    write(os.path.join(d, "system/fvSolution"), header("dictionary", "fvSolution") + """
solvers
{
%s}
SIMPLE
{
    nNonOrthogonalCorrectors 0;
    consistent      no;
    pRefPoint       (20 20 %g);
    pRefValue       0;
    residualControl { p 1e-7; U 1e-8; nuTilda 1e-8; }
}
relaxationFactors
{
    fields    { p 0.3; }
    equations { U 0.7; nuTilda 0.7; }
}
""" % ((
        # cell linking makes linked pairs stiff two-cell systems that a point smoother
        # cannot move together; algebraic pair agglomeration puts each pair into one
        # coarse cell, so the multigrid coarse correction moves the merged value; for
        # the pressure, the incomplete-Cholesky smoother also resolves the stiff links
        '    p { solver GAMG; agglomerator algebraicPair; smoother DICGaussSeidel; tolerance 1e-10; relTol 0.05; }\n'
        '    "(U|nuTilda)" { solver GAMG; agglomerator algebraicPair; smoother GaussSeidel; nSweeps 2;'
        ' tolerance 1e-12; relTol 0.1; }\n') if a.kind == "ibm" else (
        '    p { solver GAMG; smoother GaussSeidel; tolerance 1e-10; relTol 0.05; }\n'
        '    "(U|nuTilda)" { solver smoothSolver; smoother GaussSeidel; nSweeps 2; tolerance 1e-12; relTol 0.1; }\n'),
       0.5 * DZ))
    cd = header("dictionary", "controlDict") + """
application     %s;
%sstartFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         %d;
deltaT          1;
writeControl    timeStep;
writeInterval   %d;
purgeWrite      0;
writeFormat     ascii;
writePrecision  12;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;
""" % ("simpleFoam" if a.kind == "bf" else "ibmSimpleFoam",
       'libs ("libcutCellIBMTurbulence.so");\n' if a.kind == "ibm" else "", a.endTime, a.endTime)
    if a.kind == "bf":
        cd += """functions
{
    forces
    {
        type            forces;
        libs            (forces);
        writeControl    timeStep;
        writeInterval   1;
        patches         (airfoil);
        rho             rhoInf;
        rhoInf          1;
        CofR            (0.25 0 0);
    }
    wallShearStress
    {
        type            wallShearStress;
        libs            (fieldFunctionObjects);
        patches         (airfoil);
        writeControl    writeTime;
    }
}
"""
    write(os.path.join(d, "system/controlDict"), cd)
    if a.levels:
        if a.kind != "ibm":
            raise SystemExit("--levels needs --kind ibm")
        write(os.path.join(d, "system/refineMeshDict"), header("dictionary", "refineMeshDict") + """
set             refine;
coordinateSystem global;
globalCoeffs    { tan1 (1 0 0); tan2 (0 1 0); }
directions      (tan1 tan2);
useHexTopology  yes;
geometricCut    no;
writeMesh       no;
""")
        shutil.copy(os.path.join(os.path.dirname(os.path.abspath(__file__)), "refineSet.py"),
                    os.path.join(d, "system/refineSet.py"))
        shutil.copy(os.path.abspath(__file__), os.path.join(d, "system/makeCase.py"))
        hs = [a.h / 2**L for L in range(a.levels + 1)]
        dist = {a.levels: a.nFine * hs[a.levels]}
        for L in range(a.levels - 1, 0, -1):
            dist[L] = dist[L + 1] + a.nBand * hs[L]
        run = "mkdir -p 0 && blockMesh > log.blockMesh 2>&1 && "
        for L in range(1, a.levels + 1):
            run += ("postProcess -func writeCellCentres -time 0 > log.C%d 2>&1 && "
                    "python3 system/refineSet.py . refine %.10g %.10g > log.refineSet%d 2>&1 && "
                    "refineMesh -overwrite -dict system/refineMeshDict > log.refineMesh%d 2>&1 && "
                    % (L, dist[L], hs[L - 1], L, L))
        run += "cp 0.orig/* 0/ && "
        print("refinement: h %s, bands %s" % (" ".join("%.4g" % x for x in hs),
                                             " ".join("%.4g" % dist[L] for L in sorted(dist))))
    elif a.kind == "bf" and a.grid:
        run = ("plot3dToFoam -2D %g -noBlank grid.p2dfmt > log.plot3dToFoam 2>&1 && topoSet > log.topoSet 2>&1 && "
               "createPatch -overwrite > log.createPatch 2>&1 && " % DZ)
    else:
        run = "blockMesh > log.blockMesh 2>&1 && "
    if a.kind == "bf":
        run += "checkMesh > log.checkMesh 2>&1; simpleFoam > log.simpleFoam 2>&1"
    else:
        run += ("checkMesh > log.checkMesh 2>&1; ibmSetGeometry -writeWallDist > log.geo 2>&1 && "
                "postProcess -func writeCellCentres -time 0 > log.C 2>&1 && ")
        if a.np > 1:
            write(os.path.join(d, "system/decomposeParDict"), header("dictionary", "decomposeParDict")
                  + "numberOfSubdomains %d;\nmethod scotch;\n" % a.np)
            run += ("decomposePar -force > log.decomposePar 2>&1 && mpirun -np %d ibmSimpleFoam -parallel "
                    "> log.ibmSimpleFoam 2>&1 && reconstructPar -latestTime > log.reconstructPar 2>&1" % a.np)
        else:
            run += "ibmSimpleFoam > log.ibmSimpleFoam 2>&1"
    write(os.path.join(d, "Allrun"), '#!/bin/bash\ncd "${0%%/*}" || exit 1\n%s\n' % run)
    os.chmod(os.path.join(d, "Allrun"), 0o755)
    print("%s: %s, alpha %g, about %d cells" % (d, a.kind, a.alpha, nc))


if __name__ == "__main__":
    main()
