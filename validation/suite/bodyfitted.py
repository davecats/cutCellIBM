"""Body-fitted references for the transient cylinder cases, standard pimpleFoam on blockMesh O-grids
(doc/validationSuite.md §12, decisions D2 and D3):

  rot  V12 problem (Re 200, alpha = a omega/U, square domain of half width 102.4, inflow U = 1,
       outflow p = 0, slip sides), cylinder D = 1 at the origin, rotatingWallVelocity. M cells per
       quarter circle; polar ring 0.5 .. R1 = 4 (first cell a quarter of the wall spacing), outer
       blocks to the square. Statistics over the whole C_L cycles in the last 50 of 150 time units.
  dfg  DFG 2D-2 (Re 100, U_m = 1.5 parabolic): channel [0, 2.2] x [0, 0.41], cylinder D = 0.1 at
       (0.2, 0.2); polar ring to the square [0.1, 0.3]^2, eight blocks around it. Same time step,
       schemes and solver settings as V5 2D-2 at D/h = 40; nOuterCorrectors as given (D2: is one
       PISO pass first order in time without the immersed wall too?).

usage: bodyfitted.py rot <dir> <alpha> <M> <ranks>
       bodyfitted.py dfg <dir> <M> <nOuterCorrectors> <dt factor> <ranks>
       bodyfitted.py eval <dir> ...
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

import common
import meshes
from common import write_dict

DZ = 0.1
SCHEMES = """ddtSchemes { default backward; } gradSchemes { default Gauss linear; }
divSchemes { default none; div(phi,U) Gauss linear; div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; } interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
"""


def graded(a, b, length):
    """cells and blockMesh expansion ratio of a geometric distribution from a to b over length"""
    if abs(b - a) < 1e-12*length:
        return max(1, round(length/a)), 1.0
    q = (length - a)/(length - b)
    return max(1, round(1 + math.log(b/a)/math.log(q))), b/a


class Blocks:
    """2-D blockMeshDict: vertices in the plane z = 0 (copied to z = DZ), hex blocks from
    counter-clockwise quads, arc edges, patches from vertex pairs"""

    def __init__(self):
        self.v, self.blocks, self.edges, self.patches, self.quads = [], [], [], {}, []

    def vertex(self, p):
        self.v.append((float(p[0]), float(p[1])))
        return len(self.v) - 1

    def hex(self, q, n1, n2, g1=1.0, g2=1.0):
        self.quads.append(q)
        self.blocks.append((q, n1, n2, g1, g2))

    def arc(self, a, b, p):
        self.edges.append((a, b, p))

    def face(self, patch, a, b):
        self.patches.setdefault(patch, ("patch", []))[1].append((a, b))

    def write(self, d, types):
        n = len(self.v)
        V = ["(%.15g %.15g 0)" % p for p in self.v] + ["(%.15g %.15g %g)" % (p[0], p[1], DZ) for p in self.v]
        B = ["hex (%s) (%d %d 1) simpleGrading (%.12g %.12g 1)"
             % (" ".join(str(i) for i in list(q) + [k + n for k in q]), n1, n2, g1, g2)
             for q, n1, n2, g1, g2 in self.blocks]
        E = []
        for a, b, p in self.edges:
            E.append("arc %d %d (%.15g %.15g 0)" % (a, b, p[0], p[1]))
            E.append("arc %d %d (%.15g %.15g %g)" % (a + n, b + n, p[0], p[1], DZ))
        P = []
        for name, (_, faces) in self.patches.items():
            P.append("%s { type %s; faces (%s); }" % (name, types.get(name, "patch"),
                     " ".join("(%d %d %d %d)" % (a, b, b + n, a + n) for a, b in faces)))
        fb = " ".join("(%d %d %d %d) (%d %d %d %d)" % (q[0], q[3], q[2], q[1], q[0] + n, q[1] + n, q[2] + n, q[3] + n)
                      for q in self.quads)
        P.append("frontAndBack { type empty; faces (%s); }" % fb)
        write_dict(Path(d) / "system" / "blockMeshDict", "dictionary", "blockMeshDict",
                   "scale 1;\nvertices (\n%s\n);\nblocks (\n%s\n);\nedges (\n%s\n);\nboundary (\n%s\n);\n"
                   % ("\n".join(V), "\n".join(B), "\n".join(E), "\n".join(P)))


def pol(c, r, deg):
    return (c[0] + r*math.cos(math.radians(deg)), c[1] + r*math.sin(math.radians(deg)))


def fv_solution(d, n_outer):
    write_dict(Path(d) / "system" / "fvSolution", "dictionary", "fvSolution", """
solvers { p { solver PCG; preconditioner DIC; tolerance 1e-9; relTol 0.01; } pFinal { $p; relTol 0; }
          U { solver PBiCGStab; preconditioner DILU; tolerance 1e-10; relTol 0; } UFinal { $U; } }
PIMPLE { momentumPredictor yes; nOuterCorrectors %d; nCorrectors 2; nNonOrthogonalCorrectors 1;
         consistent yes; ddtCorr yes; }
""" % n_outer)


def common_files(d, nu, ranks, cofr):
    write_dict(Path(d) / "constant" / "transportProperties", "dictionary", "transportProperties",
               "transportModel Newtonian;\nnu %r;\n" % nu)
    write_dict(Path(d) / "constant" / "turbulenceProperties", "dictionary", "turbulenceProperties",
               "simulationType laminar;\n")
    write_dict(Path(d) / "system" / "fvSchemes", "dictionary", "fvSchemes", SCHEMES)
    write_dict(Path(d) / "system" / "decomposeParDict", "dictionary", "decomposeParDict",
               "numberOfSubdomains %d;\nmethod scotch;\n" % ranks)
    return ("functions { forces { type forces; libs (forces); patches (cylinder); rho rhoInf; rhoInf 1; "
            "CofR (%r %r 0); writeControl timeStep; writeInterval 1; } }" % cofr)


def rot(d, alpha, M, ranks, L=102.4, R0=0.5, R1=4.0, end=150.0):
    d = Path(d)
    b = Blocks()
    ang = [45, 135, 225, 315]
    I = [b.vertex(pol((0, 0), R0, a)) for a in ang]
    Mv = [b.vertex(pol((0, 0), R1, a)) for a in ang]
    S = [b.vertex(p) for p in ((L, L), (-L, L), (-L, -L), (L, -L))]
    hw = 2*math.pi*R0/(4*M)                         # tangential spacing at the wall
    n1, g1 = graded(0.25*hw, 2*math.pi*R1/(4*M), R1 - R0)
    n2, g2 = graded(2*math.pi*R1/(4*M), 2*L/M, L - R1)
    for k in range(4):
        a, c = k, (k + 1) % 4
        b.hex((I[a], Mv[a], Mv[c], I[c]), n1, M, g1)
        b.hex((Mv[a], S[a], S[c], Mv[c]), n2, M, g2)
        for r, V in ((R0, I), (R1, Mv)):
            b.arc(V[a], V[c], pol((0, 0), r, ang[k] + 45))
        b.face("cylinder", I[a], I[c])
    b.face("sides", S[0], S[1])
    b.face("inlet", S[1], S[2])
    b.face("sides", S[2], S[3])
    b.face("outlet", S[3], S[0])
    b.write(d, {"cylinder": "wall"})
    nu = 1/200.0
    fo = common_files(d, nu, ranks, (0.0, 0.0))
    fv_solution(d, 3)
    write_dict(d / "system" / "controlDict", "dictionary", "controlDict", """application pimpleFoam;
startFrom latestTime; startTime 0; stopAt endTime; endTime %r; deltaT %r;
adjustTimeStep yes; maxCo 0.5; maxDeltaT 0.01;
writeControl adjustableRunTime; writeInterval %r; purgeWrite 0; writeFormat ascii; writePrecision 12;
timeFormat general; timePrecision 10; runTimeModifiable false;
%s
""" % (end, 0.1*hw/(1.5 + alpha), end/3, fo))
    meshes.write_fields(d, {
        "U": ("volVectorField", "[0 1 -1 0 0 0 0]", "uniform (1 0 0)",
              {"inlet": "type fixedValue; value uniform (1 0 0);", "outlet": "type zeroGradient;",
               "sides": "type slip;",
               "cylinder": "type rotatingWallVelocity; origin (0 0 0); axis (0 0 1); omega %r; value uniform (0 0 0);"
               % (2.0*alpha)}),
        "p": ("volScalarField", "[0 2 -2 0 0 0 0]", "uniform 0",
              {"inlet": "type zeroGradient;", "outlet": "type fixedValue; value uniform 0;", "sides": "type slip;",
               "cylinder": "type zeroGradient;"})})
    meta = {"kind": "rot", "alpha": alpha, "M": M, "ranks": ranks, "n_ring": n1, "n_outer": n2}
    (d / "case.json").write_text(json.dumps(meta))
    return meta


def dfg(d, M, n_outer, dt_factor, ranks, end=16.0):
    d = Path(d)
    H, Dc, c = 0.41, 0.1, (0.2, 0.2)
    b = Blocks()
    X, Y = (0.0, 0.1, 0.3, 2.2), (0.0, 0.1, 0.3, H)
    G = [[b.vertex((x, y)) for y in Y] for x in X]          # G[ix][iy]
    hs = 0.2/M                                              # spacing on the square
    nL, nB, nT = round(0.1/hs), round(0.1/hs), round((H - 0.3)/hs)
    nR, gR = graded(hs, 0.02, 2.2 - 0.3)
    nx, ny, gx = (nL, M, nR), (nB, M, nT), (1.0, 1.0, gR)
    for ix in range(3):
        for iy in range(3):
            if (ix, iy) != (1, 1):
                b.hex((G[ix][iy], G[ix + 1][iy], G[ix + 1][iy + 1], G[ix][iy + 1]), nx[ix], ny[iy], gx[ix])
    ang = [45, 135, 225, 315]
    I = [b.vertex(pol(c, Dc/2, a)) for a in ang]
    Q = [G[2][2], G[1][2], G[1][1], G[2][1]]
    hw = math.pi*Dc/(4*M)
    nr, gr = graded(0.5*hw, hs, 0.1 - Dc/2)
    for k in range(4):
        a, e = k, (k + 1) % 4
        b.hex((I[a], Q[a], Q[e], I[e]), nr, M, gr)
        b.arc(I[a], I[e], pol(c, Dc/2, ang[k] + 45))
        b.face("cylinder", I[a], I[e])
    for iy in range(3):
        b.face("inlet", G[0][iy], G[0][iy + 1])
        b.face("outlet", G[3][iy], G[3][iy + 1])
    for ix in range(3):
        b.face("walls", G[ix][0], G[ix + 1][0])
        b.face("walls", G[ix][3], G[ix + 1][3])
    b.write(d, {"cylinder": "wall", "walls": "wall"})
    fo = common_files(d, 1e-3, ranks, c)
    fv_solution(d, n_outer)
    # the time step of V5 2D-2 at D/h = 40 (h = 0.0025, Courant 0.8 with |u| <= 1.5 U_m), times dt_factor
    dt = 0.8*0.0025/(1.5*1.5)
    dt = 0.33/math.ceil(0.33/dt)*dt_factor
    write_dict(d / "system" / "controlDict", "dictionary", "controlDict", """application pimpleFoam;
startFrom latestTime; startTime 0; stopAt endTime; endTime %r; deltaT %r;
writeControl timeStep; writeInterval %d; purgeWrite 0; writeFormat ascii; writePrecision 12;
timeFormat general; timePrecision 10; runTimeModifiable false;
%s
""" % (end, dt, int(round(end/dt/4)), fo))
    common.sh("blockMesh", d, "log.blockMesh")
    y = meshes.face_centres(d)["inlet"][:, 1]
    Uin = np.zeros((len(y), 3))
    Uin[:, 0] = 4*1.5*y*(H - y)/H**2
    meshes.write_fields(d, {
        "U": ("volVectorField", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)",
              {"inlet": "type fixedValue; value %s;" % meshes.vector_list(Uin), "walls": "type noSlip;",
               "cylinder": "type noSlip;", "outlet": "type zeroGradient;"}),
        "p": ("volScalarField", "[0 2 -2 0 0 0 0]", "uniform 0",
              {"inlet": "type zeroGradient;", "walls": "type zeroGradient;", "cylinder": "type zeroGradient;",
               "outlet": "type fixedValue; value uniform 0;"})})
    meta = {"kind": "dfg", "M": M, "nOuter": n_outer, "dt": dt, "ranks": ranks, "n_ring": nr, "n_right": nR}
    (d / "case.json").write_text(json.dumps(meta))
    return meta


def run(d):
    d = Path(d)
    ranks = json.loads((d / "case.json").read_text())["ranks"]
    common.sh("blockMesh", d, "log.blockMesh")
    common.sh("checkMesh", d, "log.checkMesh")
    common.sh("decomposePar -force", d, "log.decomposePar")
    common.sh("pimpleFoam", d, "log.pimpleFoam", np_=ranks)


def history(d):
    """time, C_D, C_L from the forces function object (U = 1 and D = 1, or mean U = 1 and D = 0.1),
    and C_M = -M_z/(U^2 a^2) about CofR (a = D/2; Mittal & Kumar's sign: positive against the rotation)"""
    d = Path(d)
    meta = json.loads((d / "case.json").read_text())
    D = 1.0 if meta["kind"] == "rot" else 0.1
    rows = []
    for f in sorted(d.glob("postProcessing/forces/*/force.dat"), key=lambda p: float(p.parent.name)):
        # v2106: "time (total) (pressure) (viscous)"
        rows.append(np.loadtxt((l.replace("(", " ").replace(")", " ") for l in f.open()), ndmin=2, comments="#"))
    h = np.vstack(rows)
    _, k = np.unique(h[:, 0], return_index=True)
    h = h[k]
    F = h[:, 1:4]
    rows = []
    for f in sorted(d.glob("postProcessing/forces/*/moment.dat"), key=lambda p: float(p.parent.name)):
        rows.append(np.loadtxt((l.replace("(", " ").replace(")", " ") for l in f.open()), ndmin=2, comments="#"))
    m = np.vstack(rows)
    _, k = np.unique(m[:, 0], return_index=True)
    m = m[k]
    assert np.array_equal(m[:, 0], h[:, 0]), "force.dat and moment.dat times differ"
    return h[:, 0], 2*F[:, 0]/(D*DZ), 2*F[:, 1]/(D*DZ), -m[:, 3]/((D/2)**2*DZ)


def evaluate(d):
    d = Path(d)
    meta = json.loads((d / "case.json").read_text())
    t, cd, cl, cm = history(d)
    out = dict(meta, t_end=float(t[-1]), cells=None)
    m = __import__("re").search(r"cells:\s+(\d+)", (d / "log.checkMesh").read_text()) if (d / "log.checkMesh").exists() else None
    out["cells"] = int(m.group(1)) if m else None
    if meta["kind"] == "rot":
        st = meshes.periodic_statistics(t, cl, cd, 50.0, cm)
        out.update(St=st.pop("f"), **st)
    else:
        import reference_dfg as dfg_ref
        s = t > t[-1] - 3.0
        out.update({k: float(v) for k, v in dfg_ref.cycle_statistics(t[s], cd[s], cl[s]).items()})
    return out


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "rot":
        print(rot(a[1], float(a[2]), int(a[3]), int(a[4])))
        run(a[1])
        print(json.dumps(evaluate(a[1])))
    elif a[0] == "dfg":
        print(dfg(a[1], int(a[2]), int(a[3]), float(a[4]), int(a[5])))
        run(a[1])
        print(json.dumps(evaluate(a[1])))
    elif a[0] == "eval":
        for x in a[1:]:
            print(json.dumps(evaluate(x)))
