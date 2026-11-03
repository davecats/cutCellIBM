#!/usr/bin/env python3
"""Body-fitted periodic cylinder (R=1, box 10x10, centre (5,5)) for the free-slip study.
usage: cylinder_bodyFitted.py <dir> <M> <wall>   wall = noSlip | slip | zeroShear
O-grid: polar ring R=1..2 (orthogonal at the wall), outer blocks to the square."""
import sys, os, math
d, M, wall = sys.argv[1], int(sys.argv[2]), sys.argv[3]
os.makedirs(d + "/system", exist_ok=True); os.makedirs(d + "/constant", exist_ok=True); os.makedirs(d + "/0", exist_ok=True)
c = 5.0; R1, R2 = 1.0, 2.0; dz = 0.1
Nr1 = max(2, round(math.log(R2/R1)/(math.pi/(2*M))))
Nr2 = max(2, M//2)
ang = [45, 135, 225, 315]
def pol(r, a): return (c + r*math.cos(math.radians(a)), c + r*math.sin(math.radians(a)))
sq = [(10, 10), (0, 10), (0, 0), (10, 0)]
pts = [pol(R1, a) for a in ang] + [pol(R2, a) for a in ang] + sq
V = ["(%.12g %.12g 0)" % p for p in pts] + ["(%.12g %.12g %g)" % (p[0], p[1], dz) for p in pts]
n = len(pts)
I = lambda k: k % 4; Mi = lambda k: 4 + k % 4; S = lambda k: 8 + k % 4
blocks, edges = [], []
for k in range(4):
    a, b = k, k + 1
    blocks.append("hex (%d %d %d %d %d %d %d %d) (%d %d 1) simpleGrading (2 1 1)" % (
        I(a), Mi(a), Mi(b), I(b), I(a)+n, Mi(a)+n, Mi(b)+n, I(b)+n, Nr1, M))
    blocks.append("hex (%d %d %d %d %d %d %d %d) (%d %d 1) simpleGrading (2.5 1 1)" % (
        Mi(a), S(a), S(b), Mi(b), Mi(a)+n, S(a)+n, S(b)+n, Mi(b)+n, Nr2, M))
    am = ang[k] + 45
    for r, f in ((R1, I), (R2, Mi)):
        p = pol(r, am)
        edges.append("arc %d %d (%.12g %.12g 0)" % (f(a), f(b), p[0], p[1]))
        edges.append("arc %d %d (%.12g %.12g %g)" % (f(a)+n, f(b)+n, p[0], p[1], dz))
cyl = " ".join("(%d %d %d %d)" % (I(k), I(k)+n, I(k+1)+n, I(k+1)) for k in range(4))
side = lambda k: "(%d %d %d %d)" % (S(k), S(k+1), S(k+1)+n, S(k)+n)
fb = " ".join("(%d %d %d %d)" % (q[0], q[3], q[2], q[1]) + " (%d %d %d %d)" % (q[0]+n, q[1]+n, q[2]+n, q[3]+n)
              for q in [(I(k), Mi(k), Mi(k+1), I(k+1)) for k in range(4)] + [(Mi(k), S(k), S(k+1), Mi(k+1)) for k in range(4)])
H = lambda o, cl="dictionary": "FoamFile { version 2.0; format ascii; class %s; object %s; }\n" % (cl, o)
open(d + "/system/blockMeshDict", "w").write(H("blockMeshDict") + "scale 1;\nvertices (\n" + "\n".join(V) + "\n);\nblocks (\n"
  + "\n".join(blocks) + "\n);\nedges (\n" + "\n".join(edges) + "\n);\nboundary (\n"
  + "cylinder { type wall; faces (%s); }\n" % cyl
  + "top { type cyclic; neighbourPatch bottom; faces (%s); }\n" % side(0)
  + "left { type cyclic; neighbourPatch right; faces (%s); }\n" % side(1)
  + "bottom { type cyclic; neighbourPatch top; faces (%s); }\n" % side(2)
  + "right { type cyclic; neighbourPatch left; faces (%s); }\n" % side(3)
  + "frontAndBack { type empty; faces (%s); }\n);\n" % fb)
# wall condition for U
if wall == "noSlip":
    ubc = "type noSlip;"
elif wall == "slip":      # du_t/dn = 0, u_n = 0: the IBM free slip (zero Laplacian wall flux)
    ubc = "type slip;"
else:                     # zero tangential traction on the circle: du_t/dr = u_t/R
    ubc = """type codedFixedValue; value uniform (0 0 0); name zeroShearWall;
        code #{
            const fvPatch& p = this->patch();
            const vectorField n(p.nf());
            const vectorField Ui(this->patchInternalField());
            const scalarField delta(p.deltaCoeffs());
            const scalar R = 1.0;
            // snGrad along the outward (into the body) normal: (ub - Ui)*delta = -ub/R
            vectorField ut(Ui - (Ui & n)*n);
            operator==(ut/(1.0 + 1.0/(R*delta)));
        #};"""
cyc = "top { type cyclic; } bottom { type cyclic; } left { type cyclic; } right { type cyclic; } frontAndBack { type empty; }"
open(d + "/0/U", "w").write(H("U", "volVectorField") + "dimensions [0 1 -1 0 0 0 0];\ninternalField uniform (0.5 0 0);\nboundaryField {\n cylinder { %s }\n %s\n}\n" % (ubc, cyc))
open(d + "/0/p", "w").write(H("p", "volScalarField") + "dimensions [0 2 -2 0 0 0 0];\ninternalField uniform 0;\nboundaryField {\n cylinder { type zeroGradient; }\n %s\n}\n" % cyc)
open(d + "/constant/transportProperties", "w").write(H("transportProperties") + "transportModel Newtonian;\nnu 1.0;\n")
open(d + "/constant/turbulenceProperties", "w").write(H("turbulenceProperties") + "simulationType laminar;\nlaminar { model Stokes; }\n")
open(d + "/constant/fvOptions", "w").write(H("fvOptions") + "momentumSource { type meanVelocityForce; selectionMode all; fields (U); Ubar (0.5 0 0); relaxation 1.0; }\n")
open(d + "/system/controlDict", "w").write(H("controlDict") + """application simpleFoam; startFrom latestTime; startTime 0; stopAt endTime; endTime 20000;
deltaT 1; writeControl timeStep; writeInterval 20000; purgeWrite 0; writeFormat binary; writePrecision 12; runTimeModifiable false;
functions { wallShear { type wallShearStress; libs (fieldFunctionObjects); patches (cylinder); writeControl onEnd; }
 forces { type forces; libs (forces); patches (cylinder); rho rhoInf; rhoInf 1; CofR (5 5 0); writeControl timeStep; writeInterval 100; } }
""")
open(d + "/system/fvSchemes", "w").write(H("fvSchemes") + """ddtSchemes { default steadyState; } gradSchemes { default Gauss linear; }
divSchemes { default none; div(phi,U) Gauss linear; div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; } interpolationSchemes { default linear; } snGradSchemes { default corrected; }
""")
open(d + "/system/fvSolution", "w").write(H("fvSolution") + """solvers { p { solver GAMG; smoother DIC; tolerance 1e-13; relTol 0.01; }
 U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-13; relTol 0.1; } }
SIMPLE { nNonOrthogonalCorrectors 1; consistent yes; pRefCell 0; pRefValue 0; residualControl { p 1e-10; U 1e-10; } }
relaxationFactors { equations { U 0.9; } fields { p 1.0; } }
""")
