"""V4 passive scalar at a static wall (tutorials/scalar/steadyCylinder_fixedValue): periodic
cylinder Re_D = 1 (Ubar 0.5), D_T = 0.1 (Pe = 10), uniform source of rate 1 per unit volume in the
disk r < 0.5 around (5, 8.5), T = 0 on the cylinder. Steady scalar content per unit source rate,
immersed (N = 41...161) against the body-fitted O-grid (validation/cylinder_bodyFitted.py, M = 64,
128) with the scalarTransport function object and a semi-implicit source in the same disk
(cells whose centres lie in it, in both)."""
import os
import re
import subprocess

import numpy as np

import common
from common import write_dict
from foamio import latest_time, read_field

ID, KIND, TIER = "V4", "case", 1
TITLE = "passive scalar at a static wall vs body-fitted"
STRESSES, RELIES = ["K11"], ["G6", "V2"]
REGRESSION = (1e-5,)

SRC = common.TUTORIALS / "scalar" / "steadyCylinder_fixedValue"
NUM = r"([-+0-9.eE]+)"


def ibm(ctx, N):
    d = common.copy_case(SRC, ctx.root / ("ibm_N%d" % N))
    common.set_keys(d / "system" / "blockMeshDict", nx=N, ny=N)
    common.set_keys(d / "system" / "controlDict", endTime=40000, writeInterval=40000)
    app = re.search(r"application\s+(\w+)", (d / "system" / "controlDict").read_text()).group(1)
    for cmd in ("blockMesh", "ibmSetGeometry", app):
        common.sh(cmd, d, "log." + cmd)
    l = [x for x in (d / ("log." + app)).read_text().splitlines() if x.startswith("scalar T:")][-1]
    c = float(re.search(r"fluid content " + NUM, l).group(1))
    r = float(re.search(r"source rate " + NUM, l).group(1))
    return {"kind": "ibm", "N": N, "content": c, "rate": r, "ratio": c/r}


def bf(ctx, M):
    d = ctx.root / ("bf_M%d" % M)
    if not (common.REUSE and (d / "system").exists()):
        subprocess.run([common.sys.executable, str(common.VALIDATION / "cylinder_bodyFitted.py"), str(d), str(M),
                        "noSlip"], check=True, capture_output=True)
        ctl = (d / "system" / "controlDict").read_text().replace("writeFormat binary", "writeFormat ascii")
        ctl = ctl.replace("endTime 20000;", "endTime 40000;").replace("writeInterval 20000;", "writeInterval 40000;")
        ctl = ctl.replace("functions {", """functions { T { type scalarTransport; libs (solverFunctionObjects); field T;
  D 0.1; nCorr 0; schemesField T;
  fvOptions { src { type scalarSemiImplicitSource; selectionMode cellSet; cellSet src; volumeMode specific;
                    injectionRateSuSp { T (1 0); } } } }
""")
        (d / "system" / "controlDict").write_text(ctl)
        t = (d / "system" / "fvSchemes").read_text().replace("div((nuEff", "div(phi,T) Gauss linear; div((nuEff")
        (d / "system" / "fvSchemes").write_text(t)
        t = (d / "system" / "fvSolution").read_text().replace(
            "solvers {", "solvers { T { solver PBiCGStab; preconditioner DILU; tolerance 1e-13; relTol 0; }")
        t = t.replace("residualControl { p 1e-10; U 1e-10; }", "residualControl { p 1e-10; U 1e-10; T 1e-10; }")
        (d / "system" / "fvSolution").write_text(t)
        cyc = "top { type cyclic; } bottom { type cyclic; } left { type cyclic; } right { type cyclic; } frontAndBack { type empty; }"
        write_dict(d / "0" / "T", "volScalarField", "T", "dimensions [0 0 0 1 0 0 0];\ninternalField uniform 0;\n"
                   "boundaryField {\n cylinder { type fixedValue; value uniform 0; }\n %s\n}\n" % cyc)
        write_dict(d / "system" / "topoSetDict", "dictionary", "topoSetDict",
                   "actions ( { name src; type cellSet; action new; source cylinderToCell; "
                   "p1 (5 8.5 -1); p2 (5 8.5 1); radius 0.5; } );\n")
    for cmd in ("blockMesh", "topoSet", "postProcess -func writeCellVolumes -time 0", "simpleFoam"):
        common.sh(cmd, d, "log." + cmd.split()[0])
    t = latest_time(str(d))
    T, _ = read_field(os.path.join(d, t, "T"))
    V, _ = read_field(os.path.join(d, "0", "V"))
    S = sum(V[i] for i in map(int, re.findall(r"\d+", (d / "constant" / "polyMesh" / "sets" / "src").read_text()
                                                .split("(", 1)[1].split(")")[0])))
    return {"kind": "bf", "M": M, "content": float(np.sum(T*V)), "rate": float(S), "ratio": float(np.sum(T*V)/S)}


def run(ctx, res):
    jobs = [(ibm, N) for N in ([41, 81] if ctx.quick else [41, 81, 161])] + [(bf, M) for M in (64, 128)]
    rows = common.pmap(lambda j: j[0](ctx, j[1]), jobs, ctx.jobs)
    res.data["rows"] = rows
    res.data["metrics"] = {"%s_%d_ratio" % (r["kind"], r.get("N", r.get("M"))): r["ratio"] for r in rows}
    ib = sorted([r for r in rows if r["kind"] == "ibm"], key=lambda r: r["N"])
    b = sorted([r for r in rows if r["kind"] == "bf"], key=lambda r: r["M"])
    res.report("body-fitted M = 64 vs 128 (content/rate)", b[0]["ratio"], b[1]["ratio"],
               "deviation %+.3f %%" % (100*(b[0]["ratio"]/b[1]["ratio"] - 1)))
    res.rel("scalar content per unit source rate, N = %d vs body-fitted M = 128" % ib[-1]["N"], ib[-1]["ratio"],
            b[-1]["ratio"], 0.01, " ".join("%.5g" % r["ratio"] for r in ib))
