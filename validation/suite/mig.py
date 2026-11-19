"""C3 driver of V6 MIG (doc/validationSuite.md §8): a sphere of diameter d in a periodic circular
pipe d_h = 1 (Rivero-Rodriguez & Scheid 2018), mean flow J = 1, in the frame of the sphere.

Geometry: pipe axis x through the origin, sphere centre (0, e, 0), e = eps the transverse position;
half domain z >= 0 with a symmetry plane through the centre; x periodic, length L, the mesh
mirror-symmetric in x about the sphere centre. Background h0 = 2 hf (cubic), one level of refinement
to hf in a box around the sphere and its gap to the wall (pipe wall and outer box at h0). The pipe is
an inverted IBM cylinder sliding at -V along its axis (`invariant yes`, C2); the sphere is no slip
turning at Omega e_z about its centre (`invariant yes`, rigid) or free slip (zero shear, stress-free).
ibmMeanVelocityForce imposes the flux A_f (J - V) through every section: the fluid-volume average
Ubar = (J - V) Vf0/Vf (Vf0, Vf: discrete fluid volumes without and with the sphere).

Force balance along x on the sphere: F_x = (pressure + viscous, IBM) + G V_B (the mean pressure
gradient on the sphere volume; G is the driving gradient along x). f = L_y/V_B in the paper's
units mu J/d_h^2, so f/Re = L_y/V_B in units rho = J = d_h = 1 (mu = 1/Re); beta from
G/G0 - 1 = beta alpha_G (G0: the empty pipe on the same mesh, alpha_G = V_B/(L pi/4)). balance()
solves F_x = 0 (and T_z = 0, rigid) for V (and Omega) by Newton with a finite-difference Jacobian,
then Broyden updates; each solve restarts from the previous one on the same mesh.

Units rho = J = d_h = 1, nu = 1/Re. All forces and volumes of the half domain are doubled."""
import json
import os
import re
import shutil
from pathlib import Path

import numpy as np

import common
import geometry
import meshes
from common import write_dict
from foamio import latest_time, read_field

J, RP = 1.0, 0.5          # mean flow, pipe radius
CYL = common.TUTORIALS / "cylinderRe1"


def ranks(cells):
    """~40k cells per rank, 4 ... 64 (MIG_MAXRANKS in the environment lowers the cap)"""
    return int(min(int(os.environ.get("MIG_MAXRANKS", 64)), max(4, 2**round(np.log2(cells/40000.0)))))


def mesh(d, hf, dsph=None, e=0.0, L=3.0, refine=True, h0=None):
    """background h0 (default 2 hf), box [-L/2, L/2] x [-W, W] x [0, W], W = 0.5 + 2 h0; nested
    levels to hf in boxes around the sphere, the box of the level with cells h_k (parent 2 h_k)
    reaching 8 h_k beyond the sphere in x and z, in y from e - R - 8 h_k to the outer box (one
    level: R + 4 h0, the pipe wall in the box near the sphere). Returns cells."""
    h0 = h0 or 2*hf
    W = RP + 2*h0
    lo, hi = (-L/2, -W, 0.0), (L/2, W, W)
    boxes = []
    if refine and dsph:
        R = dsph/2
        n = int(round(np.log2(h0/hf)))
        assert abs(h0/2**n - hf) < 1e-9*hf, "h0/hf not a power of 2"
        for k in range(1, n + 1):
            m = 8*h0/2**k
            boxes.append(((-R - m, max(-W, e - R - m), -1.0), (R + m, W + 1.0, R + m)))
    patches = {"left": ("cyclic", ["xmin"], "neighbourPatch right;"),
               "right": ("cyclic", ["xmax"], "neighbourPatch left;"),
               "sym": ("symmetryPlane", ["zmin"]), "outer": ("wall", ["ymin", "ymax", "zmax"])}
    return meshes.refined(d, lo, hi, h0, boxes, patches, two_d=False)


def setup(d, Re, hf, dsph=None, e=0.0, interface="rigid", V=0.0, Om=0.0, L=3.0, ubar=None, sphere=True, h0=None,
          **ibm):
    """case d: the pipe sliding at -V (static if V = 0), the mesh refined around a sphere of
    diameter dsph at e (none: the plain pipe mesh); sphere False: the empty pipe on that mesh; ubar:
    the imposed fluid-volume average (default J - V, the empty pipe)"""
    d = Path(d)
    if not (d / "constant" / "polyMesh" / "owner").exists():
        mesh(d, hf, dsph, e, L, h0=h0)
    pipe = {"type": "cylinder", "centre": (0, 0, 0), "axis": (1, 0, 0), "radius": RP, "invert": True}
    if V:
        pipe["motion"] = {"type": "translating", "velocity": (-V, 0, 0), "invariant": True}
    bodies = {"pipe": pipe}
    if dsph and sphere:
        s = {"type": "sphere", "centre": (0, e, 0), "radius": dsph/2,
             "wallType": "noSlip" if interface == "rigid" else "freeSlip"}
        if interface == "rigid":
            s["motion"] = {"type": "none", "angularVelocity": (0, 0, Om), "invariant": True}
        bodies["sphere"] = s
    geometry.ibm_properties(d, bodies, **ibm)
    for x in ("transportProperties", "turbulenceProperties"):
        (d / "constant" / x).write_text((CYL / "constant" / x).read_text())
    common.set_keys(d / "constant" / "transportProperties", nu=1.0/Re)
    ub = (J - V) if ubar is None else ubar
    write_dict(d / "constant" / "fvOptions", "dictionary", "fvOptions", """
momentumSource { type ibmMeanVelocityForce; selectionMode all; fields (U); Ubar (%r 0 0); relaxation 1.0; }
""" % ub)
    bU = "left { type cyclic; } right { type cyclic; } sym { type symmetryPlane; } outer { type noSlip; }"
    bp = "left { type cyclic; } right { type cyclic; } sym { type symmetryPlane; } outer { type zeroGradient; }"
    (d / "0").mkdir(exist_ok=True)
    if not (d / "0" / "U").exists():
        write_dict(d / "0" / "U", "volVectorField", "U",
                   "dimensions [0 1 -1 0 0 0 0];\ninternalField uniform (%r 0 0);\nboundaryField { %s }\n" % (ub, bU))
        write_dict(d / "0" / "p", "volScalarField", "p",
                   "dimensions [0 2 -2 0 0 0 0];\ninternalField uniform 0;\nboundaryField { %s }\n" % bp)
    (d / "system" / "fvSchemes").write_text((CYL / "system" / "fvSchemes").read_text())
    # p reference far from the sphere, on the other side of the axis
    write_dict(d / "system" / "fvSolution", "dictionary", "fvSolution", """
solvers { p { solver PCG; preconditioner { preconditioner linkedGroup; base { preconditioner GAMG;
                smoother DICGaussSeidel; nCellsInCoarsestLevel 10; tolerance 0; relTol 0.01; } }
              tolerance 1e-13; relTol 0.01; maxIter 1000; }
          U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-13; relTol 0.1; maxIter 50; } }
SIMPLE { nNonOrthogonalCorrectors 0; consistent yes; pRefPoint (%r -0.21 0.013); pRefValue 0;
         residualControl { p 1e-9; U 1e-9; } }
relaxationFactors { equations { U 0.9; } fields { p 1.0; } }
""" % (-0.4*L,))
    geometry.control_dict(d, "ibmSimpleFoam", endTime=20000, deltaT=1, writeInterval=20000)


def fluid_volume(d):
    """discrete fluid volume of the half domain, sum ibmAlpha V"""
    if not (d / "0" / "V").exists():
        common.sh("postProcess -func writeCellVolumes -time 0", d, "log.writeCellVolumes")
    a, _ = read_field(os.path.join(d, "0", "ibmAlpha"))
    v, _ = read_field(os.path.join(d, "0", "V"))
    return float(a @ v)


def geometry_only(d):
    common.sh("postProcess -func writeCellCentres -time 0", d, "log.writeCellCentres")
    common.sh("ibmSetGeometry", d, "log.ibmSetGeometry")
    return fluid_volume(d)


def solve(d, np_=None):
    cells = common.last_match((d / "log.checkMesh").read_text(), r"cells:\s+(\d+)", cast=int)
    np_ = np_ or ranks(cells)
    if np_ > 1:
        write_dict(d / "system" / "decomposeParDict", "dictionary", "decomposeParDict",
                   "numberOfSubdomains %d;\nmethod scotch;\n%s" % (np_, common.LINK_CONSTRAINT))
        common.link_zone(d)
        common.sh("decomposePar -force", d, "log.decomposePar")
        common.sh("ibmSimpleFoam", d, "log.ibmSimpleFoam", np_=np_)
        common.sh("reconstructPar -latestTime", d, "log.reconstructPar")
    else:
        common.sh("ibmSimpleFoam", d, "log.ibmSimpleFoam")
    return extract(d)


def extract(d):
    """half-domain values of the last iteration: G along x, sphere force and moment, convergence"""
    log = (d / "log.ibmSimpleFoam").read_text()
    ub = float(re.search(r"Ubar \(([-0-9.eE+]+)", (d / "constant" / "fvOptions").read_text()).group(1))
    g, change = common.driving_gradient(log)
    # converged: the residual control (1e-9; the zero-shear free slip has a floor of the U residual
    # at ~2e-10) or a stationary iteration (G to 1e-10 over 200 iterations)
    out = {"G": float(np.sign(ub)*g[-1]), "Gchange": change(200),
           "converged": "SIMPLE solution converged" in log or change(200) < 1e-10,
           "iterations": common.last_match(log, r"\nTime = (\d+)", cast=int)}
    if "IBM force on sphere" in log:
        f = common.ibm_force(log, "sphere")
        out["F"] = f["total"].tolist()
        out["Fp"], out["Fv"] = f["pressure"].tolist(), f["viscous"].tolist()
        try:
            out["T"] = common.ibm_moment(log, "sphere")["total"].tolist()
        except ValueError:
            out["T"] = [0.0, 0.0, 0.0]
    return out


def restart_from(d, src):
    t = latest_time(str(src))
    for f in ("U", "p"):
        shutil.copy(Path(src) / t / f, d / "0" / f)


class Point:
    """one configuration (interface, d, e, Re, mesh hf): solves at given (V, Omega), the empty-pipe
    reference G0, the force balance"""

    def __init__(self, root, interface, dsph, e, Re, hf, L=3.0, tag="", h0=None, **ibm):
        self.root, self.interface, self.dsph, self.e, self.Re, self.hf, self.L = Path(root), interface, dsph, e, Re, hf, L
        self.h0 = h0 if h0 and abs(h0 - 2*hf) > 1e-12 else None
        self.ibm = ibm
        self.VB = np.pi*dsph**3/6
        self.name = "%s_d%g_e%.4f_Re%g_h%d%s%s" % (interface, dsph, e, Re, round(1/hf),
                                                  "_H%d" % round(1/self.h0) if self.h0 else "", tag)
        self.runs = []

    def empty(self):
        """G0 and Vf0 of the empty pipe on the mesh of this point (static wall, Ubar = J)"""
        d = self.root / (self.name + "_pipe")
        f = d / "result.json"
        if f.exists():
            return json.loads(f.read_text())
        d.mkdir(parents=True, exist_ok=True)
        setup(d, self.Re, self.hf, self.dsph, self.e, L=self.L, sphere=False, h0=self.h0, **self.ibm)
        # the mesh with the refinement box of the sphere
        if not (d / "log.refineMesh.0").exists():
            raise RuntimeError("empty pipe without the refinement box")
        Vf0 = geometry_only(d)
        r = dict(solve(d), Vf0=Vf0)
        f.write_text(json.dumps(r))
        return r

    def run(self, V, Om=0.0, k=None):
        V, Om = float(V), float(Om)
        k = len(self.runs) if k is None else k
        d = self.root / self.name / ("s%02d" % k)
        f = d / "result.json"
        if f.exists():
            r = json.loads(f.read_text())
            # a cached solve is found by its index only: it must be the state asked for (a replay
            # along another path, e.g. a chain point run alone, would pair solves and states wrongly)
            if abs(r["V"] - V) > 1e-9*max(1.0, abs(V)) or abs(r["Om"] - Om) > 1e-9*max(1.0, abs(Om)):
                raise RuntimeError("%s: cached V, Omega %r, %r but %r, %r asked for"
                                   % (d, r["V"], r["Om"], V, Om))
        else:
            ref = self.empty()
            d.mkdir(parents=True, exist_ok=True)
            # the mesh of the empty pipe (identical), then the case with the sphere
            if not (d / "constant" / "polyMesh").exists():
                shutil.copytree(self.root / (self.name + "_pipe") / "constant" / "polyMesh", d / "constant" / "polyMesh")
                shutil.copytree(self.root / (self.name + "_pipe") / "system", d / "system")
                shutil.copy(self.root / (self.name + "_pipe") / "log.checkMesh", d / "log.checkMesh")
            setup(d, self.Re, self.hf, self.dsph, self.e, self.interface, V, Om, self.L, h0=self.h0, **self.ibm)
            Vf = geometry_only(d)
            ub = float((J - V)*ref["Vf0"]/Vf)
            (d / "constant" / "fvOptions").write_text(re.sub(r"Ubar \([^)]*\)", "Ubar (%r 0 0)" % ub,
                                                            (d / "constant" / "fvOptions").read_text()))
            prev = [x for x in self.runs if x.get("dir")]
            if prev:
                restart_from(d, prev[-1]["dir"])
            r = dict(solve(d), V=V, Om=Om, Vf=Vf, Ubar=ub, dir=str(d))
            f.write_text(json.dumps(r))
        self.runs.append(r)
        return r

    def residual(self, r):
        """full-sphere F_x/V_B (with the mean gradient) and T_z/V_B"""
        Fx = 2*r["F"][0]/self.VB + r["G"]
        res = [Fx]
        if self.interface == "rigid":
            res.append(2*r["T"][2]/self.VB)
        return np.array(res)

    def balance(self, x0, tol, dx=None, Jac=None, maxit=8):
        """Newton on x = (V, Omega) (rigid) or (V,) from x0; Jac: a Jacobian to start from (else
        finite differences with steps dx); returns the result at the balanced state and the Jacobian"""
        x = np.array(x0, float)
        n = len(x)
        r = self.run(*x)
        F = self.residual(r)
        if Jac is None:
            Jac = np.zeros((n, n))
            for i in range(n):
                xi = x.copy()
                xi[i] += dx[i]
                Jac[:, i] = (self.residual(self.run(*xi)) - F)/dx[i]
        for _ in range(maxit):
            if np.all(np.abs(F) < tol):
                break
            s = -np.linalg.solve(Jac, F)
            x = x + s
            r = self.run(*x)
            Fn = self.residual(r)
            Jac = Jac + np.outer(Fn - F - Jac @ s, s)/(s @ s)   # Broyden
            F = Fn
        ref = self.empty()
        aG = self.VB/(self.L*np.pi/4)
        out = {"interface": self.interface, "d": self.dsph, "e": self.e, "eps": self.e/(0.5*(1 - self.dsph)),
               "Re": self.Re, "hf": self.hf, "L": self.L, "V": float(x[0]), "Om": float(x[1]) if n > 1 else None,
               "f": 2*r["F"][1]/self.VB*self.Re, "fRe": 2*r["F"][1]/self.VB, "residual": F.tolist(),
               "balanced": bool(np.all(np.abs(F) < tol)), "solves": len(self.runs), "G": r["G"], "G0": ref["G"],
               "beta_aG": r["G"]/ref["G"] - 1, "beta": (r["G"]/ref["G"] - 1)/aG,
               "converged": bool(r["converged"] and ref["converged"]), "Gchange": max(r["Gchange"], ref["Gchange"]),
               "iterations": r["iterations"], "name": self.name}
        return out, Jac
