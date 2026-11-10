"""G9 stagnation flow onto a plane wall, exact NS solution (validation/reference_hiemenz.py):
Hiemenz (2-D, box [-0.6, 1.4] x [-0.5, 1]) at cut fractions f = 0.1, 0.5, 0.9 of a level wall and on walls tilted at
tan(alpha) = 1/3, 1; Homann (3-D) level and tilted. Convection at the wall (K7) and the wall
traction (K3) without curvature or separation (doc/validationSuite.md §11).

a = 1, nu = 0.01; U from the exact solution on every outer face (cut faces: at the centroid of their wet part), p zeroGradient (fixedFluxPressure is evaluated before the solver sets it) with the
exact value in a reference cell; steady ibmSimpleFoam. Errors at the fluid centroids
(alpha-weighted L2, Linf) of u and p; the wall force against the exact integral over the wall
polygon inside the box (viscous: the tangential part against the shear, the normal part
nu u_n/d reported; pressure: p0 - a^2 s^2/2, with the mean pressure error removed as in the norms). Body-fitted control:
level wall as a mesh patch."""
import math
import os

import numpy as np

import common
import geometry
import reference_hiemenz as hz
from common import write_dict
from foamio import latest_time, read_field

ID, KIND, TIER = "G9", "gate", 0
TITLE = "stagnation flow onto a plane wall (Hiemenz 2-D, Homann 3-D), exact"
STRESSES, RELIES = ["K7", "K3", "K5", "K6"], ["G1", "G2"]
REGRESSION = (0.01,)

A, NU = 1.0, 0.01
CYL = common.TUTORIALS / "cylinderRe1"


def config(kind, tan, f, N, bf=False):
    """box, origin, unit normal of the wall"""
    al = math.atan(tan)
    nrm = np.array([-math.sin(al), math.cos(al), 0.0]) if kind == "hiemenz" else \
        np.array([0.0, -math.sin(al), math.cos(al)])
    if kind == "hiemenz":
        # off-centre in x: the net wall shear of the box does not vanish by symmetry
        lo, hi = np.array([-0.6, -0.5, 0.0]), np.array([1.4, 1.0, 0.0])
        n = (2*N, int(round(1.5*N)), 1)
        hi[2] = 1.0/N
        ax = 1
    else:
        lo, hi = np.array([-0.4, -0.45, -0.4]), np.array([0.6, 0.55, 0.6])
        n = (N, N, N)
        ax = 2
    h = 1.0/N
    origin = np.zeros(3)
    if tan == 0:
        origin[ax] = -f*h       # the wall cuts the cell [-h, 0] with fluid fraction f
    else:
        origin[:2] = (0.0137, -0.0213)
    if bf:                      # body-fitted: the wall is the bottom patch
        lo[ax], origin[ax] = 0.0, 0.0
        n = tuple(int(round((hi[i] - lo[i])/h)) if i == ax else n[i] for i in range(3))
    return lo, hi, n, origin, nrm, ax


def write_case(d, kind, tan, f, N, bf, ma, conv):
    lo, hi, n, origin, nrm, ax = config(kind, tan, f, N, bf)
    sides = ["xmin", "xmax", "ymin", "ymax"] + ([] if kind == "hiemenz" else ["zmin", "zmax"])
    bottom = ["ymin", "zmin"][ax - 1]
    patches = {"outer": ("patch", [s for s in sides if not (bf and s == bottom)])}
    if bf:
        patches["wall"] = ("wall", [bottom])
    geometry.box(d, lo, hi, n, patches)
    bodies = {} if bf else {"wall": {"type": "plane", "point": origin, "normal": nrm, "wallType": "noSlip"}}
    geometry.ibm_properties(d, bodies, mergeAlpha=ma, ibmConvection=conv, nSubSamples=32 if kind == "hiemenz" else 16)
    for x in ("transportProperties", "turbulenceProperties"):
        (d / "constant" / x).write_text((CYL / "constant" / x).read_text())
    common.set_keys(d / "constant" / "transportProperties", nu=NU)
    for x in ("fvSchemes",):
        (d / "system" / x).write_text((CYL / "system" / x).read_text())
    ref = origin + 0.35*nrm + np.array([0.11, 0.07, 0.0 if kind == "hiemenz" else 0.05])
    _, pref = hz.exact(ref[None, :], A, NU, origin, nrm, kind)
    write_dict(d / "system" / "fvSolution", "dictionary", "fvSolution", """
solvers
{
    p { solver PCG; preconditioner DIC; tolerance 1e-13; relTol 0.01; }
    U { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-13; relTol 0.1; }
}
SIMPLE
{
    nNonOrthogonalCorrectors 0; consistent yes;
    pRefPoint %s; pRefValue %r;
    residualControl { p 1e-10; U 1e-10; }
}
relaxationFactors { equations { U 0.9; } fields { p 1.0; } }
""" % (common.foam_value(ref), float(pref[0])))
    geometry.control_dict(d, "ibmSimpleFoam", endTime=20000, deltaT=1, writeInterval=20000)
    return origin, nrm


def face_areas(xf, lo, hi, n, kind):
    """outward area vectors of boundary faces of the uniform box from their centres"""
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    h = (hi - lo)/np.asarray(n, float)
    S = np.zeros_like(xf)
    dims = 2 if kind == "hiemenz" else 3
    for i in range(dims):
        area = np.prod([h[j] for j in range(3) if j != i])
        S[np.abs(xf[:, i] - lo[i]) < 1e-9*(hi[i] - lo[i]), i] = -area
        S[np.abs(xf[:, i] - hi[i]) < 1e-9*(hi[i] - lo[i]), i] = area
    return S


def wet_centres(xf, lo, hi, n, origin, nrm):
    """centroid and fraction of the part of each boundary face on the fluid side of the plane
    wall: the face rectangle (axis-aligned box face) clipped by (x - origin).nrm >= 0"""
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    h = (hi - lo)/np.asarray(n, float)
    xw, th = xf.copy(), np.ones(len(xf))
    for m, c in enumerate(xf):
        i = int(np.argmin([min(abs(c[a] - lo[a]), abs(c[a] - hi[a]))/h[a] for a in range(3)]))
        j, k = [a for a in range(3) if a != i]
        poly = []
        for sj, sk in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
            q = c.copy()
            q[j] += 0.5*sj*h[j]
            q[k] += 0.5*sk*h[k]
            poly.append(q)
        dist = [float(np.dot(q - origin, nrm)) for q in poly]
        if min(dist) >= 0:
            continue
        clip = []
        for a in range(4):
            q0, q1, d0, d1 = poly[a], poly[(a + 1) % 4], dist[a], dist[(a + 1) % 4]
            if d0 >= 0:
                clip.append(q0)
            if (d0 >= 0) != (d1 >= 0):
                clip.append(q0 + d0/(d0 - d1)*(q1 - q0))
        if len(clip) < 3:
            th[m] = 0.0
            continue
        # area-weighted centroid of the fan triangles
        A, X = 0.0, np.zeros(3)
        for a in range(1, len(clip) - 1):
            t = 0.5*np.linalg.norm(np.cross(clip[a] - clip[0], clip[a + 1] - clip[0]))
            A += t
            X += t*(clip[0] + clip[a] + clip[a + 1])/3
        th[m] = A/(h[j]*h[k])
        if A > 0:
            xw[m] = X/A
    return xw, th


def boundary_values(d, kind, tan, f, N, bf, origin, nrm):
    """0/U and 0/p with the exact U on the outer faces (face centres from writeCellCentres).
    The solver's boundary flux is theta U.S (theta: fluid fraction of the face, ibmTheta). Point
    values at the face centres leave an O(h^2) imbalance; the outflow faces are scaled by
    1 + O(h^2) so that it closes (all faces are Dirichlet). Faces entirely in the solid: U = 0."""
    lo, hi, n, _, _, _ = config(kind, tan, f, N, bf)
    _, pc = read_field(os.path.join(d, "0", "C"))
    th = read_field(os.path.join(d, "0", "ibmTheta"))[1] if not bf else {}
    pu, pp = [], []
    for name, xf in pc.items():
        if name in ("frontAndBack",) or xf.ndim != 2:
            continue
        if name == "wall":
            pu.append("    wall { type noSlip; }")
            pp.append("    wall { type zeroGradient; }")
            continue
        # a cut face carries theta U_b.S: U_b at the centroid of its wet part (at the geometric
        # centre, on or behind the wall, the face would carry no flux or the wrong one)
        xw, _ = wet_centres(xf, lo, hi, n, origin, nrm) if not bf else (xf, None)
        U, _ = hz.exact(xw, A, NU, origin, nrm, kind)
        # faces in the solid get the continuation of the exact solution (divergence-free)
        theta = np.broadcast_to(th.get(name, 1.0), (len(U),)).astype(float)
        U[theta <= 0] = 0.0
        phi = theta*np.einsum("ij,ij->i", U, face_areas(xf, lo, hi, n, kind))
        out = phi > 0
        U[out] *= 1.0 - phi.sum()/phi[out].sum()
        pu.append("    %s { type fixedValue; value nonuniform List<vector> %d (%s); }"
                  % (name, len(U), " ".join("(%r %r %r)" % tuple(map(float, u)) for u in U)))
        pp.append("    %s { type zeroGradient; }" % name)
    empty = "    frontAndBack { type empty; }" if kind == "hiemenz" else ""
    write_dict(d / "0" / "U", "volVectorField", "U", "dimensions [0 1 -1 0 0 0 0];\ninternalField uniform (0 0 0);\n"
               "boundaryField\n{\n%s\n%s\n}\n" % ("\n".join(pu), empty))
    write_dict(d / "0" / "p", "volScalarField", "p", "dimensions [0 2 -2 0 0 0 0];\ninternalField uniform 0;\n"
               "boundaryField\n{\n%s\n%s\n}\n" % ("\n".join(pp), empty))


def wall_polygon(lo, hi, origin, nrm):
    """vertices (ordered) of the wall plane (x - origin).nrm = 0 inside the box [lo, hi]"""
    c = [np.array([[lo, hi][i][0], [lo, hi][j][1], [lo, hi][k][2]]) for i in (0, 1) for j in (0, 1) for k in (0, 1)]
    pts = []
    for a in range(8):
        for b in range(a + 1, 8):
            if bin(a ^ b).count("1") != 1:      # box edges join corners differing in one coordinate
                continue
            da, db = np.dot(c[a] - origin, nrm), np.dot(c[b] - origin, nrm)
            if da*db < 0:
                pts.append(c[a] + da/(da - db)*(c[b] - c[a]))
            elif da == 0:
                pts.append(c[a])
    P = np.unique(np.round(np.array(pts), 14), axis=0)
    m = P.mean(0)
    e1 = (P[0] - m)/np.linalg.norm(P[0] - m)
    e2 = np.cross(nrm, e1)
    return P[np.argsort(np.arctan2((P - m) @ e2, (P - m) @ e1))]


def wall_force_exact(kind, lo, hi, origin, nrm):
    """exact (pressure, viscous) force of the fluid on the wall inside the box, p0 = 0: the wall
    polygon fan-triangulated, mid-edge rule (exact for the quadratic wall pressure and the
    linear shear; a staircase quadrature over the box was off by 0.8 % for Homann)"""
    _, fpp0, _ = hz.solve(kind)
    tau = NU*A*math.sqrt(A/NU)*fpp0
    t = np.array([nrm[1], -nrm[0], 0.0])
    P = wall_polygon(lo, hi, origin, nrm)
    Fp, Fv = np.zeros(3), np.zeros(3)
    for i in range(1, len(P) - 1):
        a, b, c = P[0], P[i], P[i + 1]
        At = 0.5*np.dot(np.cross(b - a, c - a), nrm)
        for q in ((a + b)/2, (b + c)/2, (c + a)/2):
            r = q - origin
            # wall coordinate: s t in the Hiemenz plane (uniform in z), r - (r.n) n for Homann
            r = np.dot(r, t)*t if kind == "hiemenz" else r - np.dot(r, nrm)*nrm
            Fp += At/3*0.5*A*A*np.dot(r, r)*nrm        # -p n, p = -a^2 |r|^2/2 (n out of the body)
            Fv += At/3*tau*r
    return Fp, Fv


def one(ctx, kind, tan, f, N, bf=False, ma=0.5, conv=True):
    name = "%s_t%s_f%g_N%d%s_ma%g_c%d" % (kind, {0: "0", 1/3: "1-3", 1: "1"}[tan], f, N, "_bf" if bf else "", ma, conv)
    d = ctx.root / name
    origin, nrm = write_case(d, kind, tan, f, N, bf, ma, conv)
    common.sh("blockMesh", d, "log.blockMesh")
    (d / "0").mkdir(exist_ok=True)
    common.sh("postProcess -func writeCellCentres -time 0", d, "log.postProcess")
    common.sh("ibmSetGeometry -writeFaceCentres", d, "log.ibmSetGeometry")
    boundary_values(d, kind, tan, f, N, bf, origin, nrm)
    common.sh("ibmSimpleFoam", d, "log.ibmSimpleFoam")
    log = (d / "log.ibmSimpleFoam").read_text()
    t = latest_time(str(d))
    C, _ = read_field(os.path.join(d, "0", "C"))
    U, _ = read_field(os.path.join(d, t, "U"))
    p, _ = read_field(os.path.join(d, t, "p"))
    if bf:
        a = np.ones(len(C))
    else:
        a, _ = read_field(os.path.join(d, "0", "ibmAlpha"))
        # the cell value lives at its value point (the group centroid of a linked cell)
    xc = common.value_points(d, t, C) if not bf else C
    fl = a > 0
    Ue, pe = hz.exact(xc, A, NU, origin, nrm, kind)
    eu = np.linalg.norm(U - Ue, axis=1)
    ep = p - pe
    pmean = np.sum(a[fl]*ep[fl])/np.sum(a[fl])
    ep -= pmean
    lo, hi, _, _, _, ax = config(kind, tan, f, N, bf)
    # interior: farther than 0.1 from every outer face but a body-fitted wall (where the wall
    # meets the Dirichlet boundary the pressure is not converged, finding 2: reported apart).
    # The bottom counts too: the 45-degree wall meets it at x = -0.465, and its junction was
    # the "interior" Linf stall of finding 3 (2.8e-2, 6.9e-3, 2.1e-3, 5.1e-4 without it).
    dist = np.min([xc[:, i] - lo[i] for i in range(ax)] + [hi[i] - xc[:, i] for i in range(ax + 1)]
                  + ([xc[:, 1] - lo[1]] if kind == "homann" else [])
                  + ([xc[:, ax] - lo[ax]] if not bf else []), axis=0)
    inner = fl & (dist > 0.1)

    def norms(m):
        w = a[m]/a[m].sum()
        return (float(np.sqrt(np.sum(w*eu[m]**2))), float(eu[m].max()),
                float(np.sqrt(np.sum(w*ep[m]**2))), float(np.abs(ep[m]).max()))
    row = {"name": name, "kind": kind, "tan": tan, "f": f, "N": N, "bf": bf, "ma": ma, "conv": conv,
           "converged": "SIMPLE solution converged" in log}
    row.update(zip(("uL2", "uLinf", "pL2", "pLinf"), norms(fl)))
    row.update(zip(("uL2i", "uLinfi", "pL2i", "pLinfi"), norms(inner)))
    if not bf:
        Fp, Fv = wall_force_exact(kind, lo, hi, origin, nrm)
        F = common.ibm_force(log, "wall")
        # pressure: defined up to a constant (fixed in one reference cell); with the mean error
        # removed as in the norms (p - c on the wall pieces moves the force by -c sum Sw)
        Sw, _ = read_field(os.path.join(d, "0", "ibmSw"))
        Fpc = F["pressure"] - pmean*Sw.sum(0)
        row["eFp"] = float(np.linalg.norm(Fpc - Fp)/np.linalg.norm(Fp))
        # viscous: the tangential part (the wall shear); the normal part nu u_n/d of the discrete
        # wall flux is first order (exact: zero normal viscous stress on a no-slip wall), the
        # same on a body-fitted wall patch (doc §12 finding 2)
        ev = F["viscous"] - Fv
        row["eFvn"] = float(np.dot(ev, nrm)/max(np.linalg.norm(Fv), 1e-300))
        row["eFv"] = float(np.linalg.norm(ev - np.dot(ev, nrm)*nrm)/max(np.linalg.norm(Fv), 1e-300))
        row["Fv"], row["Fv_exact"] = F["viscous"].tolist(), Fv.tolist()
        row["Fp"], row["Fp_exact"] = F["pressure"].tolist(), Fp.tolist()
    return row


def run(ctx, res):
    N2 = [20, 40] if ctx.quick else [20, 40, 80, 160]
    N3 = [10, 20] if ctx.quick else [20, 40, 80]
    jobs = [("hiemenz", 0, f, N) for f in (0.1, 0.5, 0.9) for N in N2]
    jobs += [("hiemenz", t, 0.5, N) for t in (1/3, 1) for N in N2]
    jobs += [("hiemenz", 0, 0.5, N, True) for N in N2]
    jobs += [("hiemenz", 0, 0.1, N, False, ma, c) for N in N2 for ma, c in ((0, True), (0.5, False), (0, False))]
    jobs += [("homann", t, 0.5, N) for t in (0, 1/3) for N in N3]
    rows = common.pmap(lambda j: one(ctx, *j), jobs, ctx.jobs)
    res.data["rows"] = rows
    res.data["metrics"] = {r["name"] + "_" + k: r[k] for r in rows for k in ("uL2", "pL2i", "eFv") if k in r}
    res.le("all runs converged (SIMPLE residualControl)", sum(not r["converged"] for r in rows), 0)
    groups = sorted({(r["kind"], r["tan"], r["f"], r["bf"], r["ma"], r["conv"]) for r in rows},
                    key=lambda g: (g[0], g[3], g[1], g[2], -g[4], -g[5]))
    for g in groups:
        kind, tan, f, bf, ma, conv = g
        s = sorted([r for r in rows if (r["kind"], r["tan"], r["f"], r["bf"], r["ma"], r["conv"]) == g],
                   key=lambda r: r["N"])
        lab = "%s tan %.3g%s%s" % (kind, tan, " f %g" % f if tan == 0 and not bf else "",
                                    " body-fitted" if bf else "")
        opt = "" if (ma, conv) == (0.5, True) else " (mergeAlpha %g, ibmConvection %s)" % (ma, "on" if conv else "off")
        h = [1/r["N"] for r in s]
        default = (ma, conv) == (0.5, True)
        for k, what, tol in (("uL2", "u L2 order (all fluid cells)", 1.8), ("pL2i", "p L2 order (interior)", 1.8),
                             ("pL2", "p L2 order (all, incl. the wall/boundary junction)", None),
                             ("uLinfi", "u Linf order (interior)", None), ("pLinfi", "p Linf order (interior)", None)):
            p = common.order(h, [r[k] for r in s])
            note = " ".join("%.2e" % r[k] for r in s)
            if default and tol is not None:
                res.ge("%s, %s" % (what, lab), p, tol, note)
            else:
                res.report("%s, %s%s" % (what, lab, opt), p, note=note)
        if not bf and default:
            res.report("p Linf at the junction of the wall with the outer boundary, " + lab, s[-1]["pLinf"],
                       note=" ".join("%.2e" % r["pLinf"] for r in s))
        if not bf:
            pv = common.order(h, [r["eFv"] for r in s])
            nv = "rel. error " + " ".join("%.2e" % r["eFv"] for r in s)
            if default:
                res.ge("wall shear force (tangential) order, " + lab, pv, 1.5, nv)
            else:
                res.report("wall shear force (tangential) order, " + lab + opt, pv, note=nv)
            res.report("wall viscous force normal part / |F_v| (finest), " + lab + opt, s[-1]["eFvn"],
                       note=" ".join("%+.2e" % r["eFvn"] for r in s))
            res.report("wall pressure force rel. error (finest), " + lab + opt, s[-1]["eFp"],
                       note=" ".join("%.2e" % r["eFp"] for r in s))
