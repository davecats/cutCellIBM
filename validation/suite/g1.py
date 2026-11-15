"""G1 geometry: fluid volume and wall area of a cylinder (2-D, outside and inside) and a sphere
(3-D, outside and inside) against the exact values, closure of the wall segments, uniform field.
Also a 3-D cylinder crossing the walls z = 0 and z = L (axis z as in V14, and tilted so that its
trace on the walls is an ellipse): cells with both a wall patch face and a cut face, where the
closure sum_f theta_f S_f + S_w = 0 includes the wet parts of boundary faces (doc §11, V14)."""
import math

import numpy as np

import common
import geometry

ID, KIND, TIER = "G1", "gate", 0
TITLE = "geometry: volume, wall area, closure (cylinder 2-D, sphere 3-D, cylinder through walls 3-D, inverted)"
STRESSES, RELIES = ["K1"], []
REGRESSION = (1e-9,)

R = 1.0
C2 = (1.5123, 1.4871, 0.0)    # off the mesh lines
C3 = (1.5123, 1.4871, 1.5033)
# cylinder through the walls z = 0, L: axis z, and tilted by tan = |(0.1, 0.05)| (centre line from
# (1.3623, 1.4121) at z = 0 to (1.6623, 1.5621) at z = 3, the ellipse inside the domain)
AXES = {"cylz": (0.0, 0.0, 1.0), "cyltilt": (0.1, 0.05, 1.0)}
SCHEMES = common.TUTORIALS / "order" / "laplaceDisk" / "system"


def key(dim, N, invert, ns, shape=None):
    return "%dd_N%d_inv%d_ns%d" % (dim, N, int(invert), ns) + ("_" + shape if shape else "")


def make(ctx, dim, N, invert, ns, shape=None):
    d = ctx.root / key(dim, N, invert, ns, shape)
    L = 3.0
    if shape:
        geometry.box(d, (0, 0, 0), (L, L, L), (N, N, N))
        body = {"type": "cylinder", "centre": C3, "axis": AXES[shape], "radius": R}
    elif dim == 2:
        geometry.box(d, (0, 0, 0), (L, L, L/N), (N, N, 1))
        body = {"type": "cylinder", "centre": C2, "axis": (0, 0, 1), "radius": R}
    else:
        geometry.box(d, (0, 0, 0), (L, L, L), (N, N, N))
        body = {"type": "sphere", "centre": C3, "radius": R}
    body["invert"] = invert
    geometry.ibm_properties(d, {"body": body}, nSubSamples=ns, writeFields=False, writeForces=False)
    geometry.control_dict(d, "ibmSetGeometry")
    for f in ("fvSchemes", "fvSolution"):
        (d / "system" / f).write_text((SCHEMES / f).read_text())
    common.sh("blockMesh", d, "log.blockMesh")
    common.sh("ibmSetGeometry", d, "log.ibmSetGeometry")
    t = (d / "log.ibmSetGeometry").read_text()
    V = common.last_match(t, r"fluid volume ([-+0-9.eE]+)")
    A = common.last_match(t, r"total wall area ([-+0-9.eE]+)")
    net = np.linalg.norm(common.vec(common.last_match(t, r"net wall area vector (\([^)]*\))", cast=str)))
    g = common.last_match(t, r"uniform field, max magnitude: ([-+0-9.eE]+)")
    nabs = common.last_match(t, r"absorbed (\d+) sliver", cast=int)
    dz = L/N if dim == 2 else 1.0
    if shape:
        # slab 0 <= z <= L: every z-section is an ellipse of area pi R^2/cos, every generator has
        # length L/cos (cos = a_z/|a|): volume pi R^2 L/cos, lateral area 2 pi R L/cos
        cos = AXES[shape][2]/np.linalg.norm(AXES[shape])
        Vb, Ab = math.pi*R*R*L/cos, 2*math.pi*R*L/cos
        Vd = L**3
    elif dim == 2:
        Vb, Ab = math.pi*R*R*dz, 2*math.pi*R*dz
        Vd = L*L*dz
    else:
        Vb, Ab = 4/3*math.pi*R**3, 4*math.pi*R*R
        Vd = L**3
    Vex = Vb if invert else Vd - Vb
    return {"dim": dim, "N": N, "invert": invert, "ns": ns, "shape": shape, "h": L/N,
            "eV": abs(V - Vex)/Vb, "eA": abs(A - Ab)/Ab, "net": net/A, "grad": g, "absorbed": nabs}


def run(ctx, res):
    N2 = [20, 40] if ctx.quick else [20, 40, 80, 160]
    N3 = [16, 32] if ctx.quick else [16, 32, 64]
    jobs = [(2, N, inv, 32) for inv in (False, True) for N in N2]
    jobs += [(2, N, False, 8) for N in N2]
    jobs += [(3, N, inv, 16) for inv in (False, True) for N in N3]
    # nSubSamples 32 (the default, as in V14): the axis-z cylinder is the 2-D one stacked in z, its
    # quadrature error does not average over the surface as the sphere's does (16: 3.2e-4 at N = 16)
    jobs += [(3, N, inv, 32, s) for s in AXES for inv in (False, True) for N in N3]
    rows = common.pmap(lambda j: make(ctx, *j), jobs, ctx.jobs)
    res.data["rows"] = rows
    res.data["metrics"] = {key(r["dim"], r["N"], r["invert"], r["ns"], r.get("shape")) + "_" + k: r[k]
                           for r in rows for k in ("eV", "eA")}
    res.le("closure: max |net wall area|/area (closed bodies)",
           max(r["net"] for r in rows if not r.get("shape")), 1e-10)
    # open bodies: the net wall area is minus the sum of the wet boundary faces, i.e. the difference
    # of the two discrete (piecewise-linear) traces on z = 0, L, equal ellipses exactly: their
    # discretisation error (4.5e-5, 2.3e-5, 3.1e-7 at N = 16, 32, 64, the same with nSubSamples 16
    # and 32), not a closure error. The closure of each cell, boundary faces included, is the uniform-field check.
    res.report("net wall area/area, cylinders through walls (difference of the wall traces)",
               max(r["net"] for r in rows if r.get("shape")),
               note=" ".join("%.1e" % r["net"] for r in rows if r.get("shape")))
    res.le("closure: max |cut-cell gradient of a uniform field|", max(r["grad"] for r in rows), 1e-10)
    names = {None: None, "cylz": "cylinder through walls, axis z", "cyltilt": "cylinder through walls, tilted"}
    for dim, inv, ns, shape in [(2, False, 32, None), (2, True, 32, None), (3, False, 16, None),
                                (3, True, 16, None), (2, False, 8, None)] + \
            [(3, inv, 32, s) for s in AXES for inv in (False, True)]:
        s = [r for r in rows if (r["dim"], r["invert"], r["ns"], r.get("shape")) == (dim, inv, ns, shape)]
        lab = "%s %s, nSubSamples %d" % (names[shape] or ("cylinder" if dim == 2 else "sphere"),
                                        "inside" if inv else "outside", ns)
        pa = common.order([r["h"] for r in s], [r["eA"] for r in s])
        note = "rel. volume error %s; absorbed slivers %s" % (
            " ".join("%.2e" % r["eV"] for r in s), " ".join(str(r["absorbed"]) for r in s))
        if ns == 8:
            res.report("volume error (max rel.), " + lab, max(r["eV"] for r in s), note=note)
            res.report("wall area order, " + lab, pa,
                       note="rel. area error %s" % " ".join("%.2e" % r["eA"] for r in s))
        else:
            # sub-sampling quadrature, not a discretisation error: bounded, no order
            res.le("volume error (max rel. to the body), " + lab, max(r["eV"] for r in s), 1e-4, note)
            res.ge("wall area order, " + lab, pa, 1.8,
                   "rel. area error %s" % " ".join("%.2e" % r["eA"] for r in s))
