"""Cartesian background meshes with static refinement in nested boxes (topoSet + refineMesh,
hanging nodes), boundary data from patch face centres, sampling helpers."""
import os
from pathlib import Path

import numpy as np

import common
import geometry
from common import write_dict
from foamio import read_field


def refined(d, lo, hi, h0, boxes, patches, two_d=True, thickness=None):
    """blockMesh of lo..hi with cells of size h0, then one level of refinement (by 2 in x, y and,
    in 3-D, z) per box of `boxes` (list of (lo, hi), applied in order: nested boxes give nested
    levels). Returns the number of cells."""
    d = Path(d)
    lo, hi = np.array(lo, float), np.array(hi, float)
    n = [int(round((hi[i] - lo[i])/h0)) for i in range(3)]
    if two_d:
        hi[2] = lo[2] + (thickness or h0)
        n[2] = 1
    assert all(abs((hi[i] - lo[i])/n[i] - h0) < 1e-9*h0 for i in range(2)), "domain not a multiple of h0"
    geometry.box(d, lo, hi, n, patches)
    if not (d / "system" / "controlDict").exists():
        geometry.control_dict(d, "blockMesh")
    for f in ("fvSchemes", "fvSolution"):
        if not (d / "system" / f).exists():
            (d / "system" / f).write_text((common.TUTORIALS / "cylinderRe1" / "system" / f).read_text())
    common.sh("blockMesh", d, "log.blockMesh")
    dirs = "(tan1 tan2)" if two_d else "(tan1 tan2 normal)"
    write_dict(d / "system" / "refineMeshDict", "dictionary", "refineMeshDict", """
set refine; coordinateSystem global;
globalCoeffs { tan1 (1 0 0); tan2 (0 1 0); }
directions %s; useHexTopology true; geometricCut false; writeMesh false;
""" % dirs)
    for i, (blo, bhi) in enumerate(boxes):
        blo, bhi = np.array(blo, float), np.array(bhi, float)
        if two_d:
            blo[2], bhi[2] = lo[2] - 1, hi[2] + 1
        write_dict(d / "system" / "topoSetDict", "dictionary", "topoSetDict", """
actions ( { name refine; type cellSet; action new; source boxToCell; box %s %s; } );
""" % (common.foam_value(blo), common.foam_value(bhi)))
        common.sh("topoSet", d, "log.topoSet.%d" % i)
        common.sh("refineMesh -overwrite -dict system/refineMeshDict", d, "log.refineMesh.%d" % i)
    common.sh("checkMesh -constant", d, "log.checkMesh")
    return common.last_match((d / "log.checkMesh").read_text(), r"cells:\s+(\d+)", cast=int)


def face_centres(d):
    """{patch: (n, 3) face centres} from writeCellCentres at time 0"""
    (Path(d) / "0").mkdir(exist_ok=True)
    common.sh("postProcess -func writeCellCentres -time 0", d, "log.writeCellCentres")
    _, pc = read_field(os.path.join(d, "0", "C"))
    return pc


def vector_list(U):
    return "nonuniform List<vector> %d (%s)" % (len(U), " ".join("(%r %r %r)" % tuple(map(float, u)) for u in U))


def write_fields(d, fields, two_d=True):
    """fields: name -> (class, dimensions, internal, {patch: entries})"""
    for name, (cls, dims, internal, bcs) in fields.items():
        bc = "\n".join("    %s { %s }" % (p, e) for p, e in bcs.items())
        if two_d:
            bc += "\n    frontAndBack { type empty; }"
        write_dict(Path(d) / "0" / name, cls, name, "dimensions %s;\ninternalField %s;\nboundaryField\n{\n%s\n}\n"
                   % (dims, internal, bc))


def fluid_cells(d, t):
    """centroids, alpha, cell volume weights of the fluid cells (geometry at time 0 or t)"""
    C, _ = read_field(os.path.join(d, "0", "C"))
    geo = t if os.path.exists(os.path.join(d, t, "ibmAlpha")) else "0"
    a, _ = read_field(os.path.join(d, geo, "ibmAlpha"))
    return common.value_points(d, t, C), a


def point_value(xc, a, f, x0, h, radius=2.5, quadratic=False):
    """value at x0 (e.g. a wall point) by a least-squares fit (linear or quadratic, in x and y) over
    the fluid cells within radius*h, weighted by alpha"""
    dx = xc[:, :2] - np.asarray(x0, float)[:2]
    m = (a > 0) & (np.hypot(dx[:, 0], dx[:, 1]) < radius*h)
    X = dx[m]/h
    cols = [np.ones(len(X)), X[:, 0], X[:, 1]]
    if quadratic:
        cols += [X[:, 0]**2, X[:, 0]*X[:, 1], X[:, 1]**2]
    A = np.column_stack(cols)*np.sqrt(a[m])[:, None]
    c, *_ = np.linalg.lstsq(A, f[m]*np.sqrt(a[m]), rcond=None)
    return float(c[0])


def nested_boxes(fine_lo, fine_hi, h0, n, dom_lo, dom_hi, margin=12, wake=0.0):
    """boxes of refinement levels 1..n (coarse to fine) around the finest box: each coarser box is
    the next finer one grown by `margin` cells of its own level (x-max by (1 + wake) times that, for
    a wake), clipped to the domain"""
    boxes = [(np.array(fine_lo, float), np.array(fine_hi, float))]
    for i in range(n - 1, 0, -1):
        lo, hi = boxes[0]
        w = margin*h0/2**i
        lo, hi = lo - w, hi + w
        hi[0] += wake*w
        boxes.insert(0, (np.maximum(lo, dom_lo), np.minimum(hi, dom_hi)))
    out = []
    for lo, hi in boxes:   # snap to the mesh lines of the level being refined (h0/2^(i-1))
        out.append((lo, hi))
    return [(tuple(lo), tuple(hi)) for lo, hi in out]


def force_history(d, depth):
    """t, total force per unit depth (n, 3) from postProcessing/ibmForces (all restarts)"""
    rows = []
    for f in sorted(Path(d).glob("postProcessing/ibmForces/*/forces.dat"), key=lambda p: float(p.parent.name)):
        rows.append(np.loadtxt(f, ndmin=2))
    h = np.vstack(rows)
    return h[:, 0], (h[:, 1:4] + h[:, 4:7])/depth


def periodic_statistics(t, cl, cd, window):
    """means and r.m.s. of the fluctuations over the last `window` time units, Strouhal frequency
    (per time unit) from the up-crossings of C_L through its mean (nan if not oscillating)"""
    s = t >= t[-1] - window
    t, cl, cd = t[s], cl[s], cd[s]
    m = {"C_L": float(cl.mean()), "C_D": float(cd.mean()),
         "C_Lrms": float(np.sqrt(np.mean((cl - cl.mean())**2))), "C_Drms": float(np.sqrt(np.mean((cd - cd.mean())**2)))}
    x = cl - cl.mean()
    i = np.where((x[:-1] < 0) & (x[1:] >= 0))[0]
    tc = t[i] - x[i]*(t[i + 1] - t[i])/(x[i + 1] - x[i])
    m["f"] = float((len(tc) - 1)/(tc[-1] - tc[0])) if len(tc) >= 3 and m["C_Lrms"] > 1e-3 else float("nan")
    return m
