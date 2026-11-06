"""G2 operator orders (ibmOrderTest, tutorials/order/laplaceDisk): Laplace and convection-
diffusion with Dirichlet and Neumann walls (LD LN CD CN), the cut-cell pressure gradient of a
linear (GL, exact against V_geo grad p) and a quadratic pressure (GQ), mergeAlpha 0 and 0.5;
GL also in 3-D inside and outside a sphere. secondOrder false reported only.

GL and GQ run unlinked only: ibmOrderTest compares grad(p) cell by cell, while a linked group
has one gradient fitted over the group (cutCellGeometry::lsqGrad), so the per-cell test is not
defined for linked cells (it reports every linked wall cell as "without a least-squares
gradient"). The linked pressure gradient is exercised in the flow by G3 and V1."""
import re

import common
import geometry

ID, KIND, TIER = "G2", "gate", 0
TITLE = "operator orders: LD LN CD CN GL GQ (ibmOrderTest), 2-D and 3-D GL"
STRESSES, RELIES = ["K1", "K2", "K3", "K5"], ["G1"]
REGRESSION = (0.01,)

SRC = common.TUTORIALS / "order" / "laplaceDisk"
PAT = r"orderTest (\w+) secondOrder=(\d) nCells=(\d+) iterations=(\d+) L2=(\S+) Linf=(\S+)"
GEO = r"against V_geo grad p\(x_c\): L2=(\S+) Linf=(\S+)"


def one(ctx, prob, N, ma, so, dim=2, invert=True):
    d = ctx.root / ("%s_%dd_N%d_ma%g_so%d_inv%d" % (prob, dim, N, ma, so, int(invert)))
    common.copy_case(SRC, d)
    if dim == 2:
        common.set_keys(d / "system" / "blockMeshDict", N=N)
    else:
        L = 3.013
        geometry.box(d, (0, 0, 0), (L, L, L), (N, N, N))
        geometry.ibm_properties(d, {"ball": {"type": "sphere", "centre": (1.5, 1.5, 1.5), "radius": 1,
                                             "invert": invert}},
                                mergeAlpha=0.5, nSubSamples=16, nBisect=40, writeFields=False,
                                writeForces=False, secondOrder=True)
    common.set_keys(d / "system" / "orderTestDict", problem=prob)
    common.set_keys(d / "constant" / "ibmProperties", secondOrder=bool(so), mergeAlpha=ma)
    common.sh("blockMesh", d, "log.blockMesh")
    common.sh("ibmOrderTest", d, "log.ibmOrderTest")
    t = (d / "log.ibmOrderTest").read_text()
    m = re.findall(PAT, t)[-1]
    g = re.findall(GEO, t)
    return {"prob": prob, "dim": dim, "N": N, "ma": ma, "so": so, "invert": invert,
            "L2": float(m[4]), "Linf": float(m[5]), "it": int(m[3]),
            "geoL2": float(g[-1][0]) if g else None}


def run(ctx, res):
    Ns = [20, 40] if ctx.quick else [20, 40, 80, 160]
    jobs = [(p, N, ma, 1) for p in ("LD", "LN", "CD", "CN") for N in Ns for ma in (0, 0.5)]
    jobs += [(p, N, 0, 1) for p in ("GL", "GQ") for N in Ns]
    jobs += [(p, N, 0, 0) for p in ("LD", "LN", "CD", "CN") for N in Ns]
    jobs3 = [("GL", N, 0, 1, 3, inv) for N in ([16] if ctx.quick else [16, 32]) for inv in (True, False)]
    rows = common.pmap(lambda j: one(ctx, *j), jobs + jobs3, ctx.jobs)
    res.data["rows"] = rows
    res.data["metrics"] = {"%s_%dd_N%d_ma%g_so%d_inv%d" % (r["prob"], r["dim"], r["N"], r["ma"], r["so"],
                                                          int(r["invert"])): r["L2"] for r in rows}
    for so in (1, 0):
        for prob in ("LD", "LN", "CD", "CN", "GQ"):
            for ma in ((0, 0.5) if so and prob != "GQ" else (0,)):
                s = sorted([r for r in rows if (r["prob"], r["ma"], r["so"], r["dim"]) == (prob, ma, so, 2)],
                           key=lambda r: r["N"])
                if not s:
                    continue
                # GQ is first order with geometric noise (pair orders 0.93 1.01 0.73 at the
                # baseline): fitted over all meshes
                last = 4 if prob == "GQ" else 3
                p = common.order([1/r["N"] for r in s], [r["L2"] for r in s], last)
                lab = "%s order (L2, finest %d), mergeAlpha %g%s" % (prob, last, ma, "" if so else ", secondOrder false")
                note = "L2 " + " ".join("%.3e" % r["L2"] for r in s)
                if not so:
                    res.report(lab, p, note=note)
                else:
                    res.ge(lab, p, 0.8 if prob == "GQ" else 1.8, note)
    gl = [r for r in rows if r["prob"] == "GL" and r["so"] == 1]
    res.le("GL exact: max L2 against V_geo grad p (2-D, unlinked)",
           max(r["geoL2"] for r in gl if r["dim"] == 2), 1e-10)
    res.le("GL exact: max L2 against V_geo grad p (3-D sphere, inside and outside, unlinked)",
           max(r["geoL2"] for r in gl if r["dim"] == 3), 1e-10)
