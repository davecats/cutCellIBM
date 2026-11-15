"""V2 periodic cylinder Re_D = 1 (tutorials/cylinderRe1), drag per unit depth on N = 41...321
against body-fitted O-grids (validation/data/secondOrder/freeSlip_zeroShear.txt): no slip, and
free slip, a stress-free wall (zeroShear, the default since 2026-11-15, doc/freeSlipPlan.md),
against the body-fitted zero-shear drag. The former IBM free slip (du/dn = 0, the OpenFOAM
`slip` condition) converged to the body-fitted slip drag, 24.8 % low."""
import common

ID, KIND, TIER = "V2", "case", 1
TITLE = "periodic cylinder Re_D = 1, drag, no slip and free slip"
STRESSES, RELIES = ["K3", "K4", "K8"], ["G1", "G2", "G3"]
REGRESSION = (1e-5,)

SRC = common.TUTORIALS / "cylinderRe1"
# body-fitted O-grid, M cells per quarter circle (runs/ibmSA_study/freeslip_bf)
BF = {"noSlip": {32: 5.85196, 64: 5.84752, 128: 5.84266, 256: 5.83958},
      "slip": {32: 3.10842, 64: 3.07644, 128: 3.05820, 256: 3.04843},
      "zeroShear": {32: 4.02338, 64: 4.01318, 128: 4.00675, 256: 4.00318}}
RTOL = 0.003   # the body-fitted M = 128 and 256 differ by 0.05 %


def one(ctx, wall, N, ma):
    d = ctx.root / ("%s_N%d_ma%g" % (wall, N, ma))
    common.copy_case(SRC, d)
    common.set_keys(d / "system" / "blockMeshDict", nx=N, ny=N)
    common.set_keys(d / "constant" / "ibmProperties", mergeAlpha=ma, writeFields=False)
    t = (d / "constant" / "ibmProperties").read_text()
    (d / "constant" / "ibmProperties").write_text(t.replace("wallType  noSlip;", "wallType  %s;" % wall))
    common.set_keys(d / "system" / "controlDict", endTime=40000, writeInterval=40000)
    for cmd in ("blockMesh", "ibmSetGeometry", "ibmSimpleFoam"):
        common.sh(cmd, d, "log." + cmd)
    log = (d / "log.ibmSimpleFoam").read_text()
    F = common.ibm_force(log, "cylinder")["total"][0]/(10.0/N)
    return {"wall": wall, "N": N, "ma": ma, "F": float(F), "converged": "SIMPLE solution converged" in log}


def run(ctx, res):
    Ns = [41, 81] if ctx.quick else [41, 81, 161, 321]
    jobs = [(w, N, 0.5) for w in ("noSlip", "freeSlip") for N in Ns] + [(w, N, 0) for w in ("noSlip",) for N in Ns]
    rows = common.pmap(lambda j: one(ctx, *j), jobs, ctx.jobs)
    res.data["rows"] = rows
    res.data["metrics"] = {"%s_N%d_ma%g" % (r["wall"], r["N"], r["ma"]): r["F"] for r in rows}
    res.le("runs not converged within 40000 iterations", sum(not r["converged"] for r in rows), 0)
    for wall, ma in (("noSlip", 0.5), ("noSlip", 0), ("freeSlip", 0.5)):
        s = sorted([r for r in rows if (r["wall"], r["ma"]) == (wall, ma)], key=lambda r: r["N"])
        Fe, p = common.richardson([10.0/r["N"] for r in s], [r["F"] for r in s])
        seq = " ".join("%.5f" % r["F"] for r in s)
        lab = "%s, mergeAlpha %g" % (wall, ma)
        if wall == "noSlip":
            res.rel("drag N = %d, %s, vs body-fitted M = 256" % (s[-1]["N"], lab), s[-1]["F"], BF["noSlip"][256],
                    RTOL, seq)
            res.report("drag Richardson extrapolation, " + lab, Fe, BF["noSlip"][256], "observed order %.2f" % p)
        else:
            res.rel("drag N = %d, %s, vs body-fitted zero shear M = 256 (stress-free wall)" % (s[-1]["N"], lab),
                    s[-1]["F"], BF["zeroShear"][256], RTOL, seq)
            res.report("drag Richardson extrapolation, " + lab, Fe, BF["zeroShear"][256], "observed order %.2f" % p)
