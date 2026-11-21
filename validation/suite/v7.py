"""V7 Spalart-Allmaras plane channel, Re_b = 2e4 (Re_tau ~ 545, tutorials/turbulence/channelSA):
immersed walls against the body-fitted channel with the same model (validation/channelSA_compare.py).
Tier 1: wall-resolved (y+ ~ 0.5), the walls cutting their cells at f = 0.05, 0.3, 0.5, linked
(default) and unlinked at f = 0.3.
Tier 2 (--tier2, doc/validationSuite.md §4 and §5): the cut-fraction sweep f = 0.01 ... 0.5 linked
and unlinked, wall-resolved and on uniform wall-function meshes (h+ ~ 10, 30, 90, against the
body-fitted channel with nutUSpaldingWallFunction on the same mesh and the wall-resolved
body-fitted solution); walls inclined to the grid at tan(alpha) = 1/8, 1/4, 1/2, 1 (uniform square
cells, three lattice offsets, initialised from the body-fitted profile, 200000 iterations judged by
stationarity); the options on the wall inclined at 1/4: ibmConvection off (x mergeAlpha 0.5, 0),
wallFunction none, ibmConvectionWallLaw on (wallFunction none also aligned at f = 0.3). The tier-2
runs are REPORT until their tolerances are decided."""
import math
import re
import subprocess

import numpy as np

import common
from channelSA_compare import load
from foamio import read_field, latest_time

ID, KIND, TIER = "V7", "case", 1
TITLE = "SA channel, immersed vs body-fitted walls (wall-resolved)"
STRESSES, RELIES = ["K12", "K13", "K10"], ["V1", "V2", "V5"]
REGRESSION = (1e-4,)

MAKE = common.TUTORIALS / "turbulence" / "channelSA" / "makeCase.py"
NU = 1e-4
END = 60000
CYCLE = 1000
END_INCLINED = 200000   # uniform fine cells develop slowly from the body-fitted profile

# tier 2
HPLUS = {10: 54, 30: 18, 90: 6}         # wall-function meshes: h+ (at u_tau ~ 0.0545) -> cells per half channel
FS = (0.01, 0.05, 0.15, 0.3, 0.5)       # wall-resolved cut fractions (the top wall sees 1 - f)
FS_WF = (0.05, 0.15, 0.3, 0.5)
TILTS = ((1, 8), (1, 4), (1, 2), (1, 1))
OFFSETS = (0.1, 0.4, 0.7)
OPTIONS = {"noConv": dict(ibmConvection=False), "noLaw": dict(law="none"), "convLaw": dict(ibmConvectionWallLaw=True)}
OPT_TILT, OPT_OFFSET = (1, 4), 0.4   # ibmConvection and ibmConvectionWallLaw act only on inclined walls
                                     # (aligned: no flow towards the wall, bit-identical results)


def mesh_args(mesh):
    """'wr': graded, y+ ~ 0.5; 'h10', 'h30', 'h90': uniform cells of that h+"""
    if mesh == "wr":
        return ["--mesh", "wr"]
    n = HPLUS[int(mesh[1:])]
    return ["--mesh", "wf", "--nh", str(n), "--h1", repr(1.0/n)]


def name_of(kind, mesh="wr", f=None, ma=0.5, opt=None, tilt=None):
    if kind == "bf":
        return "bf_" + mesh
    n = "ibm_%s_f%g_ma%g" % (mesh, f, ma) if tilt is None else "inc_%s_t%d-%d_f%g_ma%g" % (mesh, tilt[0], tilt[1], f, ma)
    return n + ("_" + opt if opt else "")


def one(ctx, kind, f=None, ma=0.5, mesh="wr", opt=None, tilt=None):
    o = dict(OPTIONS.get(opt, {}))
    name = name_of(kind, mesh, f, ma, opt, tilt)
    d = ctx.root / name
    if not (common.REUSE and d.exists()):
        args = [str(d), "--kind", kind] + mesh_args(mesh)
        if kind == "ibm":
            args += ["--offset", str(f), "--mergeAlpha", str(ma), "--wallFunction", o.pop("law", "spalartAllmaras")]
        if tilt:
            args += ["--incline", str(tilt[0]), str(tilt[1]), "--init", str(ctx.root / "bf_wr")]
        subprocess.run([common.sys.executable, str(MAKE)] + args, check=True, capture_output=True)
        common.set_keys(d / "system" / "controlDict", endTime=END_INCLINED if tilt else END)
        if tilt:   # 2-D: residualControl would stop while u_tau still drifts; judged by stationarity
            fs = d / "system" / "fvSolution"
            fs.write_text(re.sub(r"residualControl\s*\{[^}]*\}", "", fs.read_text()))
        if o:
            common.set_keys(d / "constant" / "ibmProperties", **o)
    app = re.search(r"application\s+(\w+)", (d / "system" / "controlDict").read_text()).group(1)
    if (d / "system" / "blockMeshDict").exists():
        common.sh("blockMesh", d, "log.blockMesh")
    (d / "0").mkdir(exist_ok=True)
    common.sh("postProcess -func writeCellCentres -time 0", d, "log.postProcess")
    if kind == "ibm":
        common.sh("ibmSetGeometry", d, "log.ibmSetGeometry")
    try:
        common.sh(app, d, "log." + app)
    except common.RunError as e:      # tier 2: a diverged run is a result, not a driver crash
        if ctx.tier2 and name not in TIER1_NAMES:
            return {"name": name, "kind": kind, "mesh": mesh, "f": f, "ma": ma, "opt": opt, "tilt": tilt,
                    "utau": float("nan"), "converged": False, "error": str(e)[:300]}
        raise
    log = (d / ("log." + app)).read_text()
    r = load(str(d)) if tilt is None else load_inclined(d, tilt)
    # residualControl cannot trigger here (Uy and p residuals are normalised by ~0): the
    # iteration must be stationary instead
    g, change = common.driving_gradient(log)
    out = {"name": name, "kind": kind, "f": f, "ma": ma, "utau": float(r["utau"]), "G": float(g[-1]),
           "change": change(), "converged": change() <= 1e-10}
    if ctx.tier2:
        # h+ ~ 10 has no steady state, body-fitted alike (finding 22): mean and half range of
        # u_tau over the last CYCLE iterations, reported without a tolerance (user, 2026-11-21)
        u = np.sqrt(np.abs(g[-CYCLE:]))
        out.update(mesh=mesh, opt=opt, tilt=tilt, utau_mean=float(np.sqrt(np.mean(g[-CYCLE:]))),
                   utau_amp=float(np.ptp(u)/2))
    for side, m in (("bottom", r["bottom"]), ("top", ~r["bottom"])):
        o_ = np.argsort(r["yp"][m])
        for q in (5, 30, 100, 300):
            out["%s_u+(%d)" % (side, q)] = float(np.interp(q, r["yp"][m][o_], r["up"][m][o_]))
        out["%s_yp1" % side] = float(r["yp"][m][o_][0])
    return out


def load_inclined(d, tilt):
    """u_tau from the driving gradient; u+ of the velocity along the channel against the distance
    of the value point to the nearer wall (the walls at eta = -1, 1 in the wall-normal coordinate
    of makeCase.py); bottom: the side below the centre line"""
    log = [p for p in d.glob("log.*Foam")][0].read_text()
    utau = math.sqrt(float(re.findall(r"pressure gradient = (\S+)", log)[-1]))
    a = math.atan2(tilt[0], tilt[1])
    t, n = np.array([math.cos(a), math.sin(a), 0.0]), np.array([-math.sin(a), math.cos(a), 0.0])
    tl = latest_time(str(d))
    C, _ = read_field(str(d / "0" / "C"))
    U, _ = read_field(str(d / tl / "U"))
    alpha, _ = read_field(str(d / "0" / "ibmAlpha"))
    eta = common.value_points(d, tl, C) @ n
    fl = alpha > 0
    return dict(utau=utau, yp=((1 - np.abs(eta))*utau/NU)[fl], up=(U @ t)[fl]/utau, bottom=eta[fl] < 0)


TIER1 = [("bf",), ("ibm", 0.05), ("ibm", 0.3), ("ibm", 0.5), ("ibm", 0.3, 0)]
TIER1_NAMES = {name_of(*j[:1], "wr", *j[1:]) for j in TIER1}


def tier2_jobs():
    """keyword dicts of the tier-2 runs beyond tier 1 (bf_wr first: the inclined runs start from it)"""
    bf = [dict(kind="bf", mesh=m) for m in ["h%d" % h for h in HPLUS]]
    jobs = []
    for mesh, fs in [("wr", FS)] + [("h%d" % h, FS_WF) for h in HPLUS]:
        for f in fs:
            for ma in (0.5, 0.0):
                jobs.append(dict(kind="ibm", mesh=mesh, f=f, ma=ma))
        jobs.append(dict(kind="ibm", mesh=mesh, f=0.3, opt="noLaw"))
    for mesh in ["wr"] + ["h%d" % h for h in HPLUS]:
        for tilt in TILTS:
            for f in OFFSETS:
                jobs += [dict(kind="ibm", mesh=mesh, f=f, ma=ma, tilt=tilt) for ma in ((0.5, 0.0) if mesh == "wr" else (0.5,))]
        # the options on an inclined wall; ibmConvection x mergeAlpha where unlinked runs converge
        jobs += [dict(kind="ibm", mesh=mesh, f=OPT_OFFSET, ma=ma, tilt=OPT_TILT, opt="noConv")
                 for ma in ((0.5, 0.0) if mesh == "wr" else (0.5,))]
        jobs += [dict(kind="ibm", mesh=mesh, f=OPT_OFFSET, tilt=OPT_TILT, opt=o) for o in ("noLaw", "convLaw")]
    jobs = [j for j in jobs if name_of(j["kind"], j.get("mesh", "wr"), j.get("f"), j.get("ma", 0.5),
                                       j.get("opt"), j.get("tilt")) not in TIER1_NAMES]
    return bf, jobs


def run(ctx, res):
    rows = common.pmap(lambda j: one(ctx, *j), TIER1, ctx.jobs)
    if ctx.tier2:
        bf, jobs = tier2_jobs()
        # longest first: the inclined wall-resolved runs
        jobs.sort(key=lambda j: (j.get("tilt") is None, j["mesh"] != "wr"))
        rows += common.pmap(lambda j: one(ctx, **j), bf + jobs, ctx.jobs)
    res.data["rows"] = rows
    res.data["metrics"] = {r["name"] + "_utau": value_of(r) for r in rows if np.isfinite(r["utau"])}
    t1 = [r for r in rows if r["name"] in TIER1_NAMES]
    res.le("runs not stationary (driving gradient over the last 100 iterations > 1e-10)", sum(not r["converged"] for r in t1), 0)
    bf = [r for r in rows if r["name"] == "bf_wr"][0]
    for r in t1:
        if r["kind"] == "bf":
            continue
        du = max(abs(r["%s_u+(%d)" % (s, q)] - bf["%s_u+(%d)" % (s, q)]) for s in ("bottom", "top")
                 for q in (5, 30, 100, 300))
        res.rel("u_tau, f %g, mergeAlpha %g, vs body-fitted" % (r["f"], r["ma"]), r["utau"], bf["utau"], 0.005,
                "max |du+| at y+ = 5, 30, 100, 300: %.3f" % du)
    if ctx.tier2:
        tier2_report(res, rows, bf)


def tier2_report(res, rows, bfwr):
    """REPORT per tier-2 run: u_tau against the wall-resolved body-fitted solution, and on the
    wall-function meshes against the body-fitted wall-function solution on the same mesh"""
    bfm = {r["mesh"]: r for r in rows if r["kind"] == "bf" and r.get("mesh")}
    t2 = [r for r in rows if r["name"] not in TIER1_NAMES and r["kind"] == "ibm"]
    bad = [r for r in t2 if not r["converged"] and r.get("mesh") != "h10"]
    res.report("tier 2: runs not stationary or failed (h+ 10 excluded: no steady state)", len(bad),
               note=", ".join(r["name"] for r in bad) or "-")
    for r in t2 + [r for m, r in sorted(bfm.items()) if m != "wr"]:
        res.report("tier 2 u_tau %s" % r["name"], value_of(r), bfwr["utau"], note_of(r, bfwr, bfm))


def value_of(r):
    """u_tau; on h+ ~ 10 the mean over the cycle"""
    return r.get("utau_mean", r["utau"]) if r.get("mesh") == "h10" else r["utau"]


def note_of(r, bfwr, bfm):
    note = "vs wall-resolved body-fitted %+.2f %%" % (100*(value_of(r)/bfwr["utau"] - 1))
    if r["kind"] == "ibm" and r["mesh"] != "wr" and r["mesh"] in bfm:
        note += "; body-fitted wall function on this mesh %+.2f %%" % (100*(value_of(bfm[r["mesh"]])/bfwr["utau"] - 1))
    if r.get("mesh") == "h10" and "utau_amp" in r and not r["converged"]:
        note += "; mean over the last %d iterations, cycle half range %.2f %% (no steady state at h+ 10, " \
                "body-fitted alike: no tolerance)" % (CYCLE, 100*r["utau_amp"]/value_of(r))
    elif not r["converged"]:
        note += "; NOT stationary" + (" (%s)" % r["error"][:80] if r.get("error") else "")
    return note
