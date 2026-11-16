"""V11 BUB: stress-free sphere (clean spherical bubble of fixed shape) D = 1 in a uniform stream
U = 1, Re = 50 ... 300 (Blanco & Magnaudet 1995 table I, Moore 1963; references in
validation/reference_bubble.py, specification doc/validationSuite.md §10).

Mesh, domain and solver of V10 SPH (v10.setup: half domain with a symmetry plane, four nested
boxes, steady ibmSimpleFoam SIMPLEC), the sphere `freeSlip` with the zero-shear condition (the
default since 2190505). The flow is steady and axisymmetric at every Re (no separation on a clean
sphere), so D/h = 24, 32, 48 serve all Re = 50, 100, 200, 300; ranks 8, 16, 32 (68k cells each).

Quantities: C_D (judged within 2 % of table I, given to three digits); the tangential stress the
resolved field has at the surface, |(I - nn)(grad u + grad u^T) n| a/U (in units of mu U/a),
from an alpha-weighted quadratic least-squares fit of u over the fluid cells within 4 h of
points on the surface (theta 15-165 deg, four meridians); it must decrease with the mesh. The
same fit gives the surface vorticity against the identity omega = 2 u_theta/a of a stress-free
sphere (reported). Reversed flow on the axis behind the bubble is a FAIL (u_x sampled on
x = a + h ... a + 3). Once at D/h 32, Re 100: ibmConvection off, mergeAlpha 0, and zeroShear off
(the former free slip du_t/dn = 0, whose surface vorticity is u_theta/a)."""
import numpy as np

import common
import reference_bubble as ref
import v10
from foamio import latest_time

ID, KIND, TIER = "V11", "case", 2
TITLE = "BUB: stress-free sphere (clean bubble) in a uniform stream, Re 50-300"
STRESSES, RELIES = ["K4"], ["V10", "V2", "G7"]
REGRESSION = (1e-4,)

R = v10.R
DHS = [24, 32, 48]
RES = [50, 100, 200, 300]
NP = {24: 8, 32: 16, 48: 32}
THETA = np.radians(np.arange(15.0, 165.01, 5.0))      # from the front stagnation point
PHI = np.radians([0.0, 30.0, 60.0, 90.0])             # meridians in z >= 0 (y = cos, z = sin)


def wall_sample(F, h):
    """per surface point: theta (deg), tangential stress a/U |(I - nn) S n| (S = grad u + grad u^T),
    |u_t|, |omega| a / (2 |u_t|) (1 on a stress-free sphere, 1/2 for du_t/dn = 0)"""
    out = []
    for th in THETA:
        for ph in PHI:
            n = np.array([-np.cos(th), np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph)])
            v, g = F.fit(R*n, h)
            G = g[:3, :3]                                # G[i, j] = d u_i / d x_j
            t = (G + G.T) @ n
            tau = np.linalg.norm(t - (t @ n)*n)*R
            ut = np.linalg.norm(v[:3] - (v[:3] @ n)*n)
            om = np.array([G[2, 1] - G[1, 2], G[0, 2] - G[2, 0], G[1, 0] - G[0, 1]])
            out.append((np.degrees(th), tau, ut, np.linalg.norm(om)*R/(2*ut)))
    return np.array(out)


def axis(F, h):
    """min u_x on the axis behind the bubble and the recirculation length (v10.bubble)"""
    x = np.arange(R + h, R + 3.0, h/2)
    u = np.array([F.fit((xi, 0, 0), h)[0][0] for xi in x])
    return float(u.min()), float(x[u.argmin()] - R), v10.bubble(F, h)


def zs_info(d):
    """the solver's line on the zero-shear wall cells and the curvature it uses"""
    log = (d / "log.ibmSetGeometry").read_text() if (d / "log.ibmSetGeometry").exists() else ""
    log += (d / "log.ibmSimpleFoam").read_text()
    zl = [l.strip() for l in log.splitlines() if "zero-shear free slip:" in l]
    return zl[0] if zl else ""


def steady(ctx, dh, Re, tag="", **opts):
    d = ctx.root / ("re%d_dh%d%s" % (Re, dh, tag))
    ncells = v10.setup(d, dh, Re, wall="freeSlip", **opts)
    v10.decomposed(d, "ibmSimpleFoam", NP[dh])
    t, cd, cl, log = v10.forces(d)
    h = 1.0/dh
    F = v10.Field(d)
    W = wall_sample(F, h)
    umin, xmin, xs = axis(F, h)
    k = int(W[:, 1].argmax())
    eq = np.abs(W[:, 0] - 90.0) < 1e-6
    return {"Re": Re, "Dh": dh, "tag": tag, "cells": ncells, "C_D": float(cd[-1]), "C_L": float(cl[-1]),
            "iterations": int(t[-1]), "converged": "SIMPLE solution converged" in log,
            "C_Dchange": float(np.ptp(cd[-200:])/abs(cd[-1])),
            "tau_max": float(W[k, 1]), "tau_max_theta": float(W[k, 0]),
            "tau_rms": float(np.sqrt(np.mean(W[:, 1]**2))),
            "ut_max": float(W[:, 2].max()), "omega_ratio_eq": float(W[eq, 3].mean()),
            "omega_ratio_range": [float(W[:, 3].min()), float(W[:, 3].max())],
            "ux_axis_min": umin, "ux_axis_min_at": xmin, "x_s": xs, "zs": zs_info(d)}


def run(ctx, res):
    dhs, Res = (DHS[:1], RES[:2]) if ctx.quick else (DHS, RES)
    mid = DHS[1]
    jobs = [(dh, Re, {}) for dh in dhs for Re in Res]
    if not ctx.quick:   # once per case
        jobs += [(mid, 100, {"tag": "_noconv", "ibmConvection": False}),
                 (mid, 100, {"tag": "_merge0", "mergeAlpha": 0.0}),
                 (mid, 100, {"tag": "_zs0", "zeroShear": False})]
    import threading
    free, cv = [max(ctx.jobs, max(NP[j[0]] for j in jobs))], threading.Condition()

    def one(j):
        with cv:
            cv.wait_for(lambda: free[0] >= NP[j[0]])
            free[0] -= NP[j[0]]
        try:
            return steady(ctx, j[0], j[1], **j[2])
        finally:
            with cv:
                free[0] += NP[j[0]]
                cv.notify_all()
    jobs.sort(key=lambda j: (-j[0], j[1]))    # the finest (longest) first
    st = sorted(common.pmap(one, jobs, len(jobs)), key=lambda r: (r["tag"], r["Re"], r["Dh"]))
    res.data["steady"] = st
    base = [r for r in st if not r["tag"]]
    M = res.data["metrics"] = {}
    for r in st:
        for q in ("C_D", "tau_max", "ux_axis_min"):
            M["re%d_dh%d%s_%s" % (r["Re"], r["Dh"], r["tag"], q)] = r[q]
    res.le("steady runs not converged", sum(not r["converged"] for r in st), 0,
           " ".join("%d/%d:%.1e" % (r["Re"], r["Dh"], r["C_Dchange"]) for r in st if not r["converged"]))
    if base:
        res.report("zero-shear wall cells and curvature (D/h %d)" % base[-1]["Dh"], None, note=base[-1]["zs"])
    for Re in Res:
        rs = [r for r in base if r["Re"] == Re]
        if not rs:
            continue
        f, seq = rs[-1], lambda q, fm="%.5g": " ".join(fm % r[q] for r in rs)
        cref = ref.drag(Re)
        res.rel("Re %d C_D, D/h %d (Blanco & Magnaudet table I)" % (Re, f["Dh"]), f["C_D"], cref, 0.02,
                seq("C_D") + "; Moore %.4g" % float(ref.moore(Re)))
        if len(rs) >= 3:
            Ce, p = common.richardson([1/r["Dh"] for r in rs], [r["C_D"] for r in rs])
            res.report("Re %d observed order of C_D" % Re, p,
                       note="extrapolated %.5g (table I %.3g, %+.2f %%)" % (Ce, cref, 100*(Ce/cref - 1)))
        note = "%s; at theta %.0f deg; rms %s; max |u_t| %.3g" % (seq("tau_max", "%.3g"), f["tau_max_theta"],
                                                                 seq("tau_rms", "%.3g"), f["ut_max"])
        if len(rs) >= 2:
            a, b = rs[-2], rs[-1]
            ratio = b["tau_max"]/a["tau_max"]
            note += "; order %.2f" % (np.log(ratio)/np.log(a["Dh"]/b["Dh"]))
            res.le("Re %d max tangential stress (mu U/a) decreases with the mesh: D/h %d against %d"
                   % (Re, b["Dh"], a["Dh"]), ratio, 1.0, note)
        else:
            res.report("Re %d max tangential stress (mu U/a), D/h %d" % (Re, f["Dh"]), f["tau_max"], note=note)
        res.report("Re %d surface vorticity / (2 u_theta/a) at the equator, D/h %d" % (Re, f["Dh"]),
                   f["omega_ratio_eq"], 1.0, "%s; range over theta 15-165 deg %.3g .. %.3g"
                   % (seq("omega_ratio_eq", "%.4g"), *f["omega_ratio_range"]))
        res.le("Re %d recirculation length on the axis behind the bubble, D/h %d" % (Re, f["Dh"]), f["x_s"], 0.0,
               "min u_x on the axis %s (at x - a = %.3g)" % (seq("ux_axis_min", "%.4g"), f["ux_axis_min_at"]))
    ex = {r["tag"]: r for r in st if r["tag"]}
    for tag, what in (("_noconv", "ibmConvection off against on"), ("_merge0", "mergeAlpha 0 against 0.5"),
                      ("_zs0", "zeroShear off (du_t/dn = 0) against on")):
        if tag in ex:
            r0 = [r for r in base if r["Re"] == 100 and r["Dh"] == mid][0]
            e = ex[tag]
            res.report("Re 100 C_D, D/h %d: %s (rel.)" % (mid, what), e["C_D"]/r0["C_D"] - 1,
                       note="%.6g against %.6g; max tangential stress %.3g against %.3g; vorticity ratio at the "
                       "equator %.3g against %.3g; min u_x on the axis %.3g" % (e["C_D"], r0["C_D"], e["tau_max"],
                       r0["tau_max"], e["omega_ratio_eq"], r0["omega_ratio_eq"], e["ux_axis_min"]))
