"""Figures and LaTeX tables of one commit's results, for the consolidated report:
results/<commit>/figures/<ID>_<what>.{pdf,png} from the result JSONs and the case directories
they name (flow pictures of the 2-D cases: flowfield.py), results/<commit>/tables/<ID>.tex (every check: status, value, reference, tolerance,
note). The narrative (results/<commit>/report/report.tex) includes them.

    latexreport.py <commit> [--diag DIR ...]

Diagnostic runs outside the suite (V5 time step / outer correctors, V12 outer correctors) are
read from ~/runs/suite/diag_V5 and diag_V12 when present."""
import json
import sys
from pathlib import Path

import numpy as np

import common
import flowfield
import report
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

sys.path.insert(0, str(common.VALIDATION))
from foamio import latest_time, read_field   # noqa: E402

# categorical slots in fixed order (validated palette), markers as the second channel
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
M = ["o", "s", "^", "D", "v", "P", "X", "*"]
REF = "#3b3b3b"
STATUS = {"PASS": "passc", "FAIL": "failc", "REPORT": "reportc"}

plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "legend.fontsize": 7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "lines.linewidth": 1.4, "lines.markersize": 4,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e4e4e4", "grid.linewidth": 0.5, "legend.frameon": False,
                     "savefig.bbox": "tight", "savefig.dpi": 160})


class Out:
    def __init__(self, d):
        self.d = d
        (d / "figures").mkdir(exist_ok=True)
        (d / "tables").mkdir(exist_ok=True)

    def save(self, fig, name):
        for ax in fig.axes:   # log abscissae of mesh sequences: label the meshes themselves
            if ax.get_xscale() != "log":
                continue
            x = sorted({float(v) for ln in ax.get_lines() for v in np.atleast_1d(ln.get_xdata())
                        if np.isfinite(v) and v > 0})
            if 1 < len(x) <= 10:
                ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(x))
                ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: "%g" % v))
                ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        for ext in ("pdf", "png"):
            fig.savefig(self.d / "figures" / ("%s.%s" % (name, ext)))
        plt.close(fig)
        print("  figure", name)


def slope(ax, x0, y0, p, x1, label=None):
    """guide line of order p through (x0, y0) for h = 1/N plotted against N"""
    x = np.array([x0, x1])
    ax.plot(x, y0*(x/x0)**(-p), color=REF, lw=0.8, ls=":")
    ax.annotate(label or "order %g" % p, (x1, y0*(x1/x0)**(-p)), fontsize=6.5, color=REF,
                xytext=(0, -3), textcoords="offset points", va="top", ha="right")


def loglog(ax, xlabel="N", ylabel=None):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)


def series(ax, i, x, y, label, ls="-", **kw):
    ax.plot(x, y, ls, color=C[i], marker=M[i], label=label, **kw)


# ---------------------------------------------------------------------------------------- tables

def tex(s):
    s = str(s)
    for a, b in (("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"), ("#", r"\#"), ("_", r"\_"),
                 ("^", r"\^{}"), ("~", r"\~{}"), ("{", r"\{"), ("}", r"\}"), ("|", r"\textbar{}"),
                 ("<=", r"$\le$"), (">=", r"$\ge$"), ("±", r"$\pm$"), ("<", r"$<$"), (">", r"$>$")):
        s = s.replace(a, b)
    return s


def brk(s):
    """tex() with break points after underscores and slashes (long case and metric names)"""
    return tex(s).replace(r"\_", r"\_\allowbreak{}").replace("/", r"/\allowbreak{}")


# fixed column widths that add up to less than the text width (status box ~1.2 cm, 10 x 3 pt
# column separation): check, value, reference, tolerance, note
COLS = (r"@{}l p{0.36\textwidth} >{\raggedleft\arraybackslash}p{0.09\textwidth} "
        r">{\raggedleft\arraybackslash}p{0.09\textwidth} >{\raggedright\arraybackslash}p{0.10\textwidth} "
        r">{\raggedright\arraybackslash}p{0.20\textwidth}@{}")


def fmt(v):
    if v is None:
        return "--"
    if isinstance(v, float):
        return "%.4g" % v
    return str(v)


def table(out, r, name=None, keep=lambda c: True, notes=True):
    L = [r"\begin{longtable}{%s}" % COLS, r"\toprule",
         r"status & check & value & reference & tol. & note\\", r"\midrule", r"\endhead"]
    for c in r["checks"]:
        if not keep(c) or c["name"] == "baseline":
            continue
        note = (c["note"] or "") if notes else ""
        L.append(r"\st%s & %s & %s & %s & {\scriptsize %s} & {\scriptsize %s}\\" % (
            c["status"].lower(), brk(c["name"]), tex(fmt(c["value"])), tex(fmt(c["ref"])),
            brk(c["tol"] or "--"), brk(note)))
    L += [r"\bottomrule", r"\end{longtable}"]
    (out.d / "tables" / ("%s.tex" % (name or r["id"]))).write_text("\n".join(L) + "\n")


def summary_table(out, res):
    """status of the checks against exact and reference data, and separately the change against the
    stored baseline (a FAIL there is a change, intended or not, not a wrong result)"""
    L = [r"\begin{tabular}{@{}l p{0.36\textwidth} c c r r r r@{}}", r"\toprule",
         r"id & title & checks & vs.\ baseline & PASS & FAIL & REPORT & time [s]\\", r"\midrule"]
    for r in res:
        cs = [c for c in r["checks"] if c["name"] != "baseline"]
        n = {s: sum(c["status"] == s for c in cs) for s in STATUS}
        st = "\\st" + report.reference_status(r).lower()
        reg = [c for c in cs if report.is_regression(c)]
        rst = ("\\st" + report.baseline_status(r).lower()) if reg else "--"
        base = sorted({m for c in reg for m in __import__("re").findall(r"baseline (\w+)", c["name"] + " " + str(c.get("note") or ""))})
        L.append(r"\hyperref[sec:%s]{%s} & %s & %s & %s & %d & %d & %d & %.0f\\" % (
            r["id"], r["id"], tex(r["TITLE"]), st, rst + (r" {\scriptsize %s}" % ", ".join(base) if base else ""),
            n["PASS"], n["FAIL"], n["REPORT"], r["seconds"]))
    L += [r"\bottomrule", r"\end{tabular}"]
    (out.d / "tables" / "summary.tex").write_text("\n".join(L) + "\n")


# ---------------------------------------------------------------------------------------- gates

def g1(out, r):
    rows = r["data"]["rows"]
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.5))
    groups = [(2, False, 32, "2-D cylinder, outside, nSub 32"), (2, True, 32, "2-D cylinder, inside, nSub 32"),
              (2, False, 8, "2-D cylinder, outside, nSub 8"), (3, False, 16, "3-D sphere, outside, nSub 16"),
              (3, True, 16, "3-D sphere, inside, nSub 16")]
    for i, (dim, inv, ns, lab) in enumerate(groups):
        s = sorted([x for x in rows if (x["dim"], x["invert"], x["ns"]) == (dim, inv, ns)], key=lambda x: x["N"])
        if not s:
            continue
        N = [x["N"] for x in s]
        series(ax[0], i, N, [abs(x["eA"]) for x in s], lab, ls="--" if ns == 8 else "-")
        series(ax[1], i, N, [abs(x["eV"]) for x in s], lab, ls="--" if ns == 8 else "-")
    loglog(ax[0], ylabel="rel. wall-area error")
    loglog(ax[1], ylabel="rel. volume error")
    slope(ax[0], 20, 3e-3, 2, 160)
    ax[1].axhline(1e-4, color=REF, lw=0.8, ls="--")
    ax[1].annotate("tolerance 1e-4", (20, 1.15e-4), fontsize=6.5, color=REF)
    ax[0].set_title("wall area (order $\\geq$ 1.8)")
    ax[1].set_title("fluid volume (bounded)")
    ax[1].legend(loc="upper left", bbox_to_anchor=(1.02, 1))
    out.save(fig, "G1_geometry")


def g2(out, r):
    rows = [x for x in r["data"]["rows"] if x["dim"] == 2 and x["L2"] is not None]
    probs = ["LD", "LN", "CD", "CN", "GQ"]
    fig, ax = plt.subplots(1, 5, figsize=(9.5, 2.3), sharey=True)
    for a, p in zip(ax, probs):
        for i, (ma, so, lab) in enumerate(((0, 1, "unlinked"), (0.5, 1, "mergeAlpha 0.5"),
                                           (0, 0, "unlinked, secondOrder off"))):
            s = sorted([x for x in rows if x["prob"] == p and x["ma"] == ma and x["so"] == so], key=lambda x: x["N"])
            if s:
                series(a, i, [x["N"] for x in s], [x["L2"] for x in s], lab, ls="--" if so == 0 else "-")
        loglog(a)
        a.set_title(p)
        slope(a, 20, 1.2e-3 if p != "GQ" else 1e-2, 2 if p != "GQ" else 1, 160)
    ax[0].set_ylabel("L2 error")
    ax[0].legend(*ax[0].get_legend_handles_labels(), loc="upper left", bbox_to_anchor=(0, -0.25), ncol=3)
    out.save(fig, "G2_operators")


def g3(out, r):
    rows = r["data"]["rows"]
    fig, ax = plt.subplots(1, 2, figsize=(5.6, 2.4))
    for i, ma in enumerate((0.0, 0.5)):
        s = sorted([x for x in rows if float(x["ma"]) == ma], key=lambda x: x["N"])
        lab = "mergeAlpha %g" % ma
        series(ax[0], i, [x["N"] for x in s], [x["maxU"] for x in s], lab)
        series(ax[1], i, [x["N"] for x in s], [x["maxEp"] for x in s], lab)
    loglog(ax[0], ylabel="max |U| (spurious)")
    loglog(ax[1], ylabel="max |p + y| (pressure error)")
    slope(ax[0], 40, 1.5e-4, 2, 160)
    slope(ax[1], 40, 3e-2, 1, 160)
    ax[0].legend()
    out.save(fig, "G3_hydrostatic")


def g6(out, r):
    m = r["data"]["metrics"]
    k = list(m)
    v = [m[x] for x in k]
    tol = [1e-8 if x.startswith("steady") else 1e-9 for x in k]
    fig, ax = plt.subplots(figsize=(5.5, 2.6))
    y = np.arange(len(k))
    ax.barh(y, v, color=C[0], height=0.55)
    ax.scatter(tol, y, marker="|", s=120, color=REF, label="tolerance", zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([x.replace("_budget", "") for x in k])
    ax.set_xscale("log")
    ax.set_xlim(1e-13, 1e-7)
    ax.set_xlabel("budget error per step (relative)")
    ax.invert_yaxis()
    ax.legend(loc="lower right")
    out.save(fig, "G6_budgets")


def g8(out, r):
    F = r["data"]["forces"]
    lab, val = [], []
    for n in (4, 6):
        lab.append("cylinder drag, %d ranks" % n)
        val.append(abs(F["cylinder_np%d" % n][0]/F["cylinder_np1"][0] - 1))
        lab.append("channel $G$, %d ranks" % n)
        val.append(abs(F["channel_np%d" % n][0]/F["channel_np1"][0] - 1))
    fig, ax = plt.subplots(figsize=(5.0, 1.8))
    y = np.arange(len(lab))
    floor = 1e-16
    ax.barh(y, [max(v, floor) for v in val], color=C[0], height=0.55)
    for yi, v in zip(y, val):
        if v == 0:
            ax.annotate("identical (0)", (floor, yi), xytext=(4, 0), textcoords="offset points", va="center",
                        fontsize=6.5)
    ax.axvline(1e-8, color=REF, ls="--", lw=0.8)
    ax.annotate("tolerance", (1e-8, len(lab) - 0.6), fontsize=6.5, color=REF, xytext=(3, 0),
                textcoords="offset points")
    ax.set_yticks(y)
    ax.set_yticklabels(lab)
    ax.set_xscale("log")
    ax.set_xlim(floor/3, 1e-6)
    ax.set_xlabel("|parallel / serial - 1|")
    ax.invert_yaxis()
    out.save(fig, "G8_decomposition")


def g9(out, r):
    rows = r["data"]["rows"]
    base = dict(ma=0.5, conv=True)
    grp = [("Hiemenz level, f 0.1", dict(kind="hiemenz", tan=0, f=0.1, bf=False)),
           ("Hiemenz level, f 0.5", dict(kind="hiemenz", tan=0, f=0.5, bf=False)),
           ("Hiemenz level, f 0.9", dict(kind="hiemenz", tan=0, f=0.9, bf=False)),
           ("Hiemenz tan 1/3", dict(kind="hiemenz", tan=1/3, bf=False)),
           ("Hiemenz 45°", dict(kind="hiemenz", tan=1, bf=False)),
           ("Hiemenz body-fitted", dict(kind="hiemenz", bf=True)),
           ("Homann 3-D, level", dict(kind="homann", tan=0, bf=False)),
           ("Homann 3-D, tan 1/3", dict(kind="homann", tan=1/3, bf=False))]
    fig, ax = plt.subplots(1, 4, figsize=(10, 2.6))
    for i, (lab, kw) in enumerate(grp):
        s = [x for x in rows if all(abs(x[k] - v) < 1e-6 if isinstance(v, float) else x[k] == v
                                    for k, v in {**kw, **base}.items())]
        s.sort(key=lambda x: x["N"])
        if not s:
            continue
        N = [x["N"] for x in s]
        series(ax[0], i, N, [x["uL2"] for x in s], lab)
        series(ax[1], i, N, [x["pL2i"] for x in s], lab)
        series(ax[2], i, N, [x["pLinfi"] for x in s], lab)
        series(ax[3], i, N, [x["pLinf"] for x in s], lab)
    for a, t in zip(ax, ("u, L2 (all fluid cells)", "p, L2 (interior)", "p, max (interior)",
                         "p, max (all, incl. the junction)")):
        loglog(a)
        a.set_title(t)
        slope(a, 20, 1e-2, 2, 160)
    ax[0].set_ylabel("error")
    ax[3].legend(loc="upper left", bbox_to_anchor=(1.02, 1))
    out.save(fig, "G9_stagnation")


# ---------------------------------------------------------------------------------------- cases

def v1(out, r):
    rows = r["data"]["rows"]
    fig, ax = plt.subplots(1, 2, figsize=(6.2, 2.5))
    for i, (ma, so, lab) in enumerate(((0.5, 1, "linked (mergeAlpha 0.5)"), (0.0, 1, "unlinked"),
                                       (0.0, 0, "unlinked, secondOrder off"))):
        s = sorted([x for x in rows if float(x["ma"]) == ma and x["so"] == so], key=lambda x: x["N"])
        N = [x["N"] for x in s]
        series(ax[0], i, N, [x["uL2"] for x in s], lab, ls="--" if so == 0 else "-")
        series(ax[1], i, N, [x["pL2"] for x in s], lab, ls="--" if so == 0 else "-")
    loglog(ax[0], ylabel="u error, L2")
    loglog(ax[1], ylabel="p error, L2")
    slope(ax[0], 40, 8e-4, 2, 320)
    slope(ax[1], 40, 3e-3, 2, 320)
    slope(ax[1], 40, 1.5e-2, 1, 320)
    ax[0].set_title("velocity: second order")
    ax[1].set_title("pressure: below second order (open)")
    ax[1].legend(loc="upper left", bbox_to_anchor=(1.02, 1))
    out.save(fig, "V1_taylorCouette")


def v2(out, r):
    rows = r["data"]["rows"]
    fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.5))
    for i, (w, ma, lab) in enumerate((("noSlip", 0.5, "no slip, mergeAlpha 0.5"), ("noSlip", 0.0, "no slip, unlinked"))):
        s = sorted([x for x in rows if x["wall"] == w and float(x["ma"]) == ma], key=lambda x: x["N"])
        series(ax[0], i, [x["N"] for x in s], [x["F"] for x in s], lab)
    ax[0].axhline(5.83958, color=REF, ls="--", lw=0.9)
    ax[0].annotate("body-fitted M = 256", (41, 5.8405), fontsize=6.5, color=REF)
    s = sorted([x for x in rows if x["wall"] == "freeSlip"], key=lambda x: x["N"])
    series(ax[1], 2, [x["N"] for x in s], [x["F"] for x in s], "free slip, mergeAlpha 0.5")
    for v, t in ((3.04843, "body-fitted slip"), (4.00318, "body-fitted zero shear")):
        ax[1].axhline(v, color=REF, ls="--", lw=0.9)
        ax[1].annotate(t, (41, v + 0.03), fontsize=6.5, color=REF)
    for a in ax:
        a.set_xscale("log")
        a.set_xlabel("N")
        a.set_ylabel("drag")
        a.legend(loc="lower right")
    ax[0].set_title("no slip")
    ax[1].set_title("free slip (stress-free wall, zeroShear)")
    out.save(fig, "V2_cylinderRe1")


def v3(out, r):
    root = Path(r["cases"])

    def fx(d):
        dat = np.loadtxt(sorted(d.glob("postProcessing/ibmForces/*/forces.dat"))[-1])
        Cf, _ = read_field(str(d / "0" / "Cf"))
        return dat[:, 0], (dat[:, 1] + dat[:, 4])/(2*Cf[:, 2].max())

    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.5))
    h = np.load(common.VALIDATION / "ref_osc_fill.npz", allow_pickle=True)["history"]
    s = h[:, 0] > 0.3
    ax[0].plot(h[s, 0], h[s, 5], color=REF, lw=2.6, alpha=0.35, label="notebook (--fillUncovered)")
    for i, opt in enumerate(("notebook", "default")):
        t, F = fx(root / ("osc_" + opt))
        s = t > 0.3
        ax[0].plot(t[s], F[s], color=C[i], lw=0.9, label="OpenFOAM, %s options" % opt)
    ax[0].set_xlabel("t")
    ax[0].set_ylabel("$F_x$")
    ax[0].set_title("oscillating cylinder, $F_x$ (after the impulsive start)")
    ax[0].legend(loc="lower right")
    for i, (case, ref, lab) in enumerate((("osc_notebook", "ref_osc_fill.npz", "oscillating"),
                                          ("gal_notebook", "ref_gal_seamFix_fill.npz", "translating"))):
        h = np.load(common.VALIDATION / ref, allow_pickle=True)["history"]
        t, F = fx(root / case)
        n = min(len(t), len(h))
        e = (F[:n] - h[:n, 5])/np.abs(h[10:n, 5]).max()   # after the impulsive start (first 10 steps)
        ax[1].plot(t[10:n]/t[n - 1], e[10:], color=C[[0, 2][i]], lw=0.9,
                   label="%s (steps 1-2: %.1e of the peak)" % (lab, np.abs(F[:3] - h[:3, 5]).max()/np.abs(h[:n, 5]).max()))
    ax[1].set_xlabel("$t / t_{end}$")
    ax[1].set_ylabel("$(F_x - F_{x,nb})$ / max$|F_{x,nb}|$")
    ax[1].set_title("notebook options: difference after step 10 (finding 5 resolved)")
    ax[1].legend(loc="upper left", bbox_to_anchor=(1.02, 1))
    out.save(fig, "V3_moving")


def v4(out, r):
    rows = r["data"]["rows"]
    ib = sorted([x for x in rows if x["kind"] == "ibm"], key=lambda x: x["N"])
    bf = sorted([x for x in rows if x["kind"] == "bf"], key=lambda x: x["M"])
    fig, ax = plt.subplots(figsize=(3.4, 2.4))
    series(ax, 0, [x["N"] for x in ib], [x["ratio"] for x in ib], "immersed wall")
    series(ax, 1, [x["M"] for x in bf], [x["ratio"] for x in bf], "body-fitted (cells per side)")
    ax.axhspan(bf[-1]["ratio"]*0.99, bf[-1]["ratio"]*1.01, color="#e8e8e8", lw=0, label="±1 % of body-fitted M = 128")
    ax.set_xscale("log")
    ax.set_xlabel("N, M")
    ax.set_ylabel("scalar content / source rate")
    ax.legend()
    out.save(fig, "V4_scalar")


def v5_history(d, hmin):
    rows = [np.loadtxt(f) for f in sorted(Path(d).glob("postProcessing/ibmForces/*/forces.dat"),
                                          key=lambda p: float(p.parent.name))]
    h = np.vstack(rows)
    t, F = h[:, 0], (h[:, 1:4] + h[:, 4:7])/(4*hmin)
    return t, 2*F[:, 0]/0.1, 2*F[:, 1]/0.1


def v5(out, r, diag):
    import reference_dfg as dfg
    st = r["data"]["steady"]
    B = dfg.BENCH1
    fig, ax = plt.subplots(1, 3, figsize=(8.4, 2.4))
    for a, q, lab in zip(ax, ("C_D", "C_L", "dp"), ("$C_D$", "$C_L$", "$\\Delta p$")):
        Dh = [x["Dh"] for x in st]
        series(a, 0, Dh, [x[q] for x in st], "immersed")
        a.axhline(B[q], color=REF, ls="--", lw=0.9, label="reference (Nabh 1998)")
        a.set_xscale("log")
        a.set_xlabel("D/h")
        a.set_title("2D-1, " + lab)
    Ce, p = common.richardson([1/x["Dh"] for x in st], [x["C_D"] for x in st])
    ax[0].plot([st[-1]["Dh"]*2], [Ce], "x", color=C[1], label="Richardson (order %.2f)" % p)
    ax[0].legend()
    out.save(fig, "V5_2D1_convergence")

    # 2D-2: the last cycle of each run against FeatFlow level 6
    root = Path(r["cases"])
    # the one-pass runs are the diagnostic of finding 6 (code d05799c, before the sliver fix, which
    # changes V5 by 4e-4)
    one = common.RUNS / "d05799c" / "V5"
    runs = [(one / "re100_k0.nOuter1", 0.0025, "D/h 40, 1 pass, $\\Delta t$ (d05799c)"),
            (diag / "dt2", 0.0025, "D/h 40, 1 pass, $\\Delta t/2$ (d05799c)"),
            (root / "re100_k0", 0.0025, "D/h 40, 3 outer correctors"),
            (one / "re100_k1.nOuter1", 0.00125, "D/h 80, 1 pass (d05799c)"),
            (root / "re100_k1", 0.00125, "D/h 80, 3 outer correctors")]
    tr, cdr, clr, _ = dfg.bench2_cycle()
    R = dfg.bench2_summary()[6]
    fig, ax = plt.subplots(1, 3, figsize=(9.6, 2.5), gridspec_kw={"width_ratios": [1.2, 1.2, 1]})
    ax[0].plot(tr*R["St"]/0.1, clr, color=REF, lw=2.6, alpha=0.35, label="FeatFlow level 6")
    ax[1].plot(tr*R["St"]/0.1, cdr, color=REF, lw=2.6, alpha=0.35)
    dev, names, done = [], [], []
    for i, (d, hmin, lab) in enumerate(runs):
        if not d.exists():
            continue
        t, cd, cl = v5_history(d, hmin)
        complete = t[-1] > 16 - 1e-3
        s = t > t[-1] - 3.0
        c = dfg.cycle_statistics(t[s], cd[s], cl[s])
        T = 0.1/c["St"]
        w = (t >= c["t0"]) & (t <= c["t0"] + T)
        lab2 = lab + ("" if complete else " (running, t = %.1f)" % t[-1])
        ax[0].plot((t[w] - c["t0"])/T, cl[w], color=C[i], lw=1.0, label=lab2)
        ax[1].plot((t[w] - c["t0"])/T, cd[w], color=C[i], lw=1.0)
        dev.append([100*(c[q]/R[q] - 1) for q in ("St", "C_Dmax", "C_Lmax")])
        names.append(lab2)
        done.append(bool(complete))
    ax[0].set_xlabel("phase $(t - t_0)/T$")
    ax[1].set_xlabel("phase $(t - t_0)/T$")
    ax[0].set_ylabel("$C_L$")
    ax[1].set_ylabel("$C_D$")
    ax[0].set_title("2D-2: last cycle, $C_L$")
    ax[1].set_title("2D-2: last cycle, $C_D$")
    dev = np.array(dev)
    y = np.arange(len(names))
    w = 0.26
    for j, (q, lab) in enumerate((("St", "St"), ("C_Dmax", "$C_{D,max}$"), ("C_Lmax", "$C_{L,max}$"))):
        ax[2].barh(y + (j - 1)*w, dev[:, j], height=w*0.9, color=C[[6, 2, 3][j]], label=lab, zorder=3)
    ax[2].axvline(0, color=REF, lw=0.6)
    ax[2].set_yticks(y)
    ax[2].set_yticklabels([str(k + 1) for k in y])
    ax[2].invert_yaxis()
    ax[2].set_xlabel("deviation from FeatFlow [%]")
    ax[2].set_title("runs 1-%d (legend order)" % len(names))
    ax[2].legend(loc="upper left", bbox_to_anchor=(1.02, 1))
    ax[0].legend(loc="upper left", bbox_to_anchor=(0, -0.28), ncol=2)
    out.save(fig, "V5_2D2_time_accuracy")
    return list(zip(names, dev.tolist(), done))


def v7(out, r):
    from channelSA_compare import load
    root = Path(r["cases"])
    names = [("bf_wr", "body-fitted"), ("ibm_wr_f0.05_ma0.5", "immersed, f 0.05"),
             ("ibm_wr_f0.3_ma0.5", "immersed, f 0.3"), ("ibm_wr_f0.5_ma0.5", "immersed, f 0.5"),
             ("ibm_wr_f0.3_ma0", "immersed, f 0.3, unlinked")]
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6))
    bf = None
    for i, (n, lab) in enumerate(names):
        d = load(str(root / n))
        m = d["bottom"]
        o = np.argsort(d["yp"][m])
        yp, up = d["yp"][m][o], d["up"][m][o]
        if bf is None:
            bf = (yp, up)
            ax[0].plot(yp, up, color=REF, lw=2.6, alpha=0.35, label=lab)
            continue
        ax[0].plot(yp, up, color=C[i - 1], lw=0, marker=M[i - 1], ms=2.5, label=lab)
        ax[1].plot(yp, up - np.interp(yp, *bf), color=C[i - 1], marker=M[i - 1], ms=2.5, lw=0.8, label=lab)
    yy = np.logspace(0, 3, 50)
    ax[0].plot(yy, np.log(yy)/0.41 + 5.0, color=REF, ls=":", lw=0.8, label="log law (0.41, 5.0)")
    for a in ax:
        a.set_xscale("log")
        a.set_xlabel("$y^+$")
    ax[0].set_ylabel("$u^+$")
    ax[1].set_ylabel("$u^+ - u^+_{bf}$")
    ax[0].set_title("lower wall, $Re_\\tau \\approx$ %.0f" % (bf[0].max()))
    ax[1].set_title("difference to the body-fitted profile")
    ax[0].legend(loc="upper left")
    out.save(fig, "V7_channel")


def v12_history(d):
    rows = [np.loadtxt(f, ndmin=2) for f in sorted(Path(d).glob("postProcessing/ibmForces/*/forces.dat"),
                                                   key=lambda p: float(p.parent.name))]
    h = np.vstack(rows)
    t, F = h[:, 0], (h[:, 1:4] + h[:, 4:7])/0.1
    return t, 2*F[:, 0], 2*F[:, 1]


def v12(out, r, diag):
    import reference_rotatingCylinder as rot
    ref = rot.figure11()
    bf = rot.bodyFitted()
    rows = r["data"]["rows"]
    fig, ax = plt.subplots(1, 3, figsize=(9.6, 2.6))
    ok = np.isfinite(ref["C_L"]) & (ref["alpha"] <= 2.6)
    A = sorted(bf)
    for a, q, lab in ((ax[0], "C_L", "mean $C_L$"), (ax[1], "C_D", "mean $C_D$")):
        a.plot(ref["alpha"][ok], ref[q][ok], color=REF, lw=2.4, alpha=0.35, label="Mittal & Kumar 2003, fig. 11")
        a.plot(A, [bf[x][q] for x in A], lw=0, color=REF, marker="o", mfc="none", ms=8,
               label="body-fitted, M 128")
        for i, Dh in enumerate((40, 80, 160)):
            s = sorted([x for x in rows if x["Dh"] == Dh], key=lambda x: x["alpha"])
            if s:
                a.plot([x["alpha"] for x in s], [x[q] for x in s], lw=0, color=C[i], marker=M[i], ms=4,
                       label="D/h %d" % Dh)
        a.set_xlabel("$\\alpha = a\\omega/U$")
        a.set_title(lab)
    # alpha = 1, D/h = 80: history and its mean over whole shedding cycles in the last 50 time units
    t, cd, cl = v12_history(Path(r["cases"]) / "alpha1_Dh80")
    T = 1/r["data"]["metrics"]["a1_Dh80_St"]
    w0 = t[-1] - T*np.floor(50/T)
    m = [x for x in rows if (x["alpha"], x["Dh"]) == (1, 80)][0]["C_L"]
    ax[2].plot(t[t >= w0 - 10], cl[t >= w0 - 10], color=C[1], lw=0.7, label="IBM D/h 80: mean %.3f" % m)
    ax[2].plot([w0, t[-1]], [m, m], color=C[1], lw=1.6)
    ax[2].axhline(bf[1.0]["C_L"], color=REF, lw=1.2, label="body-fitted M 128: mean %.3f" % bf[1.0]["C_L"])
    ax[2].axhline(-2.265, color=REF, ls="--", lw=0.9, label="fig. 11 (digitised): -2.265")
    ax[2].set_ylim(-3.3, -1.5)
    ax[2].set_xlabel("t U/D")
    ax[2].set_ylabel("$C_L$")
    ax[2].set_title("$\\alpha$ = 1, D/h 80 (finding 7: the reference)")
    ax[0].legend(loc="lower left")
    ax[2].legend(loc="upper left", bbox_to_anchor=(1.02, 1))
    out.save(fig, "V12_rotating")


def v13(out, r):
    import reference_oscillatingCylinder as osc
    rows = r["data"]["rows"]
    fig, ax = plt.subplots(1, 3, figsize=(9.2, 2.5))
    for a, q, lab in ((ax[0], "c_d", "$c_d$"), (ax[1], "c_i", "$c_i$")):
        Dh = [x["Dh"] for x in rows]
        series(a, 0, Dh, [x[q] for x in rows], "immersed")
        a.axhline(osc.REFERENCE[q], color=REF, ls="--", lw=0.9, label="Dütsch et al. set C")
        a.axhspan(osc.REFERENCE[q]*0.98, osc.REFERENCE[q]*1.02, color="#ececec", lw=0, label="±2 %")
        e, p = common.richardson([1/x["Dh"] for x in rows], [x[q] for x in rows])
        a.plot([Dh[-1]*2], [e], "x", color=C[1], label="Richardson (order %.2f)" % p)
        a.set_xscale("log")
        a.set_xlabel("D/h")
        a.set_title(lab)
    ax[0].legend(loc="lower left")
    A, F = 5/(2*np.pi), 0.2
    d = Path(r["cases"]) / ("Dh%d" % rows[-1]["Dh"])
    rows_ = [np.loadtxt(f, ndmin=2) for f in sorted(d.glob("postProcessing/ibmForces/*/forces.dat"),
                                                    key=lambda p: float(p.parent.name))]
    h = np.vstack(rows_)
    t, Fx = h[:, 0], (h[:, 1] + h[:, 4])/0.1
    w = 2*np.pi*F
    xd, xdd = -A*w*np.cos(w*t), A*w*w*np.sin(w*t)
    cd, ci = rows[-1]["c_d"], rows[-1]["c_i"]
    ax[2].plot(t*F, Fx, color=C[0], lw=0.8, label="D/h %d" % rows[-1]["Dh"])
    ax[2].plot(t*F, -0.5*cd*xd*np.abs(xd) - 0.25*np.pi*ci*xdd, color=REF, ls="--", lw=0.8, label="Morison fit")
    ax[2].set_xlim(4, 10)
    w = t*F >= 4
    ax[2].set_ylim(1.15*Fx[w].min(), 1.15*Fx[w].max())
    ax[2].set_xlabel("t/T")
    ax[2].set_ylabel("$F_x$")
    ax[2].set_title("in-line force (fit over cycles 5-10)")
    ax[2].legend(loc="upper left", bbox_to_anchor=(1.02, 1))
    out.save(fig, "V13_oscillating")


def v15(out, r):
    import reference_tmr as tmr
    d = Path(r["cases"]) / "bf_273"
    t = latest_time(str(d))
    _, pt = read_field(str(d / t / "wallShearStress"))
    _, pc = read_field(str(d / "0" / "C"))
    x, cf = pc["plate"][:, 0], 2*np.abs(pt["plate"][:, 0])
    o = np.argsort(x)
    rows = r["data"]["rows"]
    fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.4), gridspec_kw={"width_ratios": [1.6, 1]})
    ax[0].plot(x[o], cf[o], color=C[0], label="body-fitted 273 x 193")
    inc = tmr.flatplate_cf()["FUN3D incompressible SA"]
    ax[0].plot([tmr.FLATPLATE["x_cf"]], [inc[1][0]], "o", color=REF, ms=5, mfc="none",
               label="TMR FUN3D incompressible, $C_f(0.97)$")
    ax[0].set_ylim(0.002, 0.006)
    ax[0].set_xlabel("x")
    ax[0].set_ylabel("$C_f$")
    ax[0].set_title("skin friction (IBM: no wall-shear output yet)")
    ax[0].legend()
    lab = ["body-fitted", "immersed f 0.5"]
    v = [next(x["C_D"] for x in rows if not x["ibm"]), next(x["C_D"] for x in rows if x["ibm"])]
    ax[1].axhspan(v[0]*0.99, v[0]*1.01, color="#ececec", lw=0, label="±1 % of body-fitted")
    for i, vi in enumerate(v):
        ax[1].plot([i], [vi], marker=M[i], color=C[i], ms=6, lw=0)
        ax[1].annotate("%.5f" % vi, (i, vi), xytext=(8, 0), textcoords="offset points", va="center", fontsize=6.5)
    ax[1].set_xticks([0, 1])
    ax[1].set_xticklabels(lab)
    ax[1].set_xlim(-0.5, 1.8)
    ax[1].legend(loc="lower left")
    ax[1].set_ylabel("plate $C_D$")
    ax[1].set_title("drag, 273 x 193")
    out.save(fig, "V15_flatPlate")


def main():
    commit = sys.argv[1]
    d = common.SUITE / "results" / commit
    out = Out(d)
    import report
    _, res = report.load(commit)
    summary_table(out, res)
    R = {r["id"]: r for r in res}
    for r in res:
        if r["id"] == "G9":
            table(out, r, keep=lambda c: "(mergeAlpha" not in c["name"])
            table(out, r, "G9_options", keep=lambda c: "(mergeAlpha" in c["name"])
        else:
            table(out, r)
    if "G7" in R:
        L = [r"\begin{tabular}{@{}l l r r@{}}", r"\toprule", r"tutorial & body & $F_x$ & $F_y$\\", r"\midrule"]
        for case, bodies in R["G7"]["data"]["forces"].items():
            for body, F in bodies.items():
                L.append(r"%s & %s & %.9g & %.4g\\" % (tex(case), tex(body), F[0], F[1]))
        L += [r"\bottomrule", r"\end{tabular}"]
        (d / "tables" / "G7_forces.tex").write_text("\n".join(L) + "\n")
    diag5, diag12 = common.RUNS / "diag_V5", common.RUNS / "diag_V12"
    for k, fn in (("G1", g1), ("G2", g2), ("G3", g3), ("G6", g6), ("G8", g8), ("G9", g9), ("V1", v1), ("V2", v2),
                  ("V3", v3), ("V4", v4), ("V7", v7), ("V13", v13), ("V15", v15)):
        if k in R:
            fn(out, R[k])
    flowfield.figures(out, R)
    if "V12" in R:
        v12(out, R["V12"], diag12)
    if "V5" in R:
        runs = v5(out, R["V5"], diag5)
        (d / "tables" / "V5_runs.json").write_text(json.dumps(runs, indent=1))
        L = [r"\begin{tabular}{@{}r l r r r@{}}", r"\toprule",
             r"run & & St [\%] & $C_{D,max}$ [\%] & $C_{L,max}$ [\%]\\", r"\midrule"]
        for i, (n, dv, _) in enumerate(runs):
            L.append(r"%d & %s & %+.2f & %+.2f & %+.2f\\" % (i + 1, n.replace("$\\Delta t$", "$\\Delta t$"), *dv))
        L += [r"\bottomrule", r"\end{tabular}"]
        (d / "tables" / "V5_runs.tex").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
