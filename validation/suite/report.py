"""Report of one commit's results (doc/validationSuite.md §6): results/<commit>/report.md (one
line per check: value, reference, tolerance, status) and matrix.png (components x gates/cases,
coloured by status)."""
import json
from pathlib import Path

import common

COMPONENTS = ["K%d" % i for i in range(1, 16)]
COLOURS = {"PASS": "#2e7d32", "FAIL": "#c62828", "REPORT": "#9e9e9e"}


def load(commit):
    d = common.SUITE / "results" / commit
    order = {c: i for i, c in enumerate(["G%d" % i for i in range(1, 10)] + ["V%d" % i for i in range(1, 17)])}
    res = [json.loads(f.read_text()) for f in d.glob("*.json")]
    return d, sorted(res, key=lambda r: order.get(r["id"], 99))


def fmt(v):
    if v is None:
        return "-"
    if isinstance(v, float):
        return "%.4g" % v
    if isinstance(v, list):
        return "(" + " ".join(fmt(x) for x in v) + ")"
    return str(v)


def is_regression(c):
    """change detection against the stored baseline (the regression check of every case, the
    per-tutorial force and moment comparison of G7), as opposed to checks against exact or reference data"""
    return c["name"].startswith("regression against baseline") or "F_base" in c["name"]


def reference_status(r):
    """status of the checks against exact and reference data (FAIL, PASS, or REPORT if none decides)"""
    st = [c["status"] for c in r["checks"] if c["name"] != "baseline" and not is_regression(c)]
    return "FAIL" if "FAIL" in st else "PASS" if "PASS" in st else "REPORT"


def baseline_status(r):
    """status of the change against the stored baseline, None without such a check"""
    st = [c["status"] for c in r["checks"] if is_regression(c)]
    return ("FAIL" if "FAIL" in st else "PASS") if st else None


def write(commit):
    d, res = load(commit)
    L = ["# Validation suite, commit %s" % commit, ""]
    L += ["Status: the checks against exact and reference data; the change against the stored baseline",
          "(regression) separately: a FAIL there is a change, intended or not, not a wrong result.", "",
          "| id | title | checks | vs. baseline | PASS | FAIL | REPORT | time [s] |",
          "|---|---|---|---|---|---|---|---|"]
    for r in res:
        n = {s: sum(c["status"] == s for c in r["checks"]) for s in COLOURS}
        L.append("| %s | %s | **%s** | %s | %d | %d | %d | %.0f |" % (
            r["id"], r["TITLE"], reference_status(r), baseline_status(r) or "-", n["PASS"], n["FAIL"], n["REPORT"],
            r["seconds"]))
    for r in res:
        L += ["", "## %s %s: %s" % (r["id"], r["TITLE"], r["status"]), "",
              "Stresses %s; relies on %s. Cases: `%s`." % (" ".join(r["STRESSES"]), " ".join(r["RELIES"]) or "-",
                                                          r["cases"]), "",
              "| status | check | value | reference | tolerance | note |", "|---|---|---|---|---|---|"]
        for c in r["checks"]:
            L.append("| %s | %s | %s | %s | %s | %s |" % (c["status"], c["name"], fmt(c["value"]), fmt(c["ref"]),
                                                       c["tol"] or "-", (c["note"] or "").replace("|", "/")))
        for f in r.get("figures", []):
            L += ["", "![%s](%s)" % (f, f)]
    (d / "report.md").write_text("\n".join(L) + "\n")
    matrix(d, res)
    print("wrote %s" % (d / "report.md"))


def matrix(d, res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    fig, ax = plt.subplots(figsize=(1.0 + 0.45*len(res), 0.4*len(COMPONENTS) + 1.2))
    for j, r in enumerate(res):
        for i, k in enumerate(COMPONENTS):
            if k in r["STRESSES"]:
                ax.add_patch(Rectangle((j, i), 0.92, 0.92, color=COLOURS[reference_status(r)]))
            elif k in r.get("RELIES_K", []):
                ax.add_patch(Rectangle((j, i), 0.92, 0.92, color="#e0e0e0"))
    ax.set_xlim(0, len(res))
    ax.set_ylim(len(COMPONENTS), 0)
    ax.set_xticks([j + 0.46 for j in range(len(res))])
    ax.set_xticklabels([r["id"] for r in res])
    ax.set_yticks([i + 0.46 for i in range(len(COMPONENTS))])
    ax.set_yticklabels(COMPONENTS)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("components stressed, by the checks against exact and reference data\n"
                 "(green PASS, red FAIL, grey REPORT; the change against the baseline not included)", fontsize=9)
    fig.tight_layout()
    fig.savefig(d / "matrix.png", dpi=130)
    plt.close(fig)
