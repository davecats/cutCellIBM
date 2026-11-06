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


def write(commit):
    d, res = load(commit)
    L = ["# Validation suite, commit %s" % commit, ""]
    L += ["| id | title | status | PASS | FAIL | REPORT | time [s] |", "|---|---|---|---|---|---|---|"]
    for r in res:
        n = {s: sum(c["status"] == s for c in r["checks"]) for s in COLOURS}
        L.append("| %s | %s | **%s** | %d | %d | %d | %.0f |" % (r["id"], r["TITLE"], r["status"], n["PASS"],
                                                              n["FAIL"], n["REPORT"], r["seconds"]))
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
                ax.add_patch(Rectangle((j, i), 0.92, 0.92, color=COLOURS[r["status"]]))
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
    ax.set_title("components stressed, by status (green PASS, red FAIL, grey REPORT)", fontsize=9)
    fig.tight_layout()
    fig.savefig(d / "matrix.png", dpi=130)
    plt.close(fig)
