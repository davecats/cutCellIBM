"""Run times of the immersed cases against their body-fitted counterparts, from the solver logs of the
stored case directories (the cost of the method next to its accuracy).

Per run: solver, pressure solver, host, ranks, cells (fluid cells for the immersed runs: total minus
solid), time steps or SIMPLE iterations, the master's CPU time (OpenFOAM ExecutionTime), the pressure
iterations per step. Cost per cell and step = ExecutionTime x ranks / (fluid cells x steps) in
core-microseconds, from the difference of the first and last ExecutionTime of each log (start-up
excluded); the master's CPU time stands for every rank's (balanced decompositions).

    timing.py <commit> [ID=commit ...]   print the pairs of that commit's runs (~/runs/suite/<commit>),
                                         ID from another commit's runs
    timing.py bench           print the same-host benchmark (~/runs/timing_bench)

The pairs name the body-fitted and immersed runs of the same problem. Only G9 compares the same mesh with
the same solver and settings; V7 and V15 the same grid. Elsewhere meshes, solver settings, OpenFOAM
versions, hosts or rank counts differ. The suite runs several cases at once (--jobs), so a time carries the
load of its host at that moment: the same work took 4.5-6.5 times longer in the tier-1 run on an
overloaded host (V4 vs V2). Only runs on the same host at the same time compare costs; the benchmark
(~/runs/timing_bench/bench.sh) continues the transient pairs of V5 and V12 from their final states, on
one host, alternately, twice."""
import json
import re
import sys
from pathlib import Path

import common

STEP = re.compile(r"^Time = ", re.M)
EXEC = re.compile(r"ExecutionTime = ([0-9.eE+-]+) s\s+ClockTime = ([0-9.eE+-]+) s")
PIT = re.compile(r"Solving for p, .*No Iterations (\d+)")
CONV = re.compile(r"solution converged in (\d+) iterations")
BENCH = Path.home() / "runs" / "timing_bench"


def cells(d):
    """total and solid cells: ibmSetGeometry's count if the run has one, else the mesh logs"""
    for name in ("log.ibmSetGeometry.linkZone", "log.ibmSetGeometry"):
        f = d / name
        if f.exists():
            m = re.findall(r"cells: (\d+) total, (\d+) solid", f.read_text(errors="replace"))
            if m:
                return int(m[-1][0]), int(m[-1][1])
    for f in (d / "log.checkMesh", d / "log.decomposePar"):
        if f.exists():
            m = re.findall(r"^\s+cells:\s+(\d+)|Number of cells = (\d+)", f.read_text(errors="replace"), re.M)
            if m:
                return int(m[-1][0] or m[-1][1]), 0
    f = d / "constant/polyMesh/owner"
    if f.exists():
        with open(f, errors="replace") as fh:
            m = re.search(r"nCells:\s*(\d+)", fh.read(4096))
        if m:
            return int(m.group(1)), 0
    return None, 0


def psolver(d):
    """pressure solver of system/fvSolution, short: GAMG, PCG/DIC, PCG/GAMG (base of linkedGroup)"""
    f = d / "system/fvSolution"
    if not f.exists():
        return "?"
    t = re.sub(r"//[^\n]*", "", f.read_text(errors="replace"))
    m = re.search(r'(?:^|[\s{;])(?:p|"\(?p[|\w]*\)?")\s*\{', t)
    if not m:
        return "?"
    i, depth = m.end(), 1
    while depth and i < len(t):
        depth += {"{": 1, "}": -1}.get(t[i], 0)
        i += 1
    b = t[m.end():i]
    s = re.search(r"solver\s+(\w+)", b).group(1)
    if s != "PCG":
        return s
    base = re.search(r"base\s*\{?\s*(?:preconditioner\s+)?(\w+)", b) or re.search(r"preconditioner\s+(\w+)", b)
    return "PCG/" + (base.group(1) if base else "?")


CACHE = common.RUNS / "timing_cache.json"   # parsed logs by path, size and mtime (the logs are large)
_cache = None
PARSER = 2   # raise when _parse changes: older cache entries are ignored


def parse(d, logs=None):
    """timing record of one case directory (all its solver logs as one run, or the given logs), or None"""
    global _cache
    d = Path(d).expanduser()
    if logs is None:
        logs = sorted(p for p in d.glob("log.*Foam*") if re.fullmatch(r"log\.(ibm)?(simple|pimple)Foam(\.\d+)?", p.name, re.I))
    if not logs:
        return None
    if _cache is None:
        try:
            _cache = json.loads(CACHE.read_text())
        except (OSError, ValueError):
            _cache = {}
    key = "v%d|" % PARSER + "|".join("%s:%d:%d" % (f, f.stat().st_size, f.stat().st_mtime) for f in logs)
    if key not in _cache:
        _cache[key] = _parse(d, logs)
        try:
            CACHE.write_text(json.dumps(_cache))
        except OSError:
            pass
    return _cache[key]


def _parse(d, logs):
    r = dict(dir=str(d), name=d.name, app=logs[0].name.split(".")[1], steps=0, timed=0, exec=0.0, dexec=0.0,
             clock=0.0, pit=0, host=None, ranks=1, converged=None, finished=False, p=psolver(d), version="?")
    for f in logs:   # a continuation log restarts ExecutionTime from zero
        t = f.read_text(errors="replace")
        if r["host"] is None:
            m = re.search(r"^Host\s*:\s*(\S+)", t, re.M)
            r["host"] = m.group(1).replace("ISTM-", "") if m else "?"
            m = re.search(r"^nProcs\s*:\s*(\d+)", t, re.M)
            r["ranks"] = int(m.group(1)) if m else 1
            m = re.search(r"OPENFOAM=(\d+)", t)
            r["version"] = "v" + m.group(1) if m else "?"
        r["steps"] += len(STEP.findall(t))
        e = [float(x[0]) for x in EXEC.findall(t)]
        if e:
            r["exec"] += e[-1]
            r["clock"] += float(EXEC.findall(t)[-1][1])
            r["dexec"] += e[-1] - e[0]
            r["timed"] += len(e) - 1
        r["pit"] += sum(int(x) for x in PIT.findall(t))
        m = CONV.search(t)
        if m:
            r["converged"] = int(m.group(1))
        r["finished"] = bool(re.search(r"^End\s*(Finalising parallel run\s*)?$", t[-200:], re.M))
    n, solid = cells(d)
    r["cells"] = n
    r["fluid"] = n - solid if n else None
    if r["timed"] and r["fluid"]:
        r["us"] = 1e6*r["dexec"]*r["ranks"]/(r["fluid"]*r["timed"])
    if r["steps"]:
        r["pit_step"] = r["pit"]/r["steps"]
    return r


def pairs(commit, other=None):
    """(case, problem, body-fitted dirs, immersed dirs, note) of the runs stored for commit; other: {ID: commit}
    for cases taken from another commit's runs"""
    R = common.RUNS / commit
    other = other or {}
    bf = common.RUNS / "bodyfitted"
    v2bf = Path.home() / "runs" / "ibmSA_study" / "freeslip_bf"   # the V2 references (v2.py)
    P = []
    for N in (20, 40, 80, 160):
        P.append(("G9", "Hiemenz, level wall $f$ 0.5, N = %d" % N,
                  [R / "G9" / ("hiemenz_t0_f0.5_N%d_bf_ma0.5_c1" % N)], [R / "G9" / ("hiemenz_t0_f0.5_N%d_ma0.5_c1" % N)],
                  "same mesh, solver, settings"))
    P.append(("V2", "cylinder array, no slip", [v2bf / ("m%d_noSlip" % M) for M in (64, 128, 256)],
              [common.RUNS / other.get("V2", commit) / "V2" / ("noSlip_N%d_ma0.5" % N) for N in (41, 81, 161, 321)],
              "O-grid vs Cartesian; residuals 1e-10 vs 1e-9; M 256 stopped at 20000"))
    P.append(("V4", "passive scalar, cylinder", [R / "V4" / "bf_M64", R / "V4" / "bf_M128"],
              [R / "V4" / ("ibm_N%d" % N) for N in (41, 81, 161)], "O-grid vs Cartesian"))
    P.append(("V5", "DFG 2D-2, D/h 40, 3 outer correctors", [bf / "dfg_M64_nO3"], [R / "V5" / "re100_k0"],
              "O-grid M 64 vs Cartesian, same dt"))
    P.append(("V7", "SA channel, wall-resolved", [R / "V7" / "bf_wr"], sorted((R / "V7").glob("ibm_wr_*")),
              "same grid, wall row cut at $f$"))
    for a in (0, 1, 2):
        P.append(("V12", "rotating cylinder, $\\alpha$ = %d" % a, [bf / ("rot_a%d_M%d" % (a, M)) for M in (64, 128)],
                  sorted((R / "V12").glob("alpha%d_Dh*" % a), key=lambda p: int(p.name.split("Dh")[1])),
                  "O-grid (adaptive dt) vs Cartesian refined boxes"))
    P.append(("V15", "TMR flat plate 273 $\\times$ 193", [R / "V15" / "bf_273"], [R / "V15" / "ibm_273_f0.5_ma0.5"],
              "same grid, wall row cut at $f$ 0.5"))
    return P


def collect(commit, other=None):
    out = []
    for case, prob, B, I, note in pairs(commit, other):
        b = [x for x in map(parse, B) if x]
        i = [x for x in map(parse, I) if x]
        if b and i:
            out.append(dict(case=case, problem=prob, note=note, bf=b, ibm=i))
    return out


def bench():
    """same-host benchmark: {name: [record per round]}"""
    out = {}
    for d in sorted(BENCH.glob("*_v*")):
        recs = [parse(d, [f]) for f in sorted(d.glob("log.*Foam.r*")) if re.fullmatch(r"log\.\w+Foam\.r\d+", f.name)]
        recs = [x for x in recs if x and x["finished"]]
        if recs:
            out[d.name] = recs
    return out


BENCH_LABEL = {"bf_v12": r"V12 $\alpha$ = 1, body-fitted M 64", "ibm_v12": r"V12 $\alpha$ = 1, immersed D/h 80",
               "bf_v5": "V5 2D-2, body-fitted M 64", "ibm_v5": "V5 2D-2, immersed D/h 40"}


def short(name):
    """run name without the parts its pair heading already gives"""
    return name.replace("hiemenz_t0_f0.5_", "").replace("_ma0.5_c1", "")


def fmt_n(n):
    return "{:,}".format(n).replace(",", r"\,")


def write_tex(tables, commit, other=None):
    """tables/timings.tex (pairs) and tables/timings_bench.tex (same-host benchmark)"""
    L = [r"\begin{longtable}{@{}l l l l l r r r r r r@{}}", r"\toprule",
         r"case & run & solver & $p$ & host & ranks & cells & steps & CPU [s] & $\mu$s & $p$ it.\\",
         r"\midrule", r"\endhead"]
    for p in collect(commit, other):
        L.append(r"\multicolumn{11}{@{}l}{\textbf{%s} %s \quad\textit{%s}}\\" % (p["case"], p["problem"], p["note"]))
        for k in ("bf", "ibm"):
            for r in p[k]:
                L.append(r"& %s & %s & %s & %s & %d & %s & %d & %s & %.3g & %.1f\\"
                         % (short(r["name"]).replace("_", r"\_"), r["app"], r["p"], r["host"], r["ranks"], fmt_n(r["fluid"]),
                            r["steps"], fmt_n(round(r["exec"])) if r["exec"] >= 100 else "%.3g" % r["exec"],
                            r.get("us", float("nan")), r.get("pit_step", float("nan"))))
        L.append(r"\addlinespace[2pt]")
    L += [r"\bottomrule", r"\end{longtable}"]
    (tables / "timings.tex").write_text("\n".join(L) + "\n")
    B = bench()
    if not B:
        return
    L = [r"\begin{tabular}{@{}l l l r r r r r@{}}", r"\toprule",
         r"run (cetus, v2412) & solver & $p$ & ranks & cells & steps & $\mu$s (rounds 1, 2) & $p$ it.\\", r"\midrule"]
    for n in sorted(B, key=lambda n: (n.split("_")[1][::-1], n)):   # V12 then V5, body-fitted first
        recs = B[n]
        r = recs[0]
        L.append(r"%s & %s & %s & %d & %s & %d & %s & %.1f\\"
                 % (BENCH_LABEL.get(n, n.replace("_", r"\_")), r["app"], r["p"], r["ranks"], fmt_n(r["fluid"]), r["timed"],
                    ", ".join("%.3g" % x["us"] for x in recs), r["pit_step"]))
    L += [r"\bottomrule", r"\end{tabular}"]
    (tables / "timings_bench.tex").write_text("\n".join(L) + "\n")


def show(r, k=""):
    print("  %-3s %-34s %-14s %-9s %-6s %-8s %3d ranks %9s cells %7d steps %10.1f s CPU %7.3g us %7.1f p-it"
          % (k, r["name"], r["app"], r["p"], r["version"], r["host"], r["ranks"], r["fluid"], r["steps"], r["exec"],
             r.get("us", float("nan")), r.get("pit_step", float("nan"))))


def main():
    if sys.argv[1] == "bench":
        for n, recs in bench().items():
            for r in recs:
                show(r)
        return
    other = dict(a.split("=") for a in sys.argv[2:])
    for p in collect(sys.argv[1], other):
        print("== %s %s (%s)" % (p["case"], p["problem"], p["note"]))
        for k in ("bf", "ibm"):
            for r in p[k]:
                show(r, k)


if __name__ == "__main__":
    main()
