#!/usr/bin/env python
"""Driver of the validation suite (doc/validationSuite.md §6).

    suite.py list
    suite.py run G1 G2 ...|gates|tier1 [--jobs 8] [--quick] [--tier2]
    suite.py accept G1 ...|all [--commit c]     results/<c>/<id>.json -> baseline/<id>.json
    suite.py report [--commit c]                results/<c>/report.md, matrix.png, report.html

Needs OpenFOAM-v2412 (IBM_SUITE_FOAM, default /opt/OpenFOAM/OpenFOAM-v2412) and the project
built for the commit under test (`./Allwmake`); python with numpy, scipy, matplotlib.
"""
import argparse
import importlib
import json
import shutil
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402

GATES = ["G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9"]
TIER1 = ["V1", "V2", "V3", "V4", "V5", "V7", "V12", "V13", "V15"]
ALL = GATES + ["V%d" % i for i in range(1, 17)] + ["V6G"]


def module(cid):
    return importlib.import_module(cid.lower())


def available():
    return [c for c in ALL if (common.SUITE / ("%s.py" % c.lower())).exists()]


class Ctx:
    def __init__(self, cid, args):
        self.id = cid
        self.commit = common.commit()
        self.root = Path(args.reuse) if args.reuse else common.fresh_dir(common.RUNS / self.commit / cid)
        self.jobs = args.jobs
        self.quick = args.quick
        self.tier2 = getattr(args, "tier2", False)
        self.figdir = common.SUITE / "results" / self.commit / "figures"
        self.figdir.mkdir(parents=True, exist_ok=True)


def run(ids, args):
    out = common.SUITE / "results" / common.commit()
    out.mkdir(parents=True, exist_ok=True)
    for cid in ids:
        mod = module(cid)
        ctx = Ctx(cid, args)
        res = common.Result(mod, ctx)
        t0 = time.time()
        print("== %s %s (%s)" % (cid, mod.TITLE, ctx.root), flush=True)
        try:
            mod.run(ctx, res)
            reg = getattr(mod, "REGRESSION", (1e-6,))
            if reg is not None:
                common.regression(res, common.load_baseline(cid), *reg)
        except Exception as e:  # a crash is a FAIL of the gate, not of the driver
            res.add("run", "FAIL", None, note="%s: %s" % (type(e).__name__, e))
            traceback.print_exc()
        res.seconds = time.time() - t0
        (out / ("%s.json" % cid)).write_text(json.dumps(res.to_dict(), indent=1))
        for c in res.checks:
            print("  %-6s %-60s %s %s %s" % (c["status"], c["name"][:60], _fmt(c["value"]),
                                             ("ref " + _fmt(c["ref"])) if c["ref"] is not None else "",
                                             c["tol"] or ""), flush=True)
        print("  -> %s in %.0f s" % (res.status(), res.seconds), flush=True)


def _fmt(v):
    if isinstance(v, float):
        return "%.6g" % v
    return str(v)


def accept(ids, c):
    src = common.SUITE / "results" / c
    dst = common.SUITE / "baseline"
    dst.mkdir(exist_ok=True)
    for cid in ids:
        f = src / ("%s.json" % cid)
        if f.exists():
            r = json.loads(f.read_text())
            r["accepted_from"] = c   # the results directory: a commit, or a commit with a platform suffix
            (dst / f.name).write_text(json.dumps(r, indent=1))
            print("accepted %s from %s" % (cid, c))


def expand(ids):
    out = []
    for i in ids:
        out += {"gates": GATES, "tier1": GATES + TIER1, "all": available()}.get(i.lower(), [i.upper()])
    return [i for i in out if i in available()] if any(x in ("gates", "tier1", "all") for x in
                                                         map(str.lower, ids)) else out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("action", choices=["list", "run", "accept", "report"])
    ap.add_argument("ids", nargs="*")
    ap.add_argument("--jobs", type=int, default=8, help="concurrent OpenFOAM processes")
    ap.add_argument("--quick", action="store_true", help="coarse meshes only (driver test)")
    ap.add_argument("--tier2", action="store_true", help="tier-1 modules with their tier-2 runs as well")
    ap.add_argument("--commit", default=None)
    ap.add_argument("--reuse", default=None, help="case root of an earlier run of ONE id: re-evaluate, "
                    "skipping the commands whose log ends with End (same commit only)")
    a = ap.parse_args()
    if a.action == "list":
        for cid in available():
            m = module(cid)
            print("%-4s %-6s tier %s  %s" % (cid, m.KIND, m.TIER, m.TITLE))
    elif a.action == "run":
        if a.reuse:
            assert len(expand(a.ids)) == 1, "--reuse takes one id"
            common.REUSE = True
        run(expand(a.ids), a)
    elif a.action == "accept":
        accept(expand(a.ids or ["all"]), a.commit or common.commit())
    elif a.action == "report":
        import report
        report.write(a.commit or common.commit())


if __name__ == "__main__":
    main()
