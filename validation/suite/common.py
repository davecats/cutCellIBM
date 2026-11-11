"""Common machinery of the validation suite (doc/validationSuite.md §6): case directories, running
OpenFOAM, editing dictionaries, checks with PASS / FAIL / REPORT, observed orders, results.

Case directories live outside the repository in $IBM_SUITE_RUNS/<commit>/<id>/ (never deleted
or overwritten: a second run of the same id and commit goes to <id>.1, <id>.2, ...); results as
JSON in validation/suite/results/<commit>/<id>.json, the accepted baseline in
validation/suite/baseline/<id>.json.
"""
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

SUITE = Path(__file__).resolve().parent
VALIDATION = SUITE.parent
REPO = VALIDATION.parent
TUTORIALS = REPO / "tutorials"
RUNS = Path(os.environ.get("IBM_SUITE_RUNS", Path.home() / "runs" / "suite"))
# OpenFOAM v2412 (the default since 2026-11-10, as on HoreKa; v2106_prf before): set IBM_SUITE_FOAM
# for another version, e.g. "source /opt/OpenFOAM/OpenFOAM-v2106_prf/etc/bashrc; export LD_LIBRARY_PATH=..."
FOAM_ENV = os.environ.get("IBM_SUITE_FOAM",
                          "source /opt/OpenFOAM/OpenFOAM-v2412/etc/bashrc > /dev/null 2>&1; "
                          "export LD_LIBRARY_PATH=/opt/OpenFOAM/ThirdParty-v2412/platforms/"
                          "linux64GccDPInt32/lib:$LD_LIBRARY_PATH")   # kahip etc. for decomposePar
sys.path.insert(0, str(VALIDATION))


CODE = ["src", "applications", "tutorials"]


def commit():
    """the code under test: short hash of the last commit that changed src, applications or
    tutorials, with '-dirty' if they differ from it (suite-only commits keep the key)"""
    h = subprocess.run(["git", "log", "-1", "--format=%h", "--"] + CODE, cwd=REPO, capture_output=True,
                       text=True, check=True).stdout.strip()
    d = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"] + CODE, cwd=REPO,
                       capture_output=True, text=True, check=True).stdout.strip()
    return h + ("-dirty" if d else "")


def head():
    return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=REPO, capture_output=True,
                          text=True, check=True).stdout.strip()


def fresh_dir(path):
    """path, or path.1, path.2, ... whichever does not exist yet; created"""
    path = Path(path)
    p, i = path, 0
    while p.exists():
        i += 1
        p = path.with_name("%s.%d" % (path.name, i))
    p.mkdir(parents=True)
    return p


def copy_case(src, dst, extra=()):
    """copy 0, constant, system (and the names in extra) of a case into dst (created): only the
    files tracked by git (a tutorial directory may hold the output of a local run)"""
    src, dst = Path(src), Path(dst)
    if REUSE and (dst / "system").exists():
        return dst
    dst.mkdir(parents=True, exist_ok=True)
    tracked = subprocess.run(["git", "ls-files", "--", str(src)], cwd=REPO, capture_output=True, text=True,
                             check=True).stdout.split()
    tracked = {(REPO / f).resolve() for f in tracked}
    for name in ("0", "constant", "system") + tuple(extra):
        s = src / name
        files = [s] if s.is_file() else sorted(s.rglob("*")) if s.is_dir() else []
        for f in files:
            if f.is_file() and (f.resolve() in tracked or not tracked) and "polyMesh" not in f.parts:
                o = dst / f.relative_to(src)
                o.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(f, o)
    return dst


def set_keys(path, **entries):
    """set top-level 'key value;' entries of an OpenFOAM dictionary (appended if absent)"""
    path = Path(path)
    t = path.read_text()
    for k, v in entries.items():
        v = foam_value(v)
        # at the start of a line or after another entry on the same line ("L 3; N 40;")
        pat = re.compile(r"(^|;)([ \t]*%s[ \t]+)[^;{}]*;" % re.escape(k), re.M)
        if pat.search(t):
            t = pat.sub(lambda m: m.group(1) + m.group(2) + v + ";", t, count=1)
        else:
            t += "\n%s %s;\n" % (k, v)
    path.write_text(t)


def foam_value(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (list, tuple, np.ndarray)):
        return "(" + " ".join(foam_value(x) for x in v) + ")"
    if isinstance(v, (float, np.floating)):
        return repr(float(v))
    if isinstance(v, np.integer):
        return str(int(v))
    return str(v)


def header(cls, obj):
    return "FoamFile { version 2.0; format ascii; class %s; object %s; }\n" % (cls, obj)


class RunError(RuntimeError):
    pass


REUSE = False   # suite.py run --reuse DIR: re-evaluate finished runs, skip commands already done


def done(cwd, log):
    f = Path(cwd) / log
    return REUSE and f.exists() and "\nEnd" in f.read_text()[-2000:]


def sh(cmd, cwd, log, timeout=None, np_=1):
    """run cmd (a shell string) in the OpenFOAM environment of this project; output to cwd/log.
    With np_ > 1 the command is wrapped in decomposePar / mpirun -parallel."""
    cwd = Path(cwd)
    if done(cwd, log):
        return 0.0
    if np_ > 1:
        cmd = "mpirun --oversubscribe -np %d %s -parallel" % (np_, cmd)
    script = "%s; source %s/etc/bashrc > /dev/null; %s" % (FOAM_ENV, REPO, cmd)
    t0 = time.time()
    with open(cwd / log, "w") as fh:
        r = subprocess.run(["bash", "-c", script], cwd=cwd, stdout=fh, stderr=subprocess.STDOUT,
                           timeout=timeout)
    if r.returncode != 0:
        tail = "".join((cwd / log).read_text().splitlines(True)[-15:])
        raise RunError("%s failed in %s (exit %d):\n%s" % (cmd, cwd, r.returncode, tail))
    return time.time() - t0


def pmap(fn, items, jobs):
    """fn over items with at most jobs threads (each typically runs one OpenFOAM process)"""
    with ThreadPoolExecutor(max_workers=max(1, jobs)) as ex:
        return list(ex.map(fn, items))


def last_match(text, pattern, group=1, cast=float):
    m = re.findall(pattern, text)
    if not m:
        raise ValueError("no match of %r" % pattern)
    v = m[-1]
    v = v[group - 1] if isinstance(v, tuple) else v
    return cast(v)


def vec(s):
    return np.array([float(x) for x in s.strip("()").split()])


def ibm_force(log_text, body=None):
    """last 'IBM force on <body>' line: dict pressure, viscous, total (arrays)"""
    pat = r"IBM force on %s: pressure \(([^)]*)\) viscous \(([^)]*)\) total \(([^)]*)\)" % (
        re.escape(body) if body else r"\S+")
    m = re.findall(pat, log_text)
    if not m:
        raise ValueError("no IBM force line")
    p, v, t = m[-1]
    return {"pressure": vec(p), "viscous": vec(v), "total": vec(t)}


def order(h, err, last=3):
    """observed order: least-squares slope of log err over log h, the `last` finest points"""
    h, err = np.asarray(h, float), np.asarray(err, float)
    i = np.argsort(h)[:last]
    if len(i) < 2 or np.any(err[i] <= 0):
        return float("nan")
    return float(np.polyfit(np.log(h[i]), np.log(err[i]), 1)[0])


def driving_gradient(log):
    """meanVelocityForce: the driving gradient after the pressure correction of each iteration
    (the last of its two 'pressure gradient' lines), and the largest relative change over the
    last n iterations (0 when the iteration is stationary, not just its residuals small)"""
    g = [re.findall(r"pressure gradient = ([-+0-9.eE]+)", it) for it in re.split(r"\nTime = ", log)[1:]]
    g = np.array([float(x[-1]) for x in g if x])
    return g, lambda n=100: float(np.ptp(g[-n:])/abs(g[-1]))


def richardson(h, q):
    """extrapolated value and observed order from the three finest grids (h any spacing)"""
    h, q = np.asarray(h, float), np.asarray(q, float)
    i = np.argsort(h)[:3]
    h, q = h[i], q[i]
    e21, e32 = q[1] - q[0], q[2] - q[1]
    if e21 == 0 or e32 == 0 or e21*e32 < 0:
        return float("nan"), float("nan")
    r21, r32 = h[1]/h[0], h[2]/h[1]
    if abs(r21 - r32) > 1e-6:
        return float("nan"), float("nan")
    p = math.log(e32/e21)/math.log(r21)
    return float(q[0] - e21/(r21**p - 1)), float(p)


class Result:
    """checks and data of one gate or case"""

    def __init__(self, mod, ctx):
        self.id = mod.ID
        self.meta = {k: getattr(mod, k) for k in ("TITLE", "STRESSES", "RELIES", "KIND", "TIER")}
        self.commit = ctx.commit
        self.checks = []
        self.data = {}
        self.figures = []
        self.cases = str(ctx.root)
        self.started = time.strftime("%Y-%m-%d %H:%M:%S")
        self.seconds = 0.0

    def add(self, name, status, value, ref=None, tol=None, note=""):
        self.checks.append({"name": name, "status": status, "value": _plain(value),
                            "ref": _plain(ref), "tol": tol, "note": note})

    def le(self, name, value, tol, note=""):
        """PASS if value <= tol (nan: FAIL)"""
        ok = value is not None and np.isfinite(value) and value <= tol
        self.add(name, "PASS" if ok else "FAIL", value, None, "<= %g" % tol, note)

    def ge(self, name, value, tol, note=""):
        ok = value is not None and np.isfinite(value) and value >= tol
        self.add(name, "PASS" if ok else "FAIL", value, None, ">= %g" % tol, note)

    def rel(self, name, value, ref, rtol, note=""):
        """PASS if |value/ref - 1| <= rtol"""
        ok = np.isfinite(value) and abs(value/ref - 1) <= rtol
        self.add(name, "PASS" if ok else "FAIL", value, ref, "%g %%" % (100*rtol),
                 (note + "; " if note else "") + "deviation %+.3g %%" % (100*(value/ref - 1)))

    def absolute(self, name, value, ref, atol, note=""):
        ok = np.isfinite(value) and abs(value - ref) <= atol
        self.add(name, "PASS" if ok else "FAIL", value, ref, "± %g" % atol,
                 (note + "; " if note else "") + "deviation %+.3g" % (value - ref))

    def report(self, name, value, ref=None, note=""):
        self.add(name, "REPORT", value, ref, None, note)

    def status(self):
        s = [c["status"] for c in self.checks]
        return "FAIL" if "FAIL" in s else ("PASS" if "PASS" in s else "REPORT")

    def to_dict(self):
        return {"id": self.id, **self.meta, "commit": self.commit, "suite": head(), "status": self.status(),
                "started": self.started, "seconds": round(self.seconds, 1), "cases": self.cases,
                "checks": self.checks, "data": _plain(self.data), "figures": self.figures}


def _plain(v):
    if isinstance(v, dict):
        return {str(k): _plain(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_plain(x) for x in v]
    if isinstance(v, np.ndarray):
        return _plain(v.tolist())
    if isinstance(v, (np.floating, float)):
        v = float(v)
        return v if np.isfinite(v) else None
    if isinstance(v, np.integer):
        return int(v)
    return v


def regression(res, baseline, rtol, keys=None, atol=0.0):
    """add a check per stored metric (flat dict of numbers in res.data['metrics']) against the
    baseline: |value - base| <= rtol |base| + atol; REPORT if there is no baseline"""
    m = res.data.get("metrics", {})
    if baseline is None:
        res.report("baseline", len(m), note="no stored baseline: this run can be accepted as one")
        return
    b = baseline.get("data", {}).get("metrics", {})
    worst, wname, missing = 0.0, "", []
    for k in (keys or m):
        if k not in b or b[k] is None or m.get(k) is None:
            missing.append(k)
            continue
        d = abs(m[k] - b[k])/(rtol*abs(b[k]) + atol) if (rtol*abs(b[k]) + atol) > 0 else 0.0
        if d > worst:
            worst, wname = d, k
    res.add("regression against baseline %s (%d metrics)" % (baseline.get("commit"), len(m)),
            "PASS" if worst <= 1 else "FAIL", worst, None, "<= 1 (in units of the tolerance)",
            "rtol %g, atol %g; worst %s%s" % (rtol, atol, wname or "-",
                                              "; not in baseline: %s" % ", ".join(missing) if missing else ""))


def load_baseline(cid):
    p = SUITE / "baseline" / ("%s.json" % cid)
    return json.loads(p.read_text()) if p.exists() else None


def value_points(case, t, C):
    """points where the cell values live: cell centre + ibmValueOffset (the group centroid of a
    linked cell; written when cells are linked), else + ibmCentroidOffset (fluid centroid), else
    the centre (first order)"""
    from foamio import read_field
    for name in ("ibmValueOffset", "ibmCentroidOffset"):
        f = os.path.join(str(case), t, name)
        if os.path.exists(f):
            return C + read_field(f)[0]
    return C.copy()


def linked_cells(d, log="log.ibmSetGeometry"):
    """number of linked cut cells reported by ibmSetGeometry (0 without linking)"""
    m = re.search(r"cell linking: (\d+) cells", (Path(d) / log).read_text())
    return int(m.group(1)) if m else 0


def write_dict(path, cls, obj, body):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(header(cls, obj) + body)
