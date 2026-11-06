"""G6 passive-scalar budgets (tutorials/scalar): content = injected + wall fluxes + swept-volume
defect for moving bodies, source rate = wall flux rate for steady ones, closed to the
linear-solver tolerance."""
import re

import common

ID, KIND, TIER = "G6", "gate", 0
TITLE = "passive scalar budgets, static and moving"
STRESSES, RELIES = ["K11"], ["G1", "G4"]
REGRESSION = (0.05, None, 1e-12)

CASES = ["staticCylinder_zeroFlux", "staticCylinder_fixedValue", "steadyCylinder_fixedValue",
         "oscillating_zeroFlux", "oscillating_fixedValue", "galilean_zeroFlux", "multiScalar"]
NUM = r"([-+0-9.eE]+)"


MAXITER = 1000   # OpenFOAM default maxIter of the linear solvers


def steps(log, s):
    """per time step: (budget line of scalar s, whether its linear solve hit maxIter)"""
    out = []
    for blk in log.split("\nTime = ")[1:]:
        ls = [l for l in blk.splitlines() if l.startswith("scalar %s:" % s)]
        if ls:
            its = [int(x) for x in re.findall(r"Solving for %s,.*No Iterations (\d+)" % re.escape(s), blk)]
            out.append((ls[-1], bool(its) and max(its) >= MAXITER))
    return out


def one(ctx, name):
    d = common.copy_case(common.TUTORIALS / "scalar" / name, ctx.root / name)
    app = re.search(r"application\s+(\w+)", (d / "system" / "controlDict").read_text()).group(1)
    for cmd in ("blockMesh", "ibmSetGeometry", app):
        common.sh(cmd, d, "log." + cmd)
    log = (d / ("log." + app)).read_text()
    out = []
    lines = [l for l in log.splitlines() if l.startswith("scalar ")]
    for s in sorted({l.split(":")[0][7:] for l in lines}):
        st = steps(log, s)
        ls = [l for l, _ in st]
        jumps = []
        if "source rate" in ls[-1]:   # steady: rates
            e = abs(float(re.search(r"balance " + NUM, ls[-1]).group(1))) \
                / abs(float(re.search(r"source rate " + NUM, ls[-1]).group(1)))
            kind = "steady"
        else:
            # the change of the budget error per step, relative to the content; steps whose
            # solve hit maxIter are reported apart (not converged, not a conservation defect)
            err = [float(re.search(r"budget error " + NUM, l).group(1)) for l in ls]
            cont = [abs(float(re.search(r"fluid content " + NUM, l).group(1))) for l in ls]
            e = 0.0
            for i in range(len(err)):
                de = abs(err[i] - (err[i - 1] if i else 0.0))/max(cont[i], 1e-30)
                if st[i][1]:
                    jumps.append(de)
                else:
                    e = max(e, de)
            e = max(e, abs(err[-1])/max(cont[-1], 1e-30) if not jumps else 0.0)
            kind = "transient"
        solid = float(re.search(r"in the solid " + NUM, ls[-1]).group(1))
        out.append((name + ("" if s == "T" else "/" + s), kind, e, solid, len(ls), jumps))
    return out


def run(ctx, res):
    rows = [r for rs in common.pmap(lambda n: one(ctx, n), CASES, ctx.jobs) for r in rs]
    res.data["rows"] = rows
    res.data["metrics"] = {n + "_budget": e for n, _, e, _, _, _ in rows}
    for n, kind, e, solid, nsteps, jumps in rows:
        # closed to the linear-solver tolerance: transient per step relative to the content,
        # steady relative to the source rate (last iteration)
        res.le("%s: %s budget error (relative%s)" % (n, kind, ", per step" if kind == "transient" else ""),
               e, 1e-9 if kind == "transient" else 1e-8,
               "%d steps; max |T - Twall| in the solid %.1e" % (nsteps, solid))
        if jumps:
            res.report("%s: steps whose scalar solve hit maxIter %d" % (n, MAXITER), len(jumps),
                       note="budget change there up to %.2e relative (open: doc/secondOrderStudy.md, "
                            "ibmConvection by default)" % max(jumps))
