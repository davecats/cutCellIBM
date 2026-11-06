"""G5 rotation sources (tutorials/rotatingCylinder, eccentricRotation): the wall-flux source of
a rotating body sums to zero over the domain at every step (closure) and its maximum must not
grow against the baseline."""
import re

import common

ID, KIND, TIER = "G5", "gate", 0
TITLE = "rotation sources: closure and magnitude"
STRESSES, RELIES = ["K9"], ["G1"]
REGRESSION = (0.05, ["rotatingCylinder_maxS", "eccentricRotation_maxS"])

CASES = ["rotatingCylinder", "eccentricRotation"]
PAT = r"mass source: sum over the domain (\S+) \(zero to roundoff by closure\), max \|S\| (\S+),"


def one(ctx, name):
    d = common.copy_case(common.TUTORIALS / name, ctx.root / name)
    app = re.search(r"application\s+(\w+)", (d / "system" / "controlDict").read_text()).group(1)
    for cmd in ("blockMesh", "ibmSetGeometry", app):
        common.sh(cmd, d, "log." + cmd)
    m = re.findall(PAT, (d / ("log." + app)).read_text())
    if not m:
        raise ValueError("no mass source lines in %s" % name)
    return name, max(abs(float(s)) for s, _ in m), max(float(x) for _, x in m), len(m)


def run(ctx, res):
    rows = common.pmap(lambda n: one(ctx, n), CASES, ctx.jobs)
    res.data["metrics"] = {}
    for name, sumS, maxS, n in rows:
        res.data["metrics"].update({name + "_sumS": sumS, name + "_maxS": maxS})
        res.le("%s: max over steps |sum S|" % name, sumS, 1e-12, "%d steps" % n)
        res.report("%s: max |S|" % name, maxS)
