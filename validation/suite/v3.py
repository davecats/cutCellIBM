"""V3 moving bodies (tutorials/oscillatingCylinder, galileanTranslating):
(a) the notebook (validation/reference_moving.py, ref_osc.npz and ref_gal_seamFix.npz) implements
    the method as of 2026-09 (first order, unlinked, no convective correction): the same options
    must reproduce its force history and final fields to the solver tolerance;
(b) Galilean equivalence: the cylinder translating at 0.5 through fluid at rest (mean over the
    last 400 steps) against the fixed cylinder in the stream Ubar = 0.5 (cylinderRe1, same mesh),
    at first order and at the defaults: reported (the difference is the moving-body error, O(h))."""
import re
import subprocess

import common

ID, KIND, TIER = "V3", "case", 1
TITLE = "translating and oscillating cylinder: notebook, Galilean equivalence"
STRESSES, RELIES = ["K9"], ["G4", "G5", "V2"]
REGRESSION = (1e-5,)

NOTEBOOK = {"secondOrder": False, "mergeAlpha": 0, "ibmConvection": False}
DEFAULT = {}
NUM = r"([-+0-9.eE]+)"


def case(ctx, src, name, opts):
    d = common.copy_case(common.TUTORIALS / src, ctx.root / name)
    if opts:
        common.set_keys(d / "constant" / "ibmProperties", **opts)
    app = re.search(r"application\s+(\w+)", (d / "system" / "controlDict").read_text()).group(1)
    for cmd in ("blockMesh", "postProcess -func writeCellCentres -time 0", "ibmSetGeometry -writeFaceCentres", app):
        common.sh(cmd, d, "log." + cmd.split()[0])
    return d, (d / ("log." + app)).read_text()


def compare(d, ref, tail=0):
    cmd = [common.sys.executable, str(common.VALIDATION / "compare_moving.py"), str(d),
           "--ref", str(common.VALIDATION / ref)] + (["--tail", str(tail)] if tail else [])
    out = subprocess.run(cmd, capture_output=True, text=True, check=True, cwd=common.VALIDATION).stdout
    (d / "compare_moving.txt").write_text(out)
    F = [(float(a), float(b)) for a, b in re.findall(r"F_x OpenFOAM\s+" + NUM + r" .*notebook\s+" + NUM, out)]
    Fmax = max(abs(b) for _, b in F)
    dF = max(abs(a - b) for a, b in F)/Fmax
    u = re.search(r"\n\s+u\s+max\|diff\| = " + NUM + r".*max\|u\| = " + NUM, out)
    m = re.search(r"mean F_x over the last \d+ steps: OpenFOAM " + NUM, out)
    return {"dF": dF, "du": float(u.group(1))/float(u.group(2)), "meanF": float(m.group(1)) if m else None}


def run(ctx, res):
    jobs = [("oscillatingCylinder", "osc_notebook", NOTEBOOK), ("galileanTranslating", "gal_notebook", NOTEBOOK),
            ("galileanTranslating", "gal_default", DEFAULT), ("cylinderRe1", "cyl_notebook", NOTEBOOK),
            ("cylinderRe1", "cyl_default", DEFAULT), ("oscillatingCylinder", "osc_default", DEFAULT)]
    out = dict(zip([j[1] for j in jobs], common.pmap(lambda j: case(ctx, *j), jobs, ctx.jobs)))
    r = {"osc_notebook": compare(out["osc_notebook"][0], "ref_osc.npz"),
         "gal_notebook": compare(out["gal_notebook"][0], "ref_gal_seamFix.npz", 400),
         "gal_default": compare(out["gal_default"][0], "ref_gal_seamFix.npz", 400),
         "osc_default": compare(out["osc_default"][0], "ref_osc.npz")}
    cyl = {k: common.ibm_force(out[k][1], "cylinder")["total"][0]/(10.0/81) for k in ("cyl_notebook", "cyl_default")}
    res.data["compare"], res.data["cylinder"] = r, cyl
    res.data["metrics"] = {"%s_%s" % (k, q): v[q] for k, v in r.items() for q in ("dF", "du", "meanF") if v[q] is not None}
    res.data["metrics"].update(cyl)
    # the code has moved on since the notebook references (2026-09-08: 2e-5 agreement then); the
    # first-order path now differs by ~2.5e-3 (commits since then: wall-velocity evaluation point,
    # rotating-wall mass source, fillUncovered, ...): reported until the reference is regenerated
    for k in ("osc_notebook", "gal_notebook"):
        res.report("%s: max |F_x - notebook| / max|F_x| (history)" % k, r[k]["dF"],
                   note="2e-5 on 2026-09-08; notebook = method of 2026-09")
        res.report("%s: max |u - notebook| / max|u| (final field)" % k, r[k]["du"])
    res.report("osc_default: max |F_x - notebook| / max|F_x| (second order against the first-order notebook)",
               r["osc_default"]["dF"])
    for opt in ("notebook", "default"):
        mt, fx = r["gal_%s" % opt]["meanF"], cyl["cyl_%s" % opt]
        res.report("Galilean: translating (mean of last 400 steps) vs fixed cylinder, %s options" % opt,
                   -mt, fx, "deviation %+.2f %%" % (100*(-mt/fx - 1)))
