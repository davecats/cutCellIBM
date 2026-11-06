"""G4 free stream around a translating body (tutorials/freeStreamTranslating): fluid and body
move together at u_b = 0.37; the uniform flow must be preserved to round-off while cells
change status."""
import os

import numpy as np

import common
from foamio import latest_time, read_field

ID, KIND, TIER = "G4", "gate", 0
TITLE = "free stream around a translating body"
STRESSES, RELIES = ["K9"], ["G1"]
REGRESSION = (1.0, None, 1e-12)   # absolute: round-off level

SRC = common.TUTORIALS / "freeStreamTranslating"
UB = 0.37


def run(ctx, res):
    d = common.copy_case(SRC, ctx.root / "case")
    for cmd in ("blockMesh", "ibmSetGeometry", "ibmPimpleFoam"):
        common.sh(cmd, d, "log." + cmd)
    worst = 0.0
    for t in sorted((x for x in os.listdir(d) if x.replace(".", "").isdigit() and x != "0"), key=float):
        U, _ = read_field(os.path.join(d, t, "U"))
        a, _ = read_field(os.path.join(d, t, "ibmAlpha"))
        fl = a > 0
        worst = max(worst, float(np.abs(U[fl] - np.array([UB, 0, 0])).max()))
    res.data["metrics"] = {"max_dU": worst}
    res.le("max |U - u_b| over the fluid, all written times", worst, 1e-10,
           "latest time %s" % latest_time(str(d)))
