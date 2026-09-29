# Raw results of the two fixes of 2026-09-29

Records of the runs behind the "Fixes" section of `doc/secondOrderStudy.md` and the
corresponding sections of `doc/cutCellIBM.tex`: the rotating-wall source over the open faces
(issue 1 of `doc/openIssues/openIssues.pdf`) and the linearly exact second-order pressure
gradient (issue 2). The case directories were not kept; the scripts named below regenerate them.

- `orderTest.txt`: last lines of every `ibmOrderTest` run of
  `tutorials/order/laplaceDisk/Allrun` (problems LD, LN, CD, CN, GL, GQ; switch off/on;
  N = 20...160), the three-dimensional GL check on a sphere, and GL/GQ before the fix.
- `taylorCouette.txt`: `validation/taylorCouette_error.py` for `tutorials/order/taylorCouette`
  at N = 40, 80, 160 before the fixes, with each fix alone and with both.
- `forces_tc_*.dat`: `postProcessing/ibmForces/.../forces.dat` of the runs with both fixes
  (second order) and with the source fix (first order).
- `hydrostatic.txt`: `validation/hydrostatic_error.py` for `tutorials/order/hydrostatic`
  (fluid at rest in an immersed container under a body force), before and after, switch off/on.
- `richardson.txt`: iterations, convergence flag and final `IBM force` line of the
  `validation/richardson_drag.sh` runs with both fixes (N = 41, 81, 161, 321).
