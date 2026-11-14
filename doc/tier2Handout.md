# Handout: from the tier-1 baseline at 8d416fd to tier 2 (written 2026-11-14)

Branch `feature_iddes_from_scalar_secondOrder`, pushed to origin at 690bd53 (this handout comes
after). Main tree `~/Codes/openFoam.ibm`, built with OpenFOAM v2412 at 8d416fd (the last commit
touching src, applications or tutorials, so it is the suite key). Nothing is running.

## Where things stand

- **Tier 1 passes at 8d416fd.** Every check against exact or reference data passes, and every case
  agrees with its stored baseline. Baselines: V5, V12, V13, G7 from 8d416fd (user, 2026-11-14,
  `baa77b4`, `0b5a787`), all others from e7d742a. Results in `validation/suite/results/8d416fd/`
  (untracked), cases in `~/runs/suite/8d416fd/<ID>`.
- **Consolidated report:** `validation/suite/results/8d416fd/report/report.{tex,pdf}` (32 pages;
  figures and tables by `latexreport.py 8d416fd`). The method document `doc/cutCellIBM.{tex,pdf}`
  describes the ddtCorr and coupled-face flux change (1994dc8).
- **Last code change, 8d416fd** (finding 11 in `doc/validationSuite.md`): the transient Rhie–Chow
  correction carries the old fluid fraction, `interpolate(rAU*alpha0)` with alpha0 the group mean
  in a linked group, and `ibm.syncCoupledFlux(phi)` makes the processor and cyclic face flux
  single-valued after every pressure correction (ibmPimpleFoam, ibmSimpleFoam). It removed a
  decomposition dependence of V12 α = 2, D/h = 40 (serial, old and link-constrained decomposition
  now agree to 5 digits over t = 100–150). Transient results moved 1–2 % at D/h = 40, below 1 %
  finer, all toward the references. Evidence runs: `~/runs/v12a2` (E*, F*, X*, L*, D*, N1/N3/N5).
- **linked_solvers** is merged (12178f8); its handout `doc/linkedSolversRemergeHandout.md` is
  closed (the V12 D/h 160 rerun is done, no pressure solve at maxIter).
- The worktree `~/Codes/openFoam.ibm.docs/ddtcorr` (branch `ddtcorr_alpha`, = 8d416fd) is merged
  and no longer needed; the `xlinks*` worktrees listed in the linked_solvers handout likewise. Ask the
  user before removing any of them.

## Next: tier 2 (`doc/validationSuite.md` §7 step 3)

Order fixed by the plan: the prerequisites C2 and the zero-shear free-slip wall first, each accepted
by its gates; then V14 before V10, V10 and V11 before V6, the V7 sweeps, V15 complete and V16
before V8, V9, V12 complete. Decision D5 (junction pressure before V14–V16) is done (finding 2).

### 1. C2: tangential wall velocity on a static geometry — DONE (8108c14, 2026-11-14)

Implemented and accepted as below; results in `doc/validationSuite.md` §8 (C2 entry), scripts
and runs in `~/runs/c2` (`scripts/c2_tc.py`, `c2_pipe.py`, `c2_sphcouette.py`, `c2_v12.py`; set
`IBM_TREE` to run them on another worktree). User decision 2026-11-14: spherical Couette
(exact torque) accepts the rotating part of C2 instead of the periodic-box MIG-G2, which stays a
V6 gate for C3. Side check, not adopted: the wall point of u_b is clipped to the cell's bounding
box, also for linked cells whose value point (group centroid) lies outside the cell (734 of 6018
wall cells at d/h = 16); removing the clip in second order (branch `wallpoint_fix`, worktree
`~/Codes/openFoam.ibm.docs/wallpoint`) changes TC only at N = 160 (4th digit) and makes the
sphere errors 1–3 % larger, so the code stays. Worktrees `c2run` (= 8108c14, used for the
D/h = 160 run on istmhydra) and `wallpoint` can go once the user agrees.
Still open for V6/V12 tier 2: MIG-G1 into the suite (C3), V12 α = 2.5, 3, 4 steady.

The text below is the plan as written before the implementation.

Spec: `doc/validationSuite.md` §8, "Prerequisite code changes". A body whose motion leaves it
invariant (cylinder or sphere rotating about its axis, body of revolution about its axis, cylinder
sliding along its axis) keeps its geometry fixed, has no wall flux and only supplies a wall velocity
to the no-slip traction, so ibmSimpleFoam can run it.

**User decision (2026-11-14): explicit opt-in key**, valid with any motion type:
`motion { type none; angularVelocity (0 0 2); invariant yes; }`. Checked at start-up by a small
virtual motion of the level set (at the cut points); fatal if the shape changes, and fatal for an
alpha-field body (no level set). Existing cases stay as they are: `rotatingCylinder` and G5 keep
testing the rotation sources of a moving body; V12 α ≥ 2 and MIG opt in.

What the code does today (read, not changed):
- `bodyMotion::moving()` (`src/cutCellIBM/ibmBody/bodyMotion.{H,C}`): `noMotion` is moving when
  |ω| > 0; every other type always.
- `cutCellGeometry` (`cutCellGeometry.C` ~l. 265–315): `moving_` = any body moving; only then is
  `UbPtr_` (field `ibmUb`) allocated, the mass source allocated, the geometry rebuilt every step and
  written with the solution.
- `moving()` is used for two different things, which C2 must separate:
  (a) *the geometry changes* (rebuild, mass source, old face fractions, ddtCorr rescaling,
  uncovered cells, `solidOld_`, scalar transport l. 187, ibmContinuityErrs, pEqn mass source);
  (b) *a wall velocity exists* (`Ub()`, the no-slip and blanking sources `movingWallSource`, the
  viscous force `U - Ub` in `calcForces` ~l. 3087, the wall functions
  `cutCellGeometryWallFunction.C` l. 76, 572, 636, 673, SA `ibmSpalartAllmaras.C` l. 69, the
  solvers' UEqn.H, `ibmTruncationTC`).
  Suggested: a second flag (e.g. `wallVelocity_`) for (b), true for moving or invariant bodies;
  `moving_` only for (a). Ub of an invariant body is evaluated at every time (the motion may be
  time-dependent, e.g. oscillating along the axis), without rebuilding the geometry.
- ibmSimpleFoam refuses `ibm.moving()` (`createFields.H` l. 32); with C2 it must accept invariant
  bodies (they are not moving) and pick up Ub in UEqn.H.

Acceptance (spec): MIG-G1 (empty pipe with a sliding IBM wall, exact profile) and MIG-G2 (rotating
sphere torque 8πμR³Ω) are the C2 gates in §8; in 2-D, V12 α ≥ 2 steady with ibmSimpleFoam against
the ibmPimpleFoam result at t = 150 (α = 2: C_L −5.359, C_D 0.2915, C_M 0.4645 at D/h 160) and the
body-fitted values. Also: an invariant body run with ibmPimpleFoam must equal the moving-body run of
the same case to discretisation level, and G1–G9 must stay bit-identical (no case opts in).

### 2. Zero-shear free slip on curved walls

Plan: `doc/freeSlipPlan.md` (diagnosis, reference experiment in `~/runs/ibmSA_study/freeslip_fix`,
fixes 1–10). Behind a switch first: it changes every free-slip result by about 30 %. Acceptance:
V2 free slip against the body-fitted stress-free wall (4.003 at M = 128; now −24.8 %), then V11
BUB in 3-D. G7 and V2 baselines will change (user approval).

### 3. Tier-2 cases (§4, specs in §8, §10, §11)

New modules: V14 (DFG 3D), V10 (SPH), V11 (BUB), V6 (MIG, needs C2, the free-slip fix and the
driver C3), V16 (TMR bump), V8 (NACA), V9 (AMR). Extensions of tier-1 modules: V7 sweeps (cut
fraction, tilt, wall function), V15 complete (cut-fraction sweep, wall function), V12 complete
(α = 0.5, 1.5, 1.9, 2.5, 3, 4; steady ones with C2), V5 D/h 160 for 2D-2.

## Rules (user, still in force)

- Never `rm` case or result directories (mv is fine); scripts use `cd <dir> || exit 1`.
- Commit per step with the attribution lines; push only when the user asks.
- Never add the untracked `foulingmodel_ibmdavide/`; `validation/suite/results/` stays untracked.
- Do not remove first order. Numerical correctness first; no case-specific fixes; component-wise
  validation (`doc/developerNotes.md`).
- No baseline accepted or re-accepted without the user. Keep the link constraint; body-fitted cases
  never get the linkedGroup wrapper; nVcycles stays 2.
- Ask the user decisions one at a time, with context.
- Machines istmio2, istmhydra, istmeridanus, istmcetus (shared, keep load moderate), istmdraco;
  shared NFS. One build per worktree (`~/Codes/tools/buildtree_of.sh <tree> v2412`); never rebuild a
  tree whose binaries a running job uses; remote jobs detached (`ssh host 'cd ... || exit 1;
  setsid nohup ... &'`, istmio2 rebooted once and killed an attached run). Stop runs by session or
  PID, never a machine-wide pkill.

## Commands

```
cd ~/Codes/openFoam.ibm/validation/suite
~/.venv/ibm/bin/python suite.py run gates --jobs 24            # every code change
~/.venv/ibm/bin/python suite.py run tier1 --jobs 24
~/.venv/ibm/bin/python suite.py run V12 --reuse ~/runs/suite/<key>/V12   # re-evaluate
~/.venv/ibm/bin/python suite.py report --commit <key>
~/.venv/ibm/bin/python latexreport.py <key>   # then latexmk -pdf in results/<key>/report
```
