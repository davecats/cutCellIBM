# Handout: from the tier-1 baseline at 8d416fd to tier 2 (written 2026-11-14)

Branch `feature_iddes_from_scalar_secondOrder`, pushed to origin at 690bd53 (this handout comes
after). Main tree `~/Codes/openFoam.ibm`, built with OpenFOAM v2412 at 2190505 (2026-11-15; the last
commit touching src, applications or tutorials, so it is the suite key). Nothing is running.

## Where things stand

- **Tier 1 passes at 8d416fd.** Every check against exact or reference data passes, and every case
  agrees with its stored baseline. Baselines: V5, V12, V13 from 8d416fd (user, 2026-11-14,
  `baa77b4`, `0b5a787`), G7 and V2 from 2190505 (zero shear on, user 2026-11-15, `2194276`), all
  others from e7d742a. Results in `validation/suite/results/8d416fd/`
  (untracked), cases in `~/runs/suite/8d416fd/<ID>`.
- **Consolidated report:** `validation/suite/results/8d416fd/report/report.{tex,pdf}` (35 pages;
  figures and tables by `latexreport.py 8d416fd`; the V2 figures and table and the G7 force table
  of 2190505 copied in as `V2zs_*`, `V2zs.tex`, `G7zs_forces.tex` from `latexreport.py 2190505`). The method document `doc/cutCellIBM.{tex,pdf}`
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

Order fixed by the plan: the prerequisites C2 and the zero-shear free-slip wall first (both done; zero
shear the default since 2190505, G7 and V2 re-accepted from it), each accepted by its gates; then V14 before V10, V10 and V11 before V6, the V7 sweeps, V15 complete and V16
before V8, V9, V12 complete. Decision D5 (junction pressure before V14–V16) is done (finding 2).

### 1. C2: tangential wall velocity on a static geometry — DONE (8108c14, 2026-11-14)

Implemented and accepted as below; results in `doc/validationSuite.md` §8 (C2 entry), scripts
and runs in `~/runs/c2` (`scripts/c2_tc.py`, `c2_pipe.py`, `c2_sphcouette.py`, `c2_v12.py`; set
`IBM_TREE` to run them on another worktree). User decision 2026-11-14: spherical Couette
(exact torque) accepts the rotating part of C2 instead of the periodic-box MIG-G2, which stays a
V6 gate for C3. Side check, not adopted: the wall point of u_b is clipped to the cell's bounding
box, also for linked cells whose value point (group centroid) lies outside the cell (734 of 6018
wall cells at d/h = 16); removing the clip in second order (branch `wallpoint_fix`, commit fff0589)
changes TC only at N = 160 (4th digit) and makes the sphere errors 1–3 % larger, so the code
stays. The worktrees `c2run` and `wallpoint` were removed (user, 2026-11-14). The consolidated
report `validation/suite/results/8d416fd/report/report.{tex,pdf}` has a C2 section (§19) and
the C2 entries in its summary.
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

### 2. Zero-shear free slip on curved walls — DONE, `zeroShear` default on (2026-11-14/15)

Commits f4f2366 (switch), 361c23c (mergeH, finding 12), 8d60565 (normal viscous flux), 2190505
(default on, user 2026-11-15), 5272abf (fix 9: tutorial, notebook comparison on a copy with
`zeroShear false`, `cutCellIBM.tex` §9.5, §11.6, README), 2194276 (G7 and V2 baselines
re-accepted from 2190505, user 2026-11-15). Status,
decisions and numbers: `doc/freeSlipPlan.md` "Status"; `doc/validationSuite.md` §8 prerequisites,
findings 12 and 13; report `validation/suite/results/8d416fd/report` §ZS. Runs and scripts:
`~/runs/zeroshear` (`scripts/v2zs.py`, `tc_zs.py`, `sph_zs.py`, `mms_zs.py`, `cmp_metrics.py`).
Gates at f4f2366, 361c23c, 8d60565: PASS; all in-plane metrics bit-identical to 8d416fd. At
2190505 (default on): G1–G6, G8, G9 bit-identical; G7 changes only through the free-slip bodies
(+37 %, +39 %, the no-slip neighbour −2.8 %), V2 free slip 3.9990 at N = 321 (−0.105 %, PASS within
0.3 % of the body-fitted stress-free 4.00318); results in `validation/suite/results/2190505/`
(gates, V2), cases in `~/runs/suite/2190505`.

User decisions (2026-11-14): global key `zeroShear` (default false); curvature from the level-set
Hessian (triSurface/alpha-field fall back with a warning); Robin factor (I + d H)⁻¹ for the wall
value, kept after the V2 study ("keep Robin, investigate V2": the V2 error was the missing normal
viscous flux); U_z growth fixed at the source (mergeH), not by zeroing U_z.

Results: exact rigid rotations (Taylor–Couette convex and concave, transient and steady; 3-D
spherical shell) and a manufactured solution with a varying slip converge at second order; V2
free slip linked 3.9641, 3.9895, 3.9974, 3.9990 (N = 41 … 321) against the body-fitted
stress-free limit ≈ 3.999 (−24.8 % before).

Next:
- Finding 13 (unlinked free slip with alphaMin > 0: absorbed slivers block the slip flow).
- V14 PASS (2026-11-17, finding 18: D/h 46.8 on HoreKa, c_L,min judged on the extrapolation, user decision; baseline accepted 2026-11-17).
- V11 BUB done (2026-11-17, finding 17, PASS with C_D judged on the extrapolation, user decision; baseline accepted 2026-11-17); MIG-G2 stress-free sphere in V6.
- V6 MIG done (2026-11-20, finding 21): driver `mig.py`, module `v6.py`; gates G1, G2 PASS; cases
  R1–R4, S1, S2 ran on HoreKa with 2190505 (9 single-node jobs, driver `~/runs/v6/scripts/v6_hk.py`),
  the Re 128 triplet (d/h 32, 48, 64 at eps/eps* 0.4) and its mergeAlpha 0 point with 8c6cc10
  (finding 19). β is judged on figure 4c (the paper's fit is ~1.23× below it). At eps/eps* 0.9 the
  paper's f/Re (figure 9 and its fit) does not match its own V, Ω, β, which we match; L = 4 gives
  the same f/Re as L = 3 (HoreKa job 5253243, results `~/runs/suite/8c6cc10_horeka/V6_L4`), so f/Re
  at 0.9 is REPORT (user, e9746a3). `suite.py run V6 --reuse ~/runs/suite/2190505/V6` evaluates PASS
  (key 8c6cc10, post-processing only: every point has its result.json; full fields only for the
  figure case). Report §V6 written; baseline accepted from 8c6cc10 (user, 2026-11-20). MIG cache check:
  a cached solve must be the (V, Ω) asked for (d85fb9c).
- ibmConvection fix (finding 19 in the branch's `doc/validationSuite.md`): branch
  `convection_upwind`, worktree `openFoam.ibm.docs/convup` (code 8c6cc10, docs to fe86bec), HoreKa
  clone `<ws>/convup` (build `build_convup.sh`, tier 1 `tier1_convup.sh`). Implicit upwind on the
  corrected faces, the value-point face value as the explicit correction phiC (e − v_U), −Sp(div
  phiC) for a bounded scheme. Gates: G7 and G9 fail only their regression (Rhie–Chow rAU = α_U/A
  shift, ∝ α_U). Tier 1 on HoreKa: G7, G9, V2, V5, V12, V13 fail only their regression (shifts
  falling tenfold per refinement), the rest PASS; merged 34f4e3d and the six baselines re-accepted
  from 8c6cc10 (dc26911, user, 2026-11-19). The α_U dependence of steady solutions (also in the old code, 1.6–2.6 % max|U| at
  G9 N 20) stays a documented property for now (user, 2026-11-17); remedy if wanted later:
  relaxation-independent Rhie–Chow.
- V10 done (2026-11-18, finding 20): the missing Re 150 D/h 64 ran on HoreKa (job 5246721, 17 min);
  the failures against the digitised J&P figures were the references (Fornberg 1988 tabulated);
  references amended (user), PASS, baseline accepted from 2190505 (10b73ee), report §V10.
- Worktrees made this session (ask before removing): `openFoam.ibm.docs/zsdiag` (branch
  zeroshear_diag, diagnostic keys zsRobin, zsNormalFlux), `zsdiag2` (zeroshear_diag2,
  zsImplicit, zsLagTime), `zsroot` (detached, 361c23c + zeroShearNormalFlux/Mode). N = 641 V2
  runs: `~/runs/zeroshear/n641` (Robin 3.9932, no Robin 3.9980 at 40000 iterations, without
  the normal flux).

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
- At the end of every session, update the validation documents if anything changed: the
  consolidated report (`validation/suite/results/<key>/report/report.tex`, hand-written sections
  and the summary of findings, decisions, open items and next; rebuild the PDF and copy it to
  `doc/validationReport.pdf`, linked from the README with its status line) and
  `doc/validationSuite.md` (user, 2026-11-14).
- Machines istmio2, istmhydra, istmeridanus, istmcetus (shared, keep load moderate), istmdraco;
  shared NFS. One build per worktree (`~/Codes/tools/buildtree_of.sh <tree> v2412`); never rebuild a
  tree whose binaries a running job uses; remote jobs detached (`ssh host 'cd ... || exit 1;
  setsid nohup ... &'`, istmio2 rebooted once and killed an attached run). Stop runs by session or
  PID, never a machine-wide pkill.
- Runs expected to take longer than 4 h go to HoreKa 2 (`hk2x86`, user 2026-11-17 / 2026-11-20), not
  the workstations. Legacy HoreKa ended 2026-11-20; its workspace `xt8786-ibm` (342 GB) is copied to
  `/hfs2/work/workspace/scratch/ka_xt8786-ibm`. Code `$PROJECT/cutCellIBM` (clone of a git bundle), built
  against `OpenFOAM/v2506-foss-2025b` (no v2412 module); jobs `$PROJECT/ibmjobs/` (`build.sh`, `gates.sh`:
  `module load OpenFOAM/v2506-foss-2025b`, `IBM_SUITE_FOAM="source $FOAM_BASH"`, the module's `python3`),
  suite runs `ka_xt8786-ibm/suite_hk2`. Gates G1–G9 PASS there (2026-11-20, job 137684). Notes `~/HOREKA.md`.
- Every case section of the report has a description of the flow case and a picture of the flow
  with the immersed wall (`flow3d.py` for the 3-D cases, `flowfield.py` for 2-D; user, 2026-11-17).

## Commands

```
cd ~/Codes/openFoam.ibm/validation/suite
~/.venv/ibm/bin/python suite.py run gates --jobs 24            # every code change
~/.venv/ibm/bin/python suite.py run tier1 --jobs 24
~/.venv/ibm/bin/python suite.py run V12 --reuse ~/runs/suite/<key>/V12   # re-evaluate
~/.venv/ibm/bin/python suite.py report --commit <key>
~/.venv/ibm/bin/python latexreport.py <key>   # then latexmk -pdf in results/<key>/report
```
