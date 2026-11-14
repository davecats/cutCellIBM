# Handout: merge `linked_solvers` again into `feature_iddes_from_scalar_secondOrder` (written 2026-11-12)

For the sessions that work on `feature_iddes_from_scalar_secondOrder`. That branch merged
`linked_solvers` at 316172f (merge 20260c5). Since then `linked_solvers` has received the work
of the forked session on linked groups split over processors (`cross_proc_links`) and four
follow-ups the user decided on 2026-11-12. This note says what is new, how to merge it, what
it changes in the results and what is left open. The derivations and all numbers of the first
part are in the Result section of `doc/crossProcessorLinksHandout.md`.

## What is new on `linked_solvers` since 316172f

| commit | content |
|---|---|
| b5fc7b6, 8a40eb8 | `linkedGroup` preconditioner (`src/cutCellIBM/linkedGroupPreconditioner`); `cutCellGeometry` keeps the cells of the split groups (`splitGroupCells`, `nSplitGroups`, `groupSlotSum`, `mergePenalty()`); V13 and the cylinderRe1 tutorials use it |
| ec8affd | G8 also runs cylinderRe1 on 5 slabs along x, which split 6 linked groups (the scotch runs split none) |
| 43cafca | the wrapper in V5, V12, V15 and the NACA generator (they keep the link constraint: nothing is split there, the wrapper is its base); AMR tutorial: p by PCG + linkedGroup(GAMG) instead of GAMG as the solver; V13: GAMG base |
| cccd578 | V13: two V-cycles in the GAMG base |
| b3dcfec | V13 baseline re-accepted from 43cafca (user, 2026-11-12) |
| 5d9caa5 | V5, V15, NACA: coarsest level of the GAMG preconditioner to relTol 0.01 instead of round-off (user, 2026-11-12) |
| e1c86bc ... 78e93ce, this file | documentation |
| merges of 8337bd9 and 935d9b5 | `feature_iddes_from_scalar_secondOrder` merged into the branch, so that it merges back without conflict |

The code (src, applications, tutorials) last changed in 5d9caa5 (the NACA generator; src in
8a40eb8): that is the key of the suite (`common.commit()`), also after the merge. Nothing the
suite runs differs between 43cafca and 5d9caa5 except the fvSolution of V5 (2D-1) and V15.

## How to merge

In the main tree (`~/Codes/openFoam.ibm`), when no job uses its binaries (none did on
2026-11-12 20:30; check with `readlink /proc/<pid>/exe` of the running solvers):

    git merge --no-ff linked_solvers          # no conflict as of 935d9b5 (tried with git merge-tree)
    ~/Codes/tools/buildtree_of.sh ~/Codes/openFoam.ibm v2412   # src changed: new class, cutCellGeometry

If the branch has moved on and `validation/suite/v12.py`, `v13.py`, `v5.py`, `v15.py`, `g8.py`
or `doc/validationSuite.md` conflict: keep the `linkedGroup` solver entries and the G8 row of
this branch. `fvSolution` entries of the immersed cases now read

    p { solver PCG; preconditioner { preconditioner linkedGroup; base DIC; } ... }                    // or
    p { solver PCG; preconditioner { preconditioner linkedGroup; base { preconditioner GAMG; ... } } ... }
    U { solver PBiCGStab; preconditioner { preconditioner linkedGroup; base DILU; } ... }

Body-fitted cases must not get the wrapper (it lives in libcutCellIBM). The build was checked
on v2106 and v2506 at 8a40eb8 (same V13 iteration counts); nothing in src changed afterwards.

## Results at 43cafca and 5d9caa5 (v2412, worktree `~/Codes/openFoam.ibm.docs/xlinks`)

JSON in `xlinks/validation/suite/results/43cafca` (and `8a40eb8`, `5d9caa5`), cases in
`~/runs/suite/43cafca` and `~/runs/suite/5d9caa5`. The code of the merged branch is the same, so these results are its
results; copy the JSON files if the report is to be made from the main tree.

| id | result | against 2b1adc3 / the baseline |
|---|---|---|
| G1 to G5, G8, G9, V1 | PASS | equal to 2b1adc3 to round-off (G9: 2e-9 in one observed order) |
| G6, G7 | PASS with the tolerances of 8337bd9 | same values as 2b1adc3 |
| G8 | PASS, with the new split-group run | drag of the 5-slab run equal to serial to 1.4e-12 |
| V15 | PASS | regression 0.009 of the tolerance |
| V13 | PASS (run `V13.1`, two V-cycles) | regression 0.53 of the tolerance 1e-4, see below |
| V5 | PASS at 8a40eb8, i.e. before it got the wrapper | regression 0.008; not rerun at 43cafca: no group is split, the wrapper is its base (the logs of the V5 D/h 80 case with and without it are identical) |
| V2, V3, V4, V7, V12 | not run on this branch | V12 see below |

`results/43cafca/V13.oneVcycle.json` is the V13 run with one V-cycle (cases in `V13`, not `V13.1`).

At 5d9caa5 (coarsest level to relTol 0.01 in V5 2D-1 and V15), started 2026-11-12 23:00 from
the xlinks worktree: G1 to G7 PASS when this was written; G8, G9, V1, V15 and V5 were still
running (logs `~/runs/xlinks/suite_*_5d9caa5.log`, JSON in `xlinks/validation/suite/results/5d9caa5`).
Checked directly: V5 2D-1 at D/h 80 (scotch, 8 ranks) converges in 5663 SIMPLE iterations with
6.49 pressure iterations per solve to C_D 5.58089439168e-05 (before: 5664, 6.5,
5.58089439196e-05; serial 5.5808943921e-05): the converged steady state does not depend on
the setting. V13 and the gates do not use the changed files.

### V13 changes in the fifth digit, and the baseline is the less accurate side

With the GAMG base c_d and c_i differ from the baseline (DIC-PCG) by 2e-5 to 5e-5:

| | c_d D/h 40 | c_d D/h 80 | c_d D/h 160 | c_i D/h 160 |
|---|---|---|---|---|
| baseline (DIC) | 2.1412144 | 2.0785776 | 2.0544565 | 1.4627526 |
| 43cafca (GAMG, two V-cycles) | 2.1411723 | 2.0786382 | 2.0545660 | 1.4627752 |

D/h 40 with the pressure solved to tolerance 1e-12 and relTol 1e-4: 2.1411846 (DIC) and
2.1411688 (GAMG). The GAMG value at the normal tolerances is within 6e-6 of both, the DIC value
2e-5 off: at tolerance 1e-9 and relTol 0.01 DIC-PCG leaves an error of that size, and at
D/h 160 a quarter of its solves stopped at maxIter 1000 on top of that. The check passes; **the
baseline was re-accepted from this run** with the user's approval (b3dcfec).
Iterations per solve 19, 24, 28 and 0.20, 0.58, 1.60 s per step at D/h 40, 80, 160
(istmeridanus, 8 ranks each, the three at once); no solve at the limit.

### V12 at D/h 160: the stored runs are not converged in the pressure, rerun needed

In the runs of 2b1adc3 and c8c771f (48 ranks, 6 linked groups split, DIC-PCG without the
wrapper) half of the pressure solves stop at maxIter 1000 with a residual of about 7.5e-7 for
the tolerance 1e-8. With the wrapper, same decomposition, 2000 steps: 268 iterations per
solve, none at the limit, 0.45 s per step (with a GAMG base and one V-cycle 20.5 iterations
and 0.38 s per step; two V-cycles were not tried there). In the merged state V12 has the wrapper and the link
constraint, so its solves converge; its D/h 160 numbers should come from a new run. V12
keeps the DIC base (the user did not ask for more); by the V13 experience a DIC result may
carry a solver error of order 1e-5. The rerun is with the feature-branch session (user).

**Done on the feature branch (2026-11-14, key 8d416fd).** Rerun with the wrapper and the link
constraint: none of the 1 008 000 pressure solves at maxIter (mean 158, largest 978 iterations);
alpha 2, D/h 160: C_L -5.359, C_D 0.2915, C_M 0.4645, all checks PASS. The run of 5d9caa5
exposed a decomposition dependence at D/h 40 (processor-face flux amplified by ddtCorr), fixed
in 8d416fd; V5, V12, V13 and G7 re-accepted from 8d416fd by the user (doc/validationSuite.md
section 12).

## Findings and what the user decided (2026-11-12)

1. **Coarsest level of the GAMG preconditioner: changed (5d9caa5).** V5, V15 and NACA used
   `tolerance 0; relTol 0` in the preconditioner dictionary, which solves the coarsest level
   to round-off in every application. V5 D/h 80, scotch, 8 ranks: `relTol 0.01` gives the same
   iterations (12.57 against 12.56 per solve) at 96 instead of 168 ms per SIMPLE iteration.
   V13 has it too. Results change at solver-tolerance level (see the 5d9caa5 runs below).
2. **One V-cycle is not safe** (`nVcycles 1`): faster per step, but one pressure solve per
   period of V13 at D/h 80 stalled (cccd578). Keep the default of two.
3. **The link constraint stays** in V5, V12 and NACA (user). The wrapper would handle split
   groups there (V5 D/h 80 with 12 split links: 12.99 iterations per solve, as without any).
4. **smoothSolver and GAMG as the solver cannot use the wrapper.** What smoothSolver leaves
   behind at a split group in a time-accurate run (nuTilda, scalars) is not measured; the AMR
   tutorial still solves U with smoothSolver.
5. **Fixed-size messages for the group sums** instead of PstreamBuffers gave no measurable
   gain (branch `cross_proc_links_rawexchange`).

## Rules that applied here (unchanged)

No baseline accepted without the user; no `rm` of case or result directories; one build per
worktree, never rebuild a tree whose binaries a running job uses; no push.

## Where things are

- Branches: `linked_solvers` (merged state), `cross_proc_links` (the same commits),
  `cross_proc_links_rawexchange` (finding 5).
- Worktrees that can go (`git worktree remove`): `~/Codes/openFoam.ibm.docs/xlinks_b`,
  `xlinks_v2106`, `xlinks_v2506`; `xlinks` holds the suite results named above.
- Test cases, scripts and logs: `~/runs/xlinks` (`stats.py`, `window.py`: iterations and time
  per iteration of a log; `v13corr.py`: iterations against the number of split links;
  `capstats.py`: solves at the iteration limit; `v13eval.py`: c_d, c_i of a V13 case;
  `mkcut.py`: the manual decompositions).
