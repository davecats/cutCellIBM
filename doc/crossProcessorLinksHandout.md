# Handout: linked cells across processor faces slow the pressure solve (written 2026-11-11)

For a forked session. Goal: find out why the linear solvers converge much worse in parallel when
a linked group (slave + master, `mergeAlpha` > 0) straddles a processor face, and fix it in a way
that keeps the result independent of the decomposition. The converged solution is not the
problem (scotch on 8 ranks equals serial to 1e-10); the convergence of the linear solve is.

## Context (what the main session did and is doing)

- Branch `linked_solvers` (worktree `~/Codes/openFoam.ibm.docs/linksolve`, from
  `feature_iddes_from_scalar_secondOrder` at 2b1adc3). Commit 5af1381: the immersed steady cases
  use, for p, CG preconditioned by GAMG with the default (face-area pair) agglomeration and the
  `DICGaussSeidel` smoother; U and nuTilda use `smoothSolver` with Gauss-Seidel; linked or not.
- Commit 155ea8d adds a decomposition constraint: `ibmSetGeometry -linkZone` writes the faceZone
  `ibmLinkFaces` of the link faces and `decomposeParDict` gets `preserveFaceZones`, so
  `decomposePar` keeps each slave with its master (V5 D/h 160: 3 split groups -> 0); used by V5,
  V12 and the NACA generator. Start your worktree from `linked_solvers` at or after 155ea8d. That helps static bodies decomposed once; it does not help moving
  bodies (links change every step), adaptive refinement with load balancing, or any run whose
  decomposition is not under our control. Do not duplicate it; your target is the remaining
  cross-processor links.
- Tier-1 validation of 2b1adc3 is running from the main tree (`~/Codes/openFoam.ibm`, v2412
  build) on istmcetus, istmdraco and istmhydra until about 2026-11-12. Do not rebuild or check
  out anything in the main tree, and do not edit its `validation/suite` (the suite imports case
  modules lazily). Work in your own worktree and build (see Rules).

## Evidence

V5 DFG 2D-1 (Re 20), D/h = 80 (`re20_k1`, 126k cells, 160 linked cells), ibmSimpleFoam,
8 ranks, istmeridanus, binaries of the main tree (2b1adc3). Pressure: tolerance 1e-12,
relTol 0.01, maxIter 1000.

| p solver | scotch (0 links across) | manual cut (12 links across) |
|---|---|---|
| PCG + GAMG(faceAreaPair, DICGaussSeidel) precond. | 6.5 it/solve (max 24), 0.150 s/it, converged in 5664 SIMPLE it. | 46 it/solve (max 219), 1.35 s/it |
| PCG + GAMG(algebraicPair, DICGaussSeidel) precond. | 5.9 (max 24), 0.153 s/it | 116 (max 352), 3.16 s/it |
| GAMG(faceAreaPair, DICGaussSeidel) | 32 (max 668), 0.176 s/it | 772 (stalls at maxIter 1000), 5.7 s/it |

- Cases: `~/runs/linksolver3/<solver>` (scotch) and `~/runs/linksolver3x/<solver>` (method
  manual, `constant/cellDecomposition`: processor = 4 [x >= 0.15] + y band (4 bands); x = 0.15
  is a mesh face on the left side of the cylinder, where links are horizontal, so slave and
  master fall on different ranks). The 3x runs were stopped after 327-1375 iterations; their
  logs show the per-solve iteration counts. Run script `run.sh` in each directory.
- Serial reference: `~/runs/suite/31d1827/V5/re20_k1` (C_D force total x 5.5808943921e-05,
  5660 iterations); scotch parallel 5.58089439196e-05.
- Cross-processor links also occur in the validation suite without any manual cut:
  V5 2D-1 D/h = 160 (`re20_k2`, scotch, 8 ranks) has 3 (its DIC-PCG hit maxIter 1000, see
  `v5.py` comment, which is why it uses the GAMG preconditioner); V13 (oscillating cylinder,
  8 ranks) has 1-20 over time (log lines "cell linking: ... (N over coupled faces)").
- Serial and in-processor experiments (same day) showed the incomplete-Cholesky smoother is
  what resolves a linked pair; plain Gauss-Seidel GAMG stalls on linked pressure even in serial.
  Notes: `~/.claude/.../memory/project-linking-solver.md` (if available to you) and the commit
  message of 5af1381.

## Mechanism (hypothesis to verify first)

`cutCellGeometry::linkCells` (`src/cutCellIBM/cutCellGeometry/cutCellGeometryTemplates.C`)
adds a penalty conductance c = mergePenalty (1e4) x max(|diag|) between slave and master. For a
link over a coupled face it goes into `internalCoeffs` (diagonal) and `boundaryCoeffs` (the
off-diagonal through the processor interface). OpenFOAM's smoothers and preconditioners treat
processor interfaces explicitly: Gauss-Seidel and DIC update the interface contribution from
the previous sweep (block Jacobi across ranks), and GAMG agglomerates only within a rank (also
on the coarse levels, until processor agglomeration, if any). For a pair coupled by c >> the
other coefficients, a Jacobi step leaves the pair-mean mode almost unchanged (factor about
c/(d + c) = 1 - 1e-4), and no coarse level ever merges the pair. Check this: e.g. count the
iterations against mergePenalty (1e4, 1e2, 1e1) on the 3x case, and look at where the residual
sits after a solve (the cross-processor linked cells?).

## Directions (ranked by my expectation; discuss with the user before building a big one)

1. Group correction across ranks: a preconditioner (or wrapper around the existing one) that,
   after each application, corrects each linked group by its merged (summed) equation: restrict
   the residual to the group sum, solve the 1x1 merged equation per group, add the correction to
   all group cells. The cross-rank group sums already exist (`groupSum`, send/recv slots). This
   is a deflation / two-level additive correction on exactly the slow mode. Runtime-selectable
   `lduMatrix::preconditioner` in the IBM library.
2. Penalty-free coupling of cross-rank links: keep a moderate conductance on coupled links only
   and make the converged pressure group-uniform by a deferred correction (as linkCells already
   does for momentum with `groupUniform`). Risk: per-step continuity in transient PISO/PIMPLE.
3. GAMG processor agglomeration (`processorAgglomerator`, e.g. masterCoarsest) to see whether
   merging ranks on coarse levels helps at all (config only, quick test).
4. Physical merging of cells (removeFaces before decomposePar) or a master choice that avoids
   rank boundaries: both break or bypass decomposition independence or only cover static bodies.

Acceptance: on the 3x case, p iterations per solve within ~1.5x of the scotch case; results
equal across decompositions to the G8 tolerance (serial vs 4 and 6 ranks); G1-G9, V1, V5, V13
pass in your worktree build; no change to serial results beyond round-off unless argued.

## Rules of this project (from the user)

- Never `rm` case or result directories; `mv` is fine. Never add `foulingmodel_ibmdavide/`.
  Never push. Commit per step on your own branch, ending messages with
  `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` and
  `Claude-Session: https://claude.ai/code/session_015RDCGhLjssV21QpDg7kQrE`.
- Numerical correctness first; no case-specific fixes; validate component-wise
  (`doc/developerNotes.md`). Ask the user decisions one at a time, with context.
- Do not accept or re-accept a suite baseline without the user's explicit approval.
- One build per worktree: `git worktree add -b <branch> ~/Codes/openFoam.ibm.docs/<name>
  linked_solvers`, build with `~/Codes/tools/buildtree_of.sh <tree> v2412` (refuses to write
  outside the tree). Never rebuild a tree whose binaries a running job uses.
- Machines: istmio2 (24 cores), istmeridanus (64), istmhydra (256), istmcetus (128, shared with
  ~48 users, keep load moderate), istmdraco (256); same NFS home and /opt/OpenFOAM. Check load
  first; the main session's tier-1 runs occupy cetus, draco and hydra. Launch remote jobs with
  `~/Codes/tools/launch.sh <log> ssh <host> <script>`. Scripts must `cd <dir> || exit 1`; stop
  runs by their session (`pkill -s <SID>` of the run script), never `pkill -x <solver>`
  machine-wide.
- OpenFOAM v2412 environment: `source /opt/OpenFOAM/OpenFOAM-v2412/etc/bashrc`, ThirdParty lib
  `/opt/OpenFOAM/ThirdParty-v2412/platforms/linux64GccDPInt32/lib` on LD_LIBRARY_PATH for
  decomposePar, then `source <tree>/etc/bashrc`. The code must keep building on v2106-v2506.

## Result (2026-11-11, branch `cross_proc_links`, worktree `~/Codes/openFoam.ibm.docs/xlinks`)

**Mechanism: confirmed, with one addition.** On the 3x case (first 150 SIMPLE iterations,
PCG + GAMG(faceAreaPair, DICGaussSeidel), 8 ranks) the pressure iterations per solve follow
the penalty: 36.3, 17.9, 13.8, 10.5 for `mergePenalty` 1e4, 1e2, 1e1, 1 (scotch, no split
group: 12.7). The addition: with GAMG the slow mode exists on every level, because the two
clusters that hold the cells of a split link are again coupled by c. A correction of the
group means on the fine level only (direction 1 as written) therefore repairs a
single-level preconditioner (DIC-PCG: 744 -> 548, scotch 522) but not GAMG (36 -> 50).

**Fix: the preconditioner `linkedGroup`** (`src/cutCellIBM/linkedGroupPreconditioner`, the
class header has the derivation), a wrapper around any base preconditioner:

    preconditioner { preconditioner linkedGroup; base { preconditioner GAMG; smoother DICGaussSeidel; ... } }
    preconditioner { preconditioner linkedGroup; base DIC; }       // or DILU with PBiCGStab

1. The base preconditioner gets a copy of the matrix in which the links over coupled faces
   are softened to 1 x diagonal (`softPenalty`), so no level sees a stiff coupling across
   processors. Fields uniform in a split group do not feel the conductance, so for them the
   softened matrix is the right one.
2. The residual of a split group enters the base as that of its merged equation (the sum
   over its cells, handed to them in proportion to their coefficients to the outside); the
   result is averaged over the group with the same weights. With these weights the softened
   conductance no longer matters (16.7 iterations for `softPenalty` 0.01 to 1 at that stage).
3. What is left, the differences inside a split group, takes half a Jacobi step (exact for a
   pair: its eigenvalue of D^-1 A is 2). A full step was the last 25 % (16.7 -> 13.0).

It is symmetric if the base is, changes the path of the linear solver and not the equations,
and is the base preconditioner when no group is split. `cutCellGeometry` now keeps the cells
of the split groups (`splitGroupCells`, also over cyclic faces) and exposes the exchange of
the group sums (`groupSlotSum`).

**Numbers** (V5 D/h 80, 8 ranks, the manual cut of `~/runs/linksolver3x`; cases, scripts and
logs in `~/runs/xlinks`):

| decomposition, p preconditioner | split links | p iterations per solve | SIMPLE iterations, C_D |
|---|---|---|---|
| first 150 iterations: | | | |
| scotch, GAMG | 0 | 12.71 | |
| scotch, linkedGroup(GAMG) | 0 | 12.71 (identical log) | |
| manual cut shifted by one cell (x >= 0.14875), GAMG | 0 | 12.99 | |
| manual cut, GAMG | 12 | 36.3 (max 219) | |
| manual cut, linkedGroup(GAMG) | 12 | 12.99 (max 23) | |
| manual cut at 0.15125 / 0.1525, GAMG | 6 / 4 | 33.8 / 27.7 | |
| manual cut, DIC (40 iterations) | 12 | 744 (max 1000) | |
| manual cut, linkedGroup(DIC) | 12 | 529 | |
| scotch, DIC | 0 | 522 | |
| to convergence (residualControl 1e-9): | | | |
| manual cut shifted, GAMG | 0 | 6.66 | 5663, 5.58089439221e-05 |
| manual cut, linkedGroup(GAMG) | 12 | 6.66 | 5664, 5.58089439158e-05 |
| manual cut, GAMG (stopped at 1398) | 12 | 46.7 | |
| serial reference / scotch (handout above) | | 6.5 | 5660 / 5664, 5.5808943921e-05 / ...196e-05 |

Acceptance: 1.00x the iterations of the same decomposition without split links (target
1.5x of scotch: 1.02x); C_D equal to the serial value to 9e-11.

Moving body, V13 D/h 40 (scotch, 8 ranks, 0 to 10 split links during the motion; first two
periods), iterations per solve by the number of split links:

| split links | DIC-PCG p: before | with linkedGroup(DIC) | PBiCGStab-DILU Ux: before | with linkedGroup(DILU) |
|---|---|---|---|---|
| 0 | 204 | 204 | 3.1 | 3.1 |
| 1 | 276 | 206 | 8.5 | 4.3 |
| 3 | 333 | 207 | 30 | 4.6 |
| 5 | 407 | 206 | 44 | 4.7 |
| 7 | 453 (max 917) | 206 | 51 | 4.8 |
| 10 | 535 (max 1000 = maxIter) | 206 (max 417) | 57 | 4.9 |

(In the stored V13 runs, e.g. `~/runs/suite/6692870/V13/Dh80`, the pressure solve hit
maxIter 1000 at 3 or more split links, i.e. stopped short of its tolerance; `v13corr.py` in
`~/runs/xlinks` makes this table from a log.)

**Cost.** No measurable cost per iteration: 288 ms per SIMPLE iteration with 12 split links
and the wrapper against 293 ms for the shifted cut without split links, at equal iteration
counts (istmdraco, runs alternated). The group sums use PstreamBuffers; fixed-size
non-blocking messages instead gave nothing measurable (8 ranks 288 against 289 ms, 32 ranks
462 against 457 ms) and stay on the branch `cross_proc_links_rawexchange`.

**Decomposition independence and gates** (v2412 build of 8a40eb8, results in
`validation/suite/results/8a40eb8`): G1 to G9 and V1 give the values of 2b1adc3 to
round-off (G1 to G8 <= 2e-15; G9 2e-9 in one observed order); the two regression failures
of G6 and G7 are those of 2b1adc3 with the same digits. G8 did not contain a split group
(scotch on 4 and 6 ranks keeps the 36 linked cells of cylinderRe1 with their masters; the
wall-resolved channel has no linked cell): it now also runs cylinderRe1 on 5 slabs along x,
which split 6 links, and requires at least one: drag equal to serial to 1.4e-12. V13 uses
the wrapper for p and U, the cylinderRe1 tutorials for p. The code builds and gives the same
V13 iteration counts on v2106 and v2506.

**V13 and V5 of the suite** (8a40eb8, on istmeridanus): both PASS. V13 with the wrapper
for p and U: c_d and c_i at D/h = 160 within 4e-7 and 8e-7 of the run of 2b1adc3, regression
against the baseline at 0.04 of its tolerance. Iterations per solve at D/h = 160 for 1 to 20
split links: DIC-PCG 545 to 576 (613 to 737 in the run of 2b1adc3; 550 with none),
PBiCGStab-DILU for Ux 6.4 to 7.6 (13 to 153; 6.3 with none); at D/h = 80 324 to 365 (340
with none) and 5.5 to 7.0 (4.7). At D/h = 160 the DIC-PCG pressure solve reaches maxIter
1000 also without any split link (mean 550): the plain DIC preconditioner on that mesh, not
the links. V5 does not use the wrapper (decomposition constraint): PASS, values within 5e-7
of the run of 2b1adc3 (the solver settings of 5af1381), regression at 0.008 of its tolerance.

**Not covered.**
- Solvers without a preconditioner cannot use it: GAMG as the solver (the AMR tutorial) and
  smoothSolver. For the pressure that means slow convergence with split groups. For U with
  smoothSolver the iteration count does not react to split links (8.8 sweeps with and
  without), although a Gauss-Seidel sweep cannot move the mean of a split group any better
  than the smoothers of the pressure: the solver stops on the global residual. What this
  leaves behind was not measured. It does not show in the steady V5 run (5664 SIMPLE
  iterations with 12 split links against 5663 without), and for U the pressure corrector
  rebuilds the velocity from A and H of the unlinked equation, so the predictor matters
  little; a transported quantity without a corrector (nuTilda, a scalar) solved by
  smoothSolver in a time-accurate run would be the case to check. PBiCGStab with
  linkedGroup(DILU), as in V13, avoids the question.
- V5, V12, V15 and NACA keep the decomposition constraint and their plain preconditioners
  (nothing is split there). The wrapper would make the constraint unnecessary; that is a
  decision for the user.
