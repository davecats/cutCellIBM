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
