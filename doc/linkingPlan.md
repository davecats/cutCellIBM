# Handout: linking always on, decomposition independent (written 2026-11-03)

Branch `feature_iddes_from_scalar_secondOrder`. Last commits: 0d76f2f (linking without sliver
absorption), 507da62 (cut-cell pressure = continuity error). Nothing pushed.

**Goal.** Make cell linking (mergeAlpha > 0) the default and replace α-based sliver absorption,
with links that cross processor (and cyclic) faces, so that the result does not depend on the
decomposition: same link groups, same fields to solver tolerance, serial or any number of
processors.

**Standing rules.** No `rm -rf` of cases or results (build fresh directories). Commit per step
on this branch, never push unless asked. Default results stay byte-identical unless a change is
deliberate and reported (Phase C breaks this on purpose, first order included). The untracked
`foulingmodel_ibmdavide/` is not ours, never add it.

**Environment.**
`source /etc/profile.d/modules.sh; module load cae/openfoam+/v2106_prf;
source /opt/OpenFOAM/OpenFOAM-v2106_prf/etc/bashrc; source ~/Codes/openFoam.ibm/etc/bashrc`,
build with `./Allwmake`. OpenFOAM sources: `/opt/OpenFOAM/OpenFOAM-v2106_prf/src`.

## Status (2026-11-03, executed)

Commits: 5baba93 (A), 7f50469 (B core), a3c5c06 (B6), 0fe162c (C1), 5300cfc (C3 docs), 8aebbce
(flow-rate forcing), e99057b (SA-neg), b81d779 (solid-row diagonal). Nothing pushed.

**Phase A** (5baba93): A1 scalar linked (`ibm_.linkCells(TEqn)` after `relax()`, before the
fvOptions constraints, as ν̃: the constraints must see the linked row); A2 fluid cells with no
open internal/coupled face blanked (none occurred in any test); A4 masters from the original
link graph; A5 tie-break by face centre. Regression: default cases byte-identical; linked TC
changes its tie choices (force asymmetry 1e-8 rotated by 90°).
A1 check, static-scalar tutorial, wall flux / injected (unlinked extrapolates to about 0.3995):

| | N81 | N161 | N321 |
|---|---|---|---|
| unlinked | 0.39203 | 0.39594 | 0.39771 |
| linked flow, scalar unlinked (old) | 0.38283 | 0.39370 | 0.39586 |
| linked flow and scalar (new) | 0.40121 | 0.39887 | 0.39958 |

A3 (TC second order, inner cut-cell pressure rms, Phase-A binary): unlinked 0.02 / 1e-3 / 1e-9:
N80 1.67e-2 / 1.67e-2 / 1.61e-2, N160 1.35e-2 / 1.21e-2 / 1.21e-2, N320 9.40e-3 / 4.90e-3 /
4.90e-3; linked 0.02 / 1e-3 / 1e-9: N80 1.36e-2 / 1.36e-2 / 1.29e-2, N160 8.20e-3 / 6.53e-3 /
6.53e-3, N320 5.89e-3 / – / 3.97e-3. 1e-3 = 1e-9 except at N80 (the α = 1/1024 cells).
`alphaMin 0` was broken: cells with α = 0 kept open faces (the test was α < alphaMin), unlinked
N320 inner rms 5.4e-2; fixed (cells without fluid are always closed), so 0 ≡ 1e-9 here.
Decision: default alphaMin 0 with linking, 0.02 without.

**Phase B** (7f50469, a3c5c06): as designed (B1–B6), plus four decomposition dependences found
by the tests and fixed:
1. The penalty used `eqn.diag()`, which lacks the processor-face coefficients (they live in
   `internalCoeffs`): cells next to a processor face got another conductance. Now
   `coupledDiag(eqn)` (diag + coupled internalCoeffs), also in the SA wall-law blend.
   It also fixed a link-order dependence (the diagonal grew while links were added).
2. Steady linked runs never converged: the corrector makes U uniform in a group, which no
   linked row satisfies, so the U residual stalled (2e-3) and the fixed point depended on the
   partial U solve (smoother, decomposition; 1.7e-6 between relTol 0.1 and a full solve).
   Now `linkCells(UEqn, true)` moves each group cell's residual minus its volume share of the
   group residual to the source: U converges (1e-11), faster (900–2200 iterations instead of
   not converging in 3000), independent of solver and U penalty. (The pressure penalty stays:
   mergePenalty 100 vs 1e4 changes the drag by 8e-7.)
3. Moving bodies: uncovered cells kept the solid pressure 0, a jump of the size of the
   pressure level in grad p; the level is set by the reference cell, which differs between
   decompositions (unlinked oscillating cylinder 8.5e-6 serial/parallel). Now
   `fillUncovered(p)` (mean of the fluid neighbours). Serial moving results change
   deliberately (oscillating drag at t = 10: 1.5e-5 relative).
4. Wall-law samples (B6): remote samples looked up on the other processors; the searches need
   `tetBasePtIs()`/`cellTree()` built on all ranks first (collective), and the surface level
   set communicates (one call per rank).
Results serial vs parallel: see the commit messages and `cutCellIBM.tex`, Section "Cell
linking" (cylinder 2×2/4×2/scotch 8/16 and across the cyclic, 5e-12; TC N80 second order with
10 links over processor faces, p 1e-10; oscillating cylinder 1e-12; linked scalar 1e-11;
inclined WF channel with all samples remote 2e-9; NACA WF samples identical on 1/4/8 ranks).
Not done: debug utility `-linkDebug` (the debug fields `ibmLinkRole`, `ibmLinkMaster`,
`ibmWallSampleDist` written with the geometry served instead); B5's extension to all faces
inside a group (only link faces are internal to the merged CV, as before).

**Phase C** (0fe162c): defaults mergeAlpha 0.5, alphaMin 0 (0.02 without linking); tutorials
and generators updated. Group wall distance of alpha-field bodies fixed (it averaged the
meaningless distances of non-wall cells; cylinderRe1_alphaField 0.649 → 0.705, level set 0.7195).
All 22 tutorials run; first-order drags change, e.g. cylinderRe1 0.6988 → 0.7195 (+3.0 %; the
second-order Richardson limit is about 0.7205).

**C2 reruns** (details in `doc/secondOrderStudy.md`, "Linking by default", and
`doc/iddesPlan.md` §12):
- Richardson drag: no-slip first order is much closer to the limit (0.1 % at N = 161 instead of
  1.6 %) but not monotone; no-slip second order is unchanged.
- Free-slip drag is 1–3 % lower than unlinked with both treatments and still rising at N = 321.
  Explained: the IBM free slip lacks the curvature shear (∂u_t/∂n = 0 instead of κ u_t, about 24 %
  of the drag); with an experimental zero-shear term the spread drops to 0.06–0.3 % and all series
  converge to the body-fitted value (`doc/secondOrderStudy.md`, "Free slip on curved walls"). Fix
  plan: `doc/freeSlipPlan.md` (postponed).
- TC second order: inner cut-cell rms 3.98e-3 at N = 320 (old default 9.40e-3).
- TC first order: the cut cells do not converge (rms about 0.30), and the first full layer is
  about 3× worse than unlinked. Explained (`doc/secondOrderStudy.md`, "Taylor–Couette first
  order with linking"): the first-order face-centre continuity flux, which a merged group roughly
  doubles; the masters carry the group pressure. Velocity and torque are as accurate as unlinked.
  Only the second-order face values fix it. Consequence: `secondOrder` is the default
  (`doc/secondOrderStudy.md`, "Second order is the default").
  Alpha-field bodies were not second order (no centroids, layer wall distance), not because of
  their geometry; fixed with centroids and a wall distance from the interpolant
  (`doc/secondOrderStudy.md`, "Alpha-field bodies at second order").
  First order can now be removed (decided 2026-11-04, not started; list in
  `doc/secondOrderStudy.md`, "Second order is the default").
- Inclined wall-function channel: linked errors of the same size as before; unlinked unchanged.
- The ibmMeanVelocityForce correction had to become group-uniform (8aebbce).
- NACA wall functions serial vs scotch 8/16 needed SA-neg (e99057b) to converge at all. Converged,
  the forces agree to 1e-6, limited by the iteration. One exact iteration from the same state:
  p 4e-11, U 1e-10.
- That test found the solid-row diagonal of the pressure matrix taken without the coupled-face
  coefficients (b81d779; it changes only the residual normalisation).
- NACA L5–L7 wall-resolved, linked by default: C_d +5.41 / +1.27 / +0.70 % against 897 (before
  +5.28 / +1.26 / +0.70 %). The runs clipped ν̃ throughout; with SA-neg they converge to the
  same values.
- Why linked L5/L6 were high, and the 0.13 points (2026-11-04): sliver absorption is the 0.13; the
  excess is the convective face value of slaves taken from the cell centre. `ibmConvection true`
  (default off): linked +1.13 / +0.47 / +0.61 % (`doc/secondOrderStudy.md`, "Convective face values").

## Why (results that led here)

- The cut-cell pressure error is a continuity truncation error, δp ≈ c·ν·α·div_h(u), c ≈ −1.2…−1.6
  at second order (`doc/secondOrderStudy.md`, "The cut-cell pressure is a continuity error").
  At second order its non-converging part comes from sliver absorption: a cell with α < alphaMin
  is blanked and all its faces closed (`cutCellGeometry.C`, block "Absorb slivers (C1)"),
  losing an O(shear·(θh)²) flux.
- Taylor–Couette N320, second order, inner cut cells, pressure rms / max: unlinked alphaMin 0.02
  9.4e-3 / 4.4e-2; unlinked 1e-9 4.9e-3 / 1.1e-2; linked 0.02 5.9e-3 / 1.9e-2; **linked 1e-9
  4.0e-3 / 9.3e-3** (converging, best variant). Data: `validation/data/secondOrder/taylorCouette_truncation.txt`.
- Scale: TC has ν = 1, Re 0.5; the errors scale with the wall shear stress (2.25), not with the
  pressure range (0.08). At high Re δCp ≈ (δp/τ_w)·Cf, ~1e-4. Absorption also distorts Neumann
  tests and moves the free-slip drag by 0.4–0.9 %.
- A quadratic continuity flux would not replace this: it removes only the O(h) face part (77 %
  of the variance), not the closed-face flux of absorbed slivers, not the viscous part, not first
  order, and wall functions need linking anyway (unlinked they diverge). Deferred.
- With linking, α-based absorption is obsolete: slaves use the group wall distance, the group
  is one control volume, newly uncovered cells are born as slaves (oscillating cylinder jitter
  linkSO 5.4e-3 vs baseSO 9.9e-3). It stays needed only as the fallback for cells without a
  master and for cells without an open face; cross-processor links remove the first case.

## Phase A: small fixes (independent of parallel)

A1. **Passive scalars are not linked, but see the linked wall coefficient (bug).** `calcLinks`
    sets `wallCoeff_` of slaves from the group wall distance (`cutCellGeometry.C` ~line 590, and
    `calcValueOffset` for masters); `ibmScalarTransport::updateWall` (`ibmScalarTransport.C:151`)
    uses `ibm_.wallCoeff()` while T is solved unlinked. Fix: `ibm_.linkCells(TEqn);` after
    `fvOptions_.constrain(TEqn)`, before `TEqn.solve()` (as ν̃ in `ibmSpalartAllmaras.C:155`).
    The penalty is antisymmetric: conservative, the budget stays consistent. Test:
    `tutorials/scalar/*` with mergeAlpha 0.5 (budget balance, wall flux vs unlinked), defaults
    byte-identical (mergeAlpha 0).
A2. **Guard for fluid cells without an open face.** α > 0 but θ = 0 on all faces gives an empty
    pressure row; today only caught if α < alphaMin. Blank any such cell (same code path as
    absorption: α = 0, faces closed), independent of α; count it in the log.
A3. **alphaMin study.** TC N80/160/320, second order, unlinked and linked (mergeAlpha 0.5), with
    alphaMin 0, 1e-3, 1e-9 (α is sub-sample quantised, smallest 1/1024). Scripts:
    `validation/taylorCouette_truncation.py`, the psum/classes tables of the data file. Expect
    1e-3 ≈ 1e-9; decide the default (0 with linking, small value as fallback without).
A4. **Order-independent masters.** `calcLinks` resets the link of a cell whose chain exceeds 10
    hops while looping, so later cells may stop at it: the result depends on cell numbering,
    hence on the decomposition. Define orphans from the original link graph (chain length > 10)
    before changing anything. Chains have no cycles (α strictly increases along a link), so the
    hop limit only bounds group size.
A5. **Deterministic tie-break in `outwardFace`.** Equal scores are common on symmetric geometry;
    today the first face in the cell's face list wins, and processor faces are renumbered in
    parallel. Compare scores with a relative tolerance (1e-10) and break ties by the face centre
    (lexicographic x, y, z), which is the same serial and parallel.

## Phase B: links across processor and cyclic faces

Inventory of everything that uses the groups (all assume today that the master is local and the
link face internal):

| Where | What |
|---|---|
| `cutCellGeometry.C` `outwardFace` (473), `calcLinks` (504) | link face, masters, group wall distance, slave wall coefficients |
| `cutCellGeometry.C` `calcValueOffset` (~1247) | group centroid → `valueOffset_`, master wall coefficient |
| `cutCellGeometry.C` `insideGroup` (~1306) | face inside a group (internal faces only) |
| `cutCellGeometry.C` `calcCentroidGeometry` (~1320–1410) | face value points of group faces, stats |
| `cutCellGeometry.C` gradient (~1600–1670) | least squares accumulated at the master, copied to slaves |
| `cutCellGeometryWallFunction.C` (146, 223) | sample cell's group, sample distance from `linkWallDist_` |
| `ibmWallDist/ibmLevelSetPatchDistMethod.C:67` | wall distance of group cells |
| `cutCellGeometryTemplates.C` `linkCells` (196), `groupAverage`, `groupMean`, `mergeA` (303), `mergeH` (343) | penalty, group averages, merged A and H |
| solvers `ibmPimpleFoam`/`ibmSimpleFoam` UEqn.H, pEqn.H; `ibmSpalartAllmaras.C` 92, 112, 155 | callers |

Design:

B1. **Candidates over coupled faces.** `outwardFace` also considers coupled boundary faces
    (processor, cyclic, processorCyclic) with θ > 0, using the neighbour α from
    `syncTools::swapBoundaryCellList`. Same score, same tie-break (A5). Score symmetric serial vs
    parallel because θ is already synchronised (C2 block).
B2. **Group identity.** A group id = global index of the master (`globalIndex(nCells)`), found by
    fixed-point iteration: id(c) = own global index if not a slave, else id(next(c)); next(c) is
    local or across a coupled face (swap per sweep, `reduce` a changed flag, ≤ 10 sweeps; orphans
    by A4). For comparisons between runs also store a decomposition-independent label: the
    master's cell centre (write it as a debug field `ibmLinkMasterC`).
B3. **One primitive: `groupReduce(Field<Type>&)`**, every member gets the sum over its group.
    Built once per geometry update (moving bodies: every step) as a `mapDistribute` from the
    global group ids (constructor with `globalIndex`, element list, `compactMap`): gather the
    contributions to the master's processor, sum, distribute back. If the summing reverse
    distribute is awkward in v2106, use `PstreamBuffers` with per-processor send lists built in
    `calcLinks`; the groups are few cells, speed is irrelevant. Then:
    - group wall distance, group centroid: `groupReduce(αV)`, `groupReduce(αV·y)`, `groupReduce(αV·x)`;
    - `groupAverage`, `groupMean`: sums of weighted values / weights;
    - gradient least squares: `groupReduce` of A (symmTensor), r, nUsed; every member solves the
      same system (no copy from the master needed);
    - `mergeA`, `mergeH`: per-cell contributions, then `groupReduce` (B5).
B4. **Penalty in `linkCells` over coupled faces.** Easiest correct form: build the penalty as an
    uncorrected Laplacian matrix with a conductance only on link faces, copying
    `gaussLaplacianScheme::fvmLaplacianUncorrected` (v2106
    `src/finiteVolume/finiteVolume/laplacianSchemes/gaussLaplacianScheme/gaussLaplacianScheme.C`):
    internal faces via `upper`/`negSumDiag`, coupled patches via
    `internalCoeffs = c·pvf.gradientInternalCoeffs(1)`, `boundaryCoeffs = −c·pvf.gradientBoundaryCoeffs(1)`.
    Then `eqn -= sgn·penalty`. c = mergePenalty·max(|diag_i|, |diag_j|) and the sign need the
    other side's diagonal: swap `diag` first. Handles vectors, asymmetric matrices, processors
    and cyclics in one path.
B5. **`mergeA`, `mergeH` over coupled faces.** The rule that is the same serial and parallel:
    for every face whose two sides are in the same group (`insideGroup` extended to coupled
    faces by comparing group ids), each side adds its own row's coupling over that face (A: the
    off-diagonal; H: off-diagonal × neighbour value; on a coupled face these are the patch
    `boundaryCoeffs` and `patchNeighbourField`, component average for A as in
    `fvMatrix::A()`). Signs: verify by test, not by reasoning (B7). Note: today only link faces
    are removed; faces between two slaves of one group stay as external couplings. Applying the
    rule to all faces inside a group is more consistent but changes serial linked results:
    make it a separate, reported step.
B6. **Wall-law samples are processor-local too (second decomposition dependence).** The sample
    point at `wallSampleDistance` h (default 1.5, on in the NACA and channel cases) is located
    with `mesh_.findCell` on the local mesh only (`cutCellGeometryWallFunction.C` ~140–210); a
    sample on another processor leaves the wall cell sampling itself, so a wall-function run
    changes with the decomposition even without linking. Fix together with B3: a
    `mapDistribute` from wall cells to the (remote) sample cells, found with a global search
    (e.g. `meshSearch` per processor plus a reduce over candidates, or `indexedOctree` of the
    processor bounding boxes); `group(c)` and the sample distance use the group id and the
    reduced group distance.
B7. **Tests (serial vs parallel = the acceptance criterion).**
    - Debug utility or `ibmOrderTest -linkDebug`: write ibmLinkMasterC, merged A, merged H,
      the penalty matrix applied to a smooth field, the group gradient; `reconstructPar` and compare
      with serial: link sets identical, values to round-off.
    - Decompositions: `simple` split through the body centre (2, 4, 8), `scotch` 8 and 16; the
      cut band must cross processor boundaries in many places. Log: orphans 0, link count equal.
    - Runs: cylinderRe1_parallel (drag to ~1e-9), TC N80/N160, oscillatingCylinder (moving, links
      rebuilt every step), a scalar tutorial, NACA L5 with wall functions on 8+ processors, and a
      cyclic case (channelSA) for links over cyclics.
    - Existing reference: linked cylinder on 4 procs, drag 0.719518608938 serial vs
      0.719518989557, 36 links each (the remaining 5e-7 difference should drop to solver tolerance).

## Phase C: defaults and reruns (deliberate change of all results)

C1. Defaults: mergeAlpha 0.5, alphaMin per A3 (0 with linking). Tutorials' `ibmProperties` set
    `alphaMin 0.02` explicitly: update them.
C2. Rerun and compare: cylinder Richardson (`validation/richardson_drag.sh`), TC series,
    NACA 3° L5–L7 (linked was +5.5 % vs +2.9 % at L5, 0.09 % spread at L7; linked wall-resolved
    NACA needs U/ν̃ relaxation 0.7), moving tutorials, scalar tutorials, channelSA.
C3. Report the first-order changes explicitly; update `doc/secondOrderStudy.md`,
    `doc/iddesPlan.md`, `doc/cutCellIBM.tex` (methods: linking, absorption, parallel).

Later, not in this plan: quadratic continuity and viscous fluxes, traction-form free slip,
full IDDES.
