# Adaptive mesh refinement for cut-cell IBM cases

Branch `feature_amr`. The goal is the simplest refinement that produces valid IBM meshes
automatically: start from a uniform hexahedral `blockMesh`, and the solver refines it
around the bodies at its first time step (or iteration). Optionally it then follows the
solution. The same mesh class works in 2-D and 3-D, in `ibmSimpleFoam` and
`ibmPimpleFoam`, in serial and in parallel, for static and moving bodies.

## 1. What it guarantees

- **Cut cells are at the finest level** (`maxLevel`), and so are `nCutCellLayers` finest
  cell sizes of cells on either side of every wall. Every cut cell is therefore
  surrounded by finest cells: the default of 2 layers covers the diagonal neighbours.
  After every change the solver reports either `all cut cells and their face neighbours
  at level N` or a warning that counts the offending cells (a moving body that has left
  its band: increase `nCutCellLayers` or decrease `refineInterval`).
- **2:1 balance across faces**, kept by the refinement engine; the transition from the
  finest band to the background is graded automatically.
- **Nothing changes without `constant/dynamicMeshDict`**: the solvers create a
  `staticFvMesh`, and the results are byte-identical to the code before AMR (9 tutorials,
  including moving bodies and scalars; see section 7).

## 2. Usage

```
blockMesh                     # uniform hexahedral background mesh
ibmPimpleFoam                 # or ibmSimpleFoam; refines at the first step
```

`constant/dynamicMeshDict` (re-read at every adaptation, so it can be changed during a
run):

```
dynamicFvMesh   ibmRefineFvMesh;

ibmRefineFvMeshCoeffs
{
    maxLevel        4;        // finest level: the cut cells (h = h0/2^maxLevel)
    nCutCellLayers  2;        // finest cell sizes on either side of the walls
    refineInterval  10;       // time steps (SIMPLE iterations) between adaptations;
                              // 0: only the first update (a static refined mesh)
    maxCells        1000000;  // no solution-based refinement beyond (optional)
    unrefine        true;     // coarsen where the targets have dropped (optional)
    nBufferLayers   1;        // cells around the solution-based regions (optional)

    regions                   // optional
    {
        near
        {
            type        distance;   // cells closer to a body than "distance"
            distance    0.5;
            level       3;
        }
        wake
        {
            type        field;      // lowerBound <= value <= upperBound
            field       vorticity;  // any registered field; vectors and tensors
            lowerBound  0.5;        // by their magnitude
            upperBound  1e15;
            level       3;
        }
    }

    correctFluxes   ((phi U));  // optional (this default): flux and velocity

    enableBalancing  true;    // parallel load balancing (optional, default false)
    allowedImbalance 0.2;     // balance when max/mean cells per processor - 1
                              // exceeds it (optional, this default)
    balanceMethod    scotch;  // optional (this default); coefficients of other
                              // methods from system/decomposeParDict
}
```

- `distance` regions use the signed distance to the bodies (the level sets), so they
  work for all level-set bodies (analytic shapes and STL surfaces).
- `field` regions use whatever field is in the database: a solution field (`U`, `p`,
  `nut`), or the result of a function object run every time step, e.g.
  `vorticity`, `Q`, `grad(U)`, `mag(...)`:
  ```
  functions
  {
      vorticity { type vorticity; libs (fieldFunctionObjects); executeControl timeStep; writeControl writeTime; }
  }
  ```
  A field that does not exist (yet) is skipped, and no cell is coarsened in that
  adaptation (a restart: the function object has not run yet).
- The level of a region is capped at `maxLevel`. The target of a cell is the highest
  level any criterion gives it.
- The refined mesh, `cellLevel`, `pointLevel` and the refinement history are written with
  every time directory (`<time>/polyMesh`), and a restart continues from them. The field
  `refinementLevel` is written for post-processing.
- In parallel, `enableBalancing true` redistributes the mesh and the fields across the
  processors when an adaptation has made the cell counts uneven (section 3, load
  balancing). Without it, or with it off, nothing changes.

## 3. How it works

`ibmRefineFvMesh` (`src/cutCellIBMAMR/ibmRefineFvMesh`) is a `dynamicFvMesh`. Each
adaptation has the following steps.

1. **Targets.** For every cell, `r` is the distance from its centre to its farthest
   point in the solved directions, `d` the distance of its centre to the bodies and
   `h_f` the finest cell size, `2 r/√D/2^(maxLevel − level)`. A cell with
   `d < r + nCutCellLayers h_f` gets the target `maxLevel`. This holds for every cell
   the wall passes through (the wall is within `r` of the centre), at any level, so the
   refinement converges onto the wall from a coarse background. Regions raise the target.
2. **Refinement passes.** Cells below their target are split, one level per pass, plus
   the cells the 2:1 balance needs. Passes repeat until no cell is below its target (at
   most `maxLevel + 1` passes). The first adaptation of a run therefore builds the whole
   mesh. The fields are mapped (children inherit the parent value).
3. **Coarsening, one level per adaptation.** A family is merged when none of its cells
   is needed at its level. The "keep" band is wider than the refinement band (`3r`
   instead of `r`, one more buffer layer for field regions). A merged parent would
   otherwise fall inside its own refinement band and be split again at the next
   adaptation.
4. **Rebuild.** The cut-cell geometry is recomputed on the new mesh
   (`cutCellGeometry::updateMesh`). The flux is kept where the faces did not change; on
   new, split and merged faces it is the cut-cell flux of the mapped velocity, at every
   stored time level (`correctFluxes`). The wall distance of the turbulence models is
   updated last, because `ibmLevelSet` waits for the geometry. The solver blanks the
   cells that have become solid and finds the pressure reference cell again.

The refinement engine itself, `hexRef3D` (3-D, 1 → 8) and `hexRef4` (2-D one cell thick,
1 → 4), with the 2:1 balance and the refinement history for unrefinement and restarts, is
taken unchanged from multiDimAMR (`src/cutCellIBMAMR/hexRef`, see its README).
The criteria, the balancing class (`dynamicMultiDimRefineBalancedFvMesh`, its modified
`fvMeshDistribute` and its decomposition constraint) and the axisymmetric variant of
multiDimAMR are not used; the load balancing above uses OpenFOAM's `fvMeshDistribute`
and decomposition methods directly.

### Load balancing (parallel)

Refinement around a body concentrates cells on the processors that hold it. With
`enableBalancing true`, every adaptation that changed the mesh is followed by a balance
check, before the cut-cell geometry is rebuilt:

1. **Imbalance.** `max/mean − 1` of the cells per processor is logged
   (`load balance: imbalance 3.13 (cells per processor 160 to 1762)`). Nothing happens
   at or below `allowedImbalance`.
2. **Decomposition by refinement trees.** All cells that descend from one cell of the
   background mesh form a cluster (the top ancestor in the refinement history), and the
   decomposition method (`balanceMethod`, default `scotch`) partitions the graph of the
   clusters, weighted by their numbers of cells. A family therefore never straddles a
   processor boundary, so the refinement history is distributed whole and every cell can
   still be coarsened later. This is what OpenFOAM's own history constraint for `hexRef8`
   and multiDimAMR's balancing do as well. The price is granularity: a cluster holds up to
   `(2^D)^maxLevel` cells (256 in 2-D at `maxLevel 4`), so the background mesh needs
   many cells per processor.
3. **Only if it helps.** The new imbalance is known before anything moves. The mesh is
   redistributed only if the new decomposition is within `allowedImbalance` or at least
   halves the imbalance; otherwise the log says `kept` (no repeated redistribution for a
   few percent).
4. **Redistribution.** OpenFOAM's `fvMeshDistribute` moves the mesh and all registered
   volume and surface fields, old time levels included; the refinement levels and the
   history follow as a whole (`hexRef::distribute`, from multiDimAMR). The values on the
   new processor patches are evaluated (`fvMeshDistribute` leaves them zero). Faces
   between processors become internal and the reverse, but no face is created or changed:
   the flux of the faces the adaptation left alone is carried over unchanged.
5. **Rebuild**, once, after the redistribution: the redistribution comes right after the
   topology change of the adaptation, and the flags of the faces the adaptation changed
   travel with the mesh. Then, as without balancing: cut-cell geometry, fluxes on the
   changed faces, wall distance; the solver blanks and finds the reference cell again.

`fvMeshDistribute` needs a patch that is neither empty nor coupled. In a domain periodic
in every solved direction, the first balancing adds one without faces (`loadBalancing`,
`zeroGradient` on the volume fields), and it is written with the mesh.

The redistribution is exact: per field, the sum over the cells of |value| (volume fields,
old time levels included) and of |net flux| (`phi`, `phi_0`) is the same before and after
to the last digit (`DebugSwitches { ibmRefineFvMesh 1; }` prints these sums and the
largest relative difference of the fluxes on the two sides of the coupled faces). The
one exception is a flux mapping error of `fvMeshDistribute` on new `processorCyclic`
faces, which is repaired (section 5). With linear solver tolerances of 1e-14, the first
time step of the Re 100 cylinder gives the same drag to 12 digits and lift to 4e-12 in
serial, on 6 processors without balancing and on 6 processors right after a
redistribution.

### Details that matter

- **Flux after a change.** Recomputing the flux on all faces after a change throws away the
  pressure-consistent (Rhie–Chow) flux of the unchanged faces. Every adaptation then kicks
  the drag by about 2 % (Re 100 cylinder; the median per-step change of C_d at
  adaptation steps was 1.1e-2 against 4e-4 otherwise). Rebuilding the flux only on the
  changed faces brings the kick down to 8e-5, close to an ordinary step (3e-5).
- **Parallel.** OpenFOAM computes the solved directions (`solutionD`) on demand with a
  global reduction, and a change of topology resets them. The cut-cell sampling asks for
  them only on processors that hold cut cells, so the collectives got out of step after
  an unrefinement (MPI_ERR_TRUNCATE, FPE). `ibmRefineFvMesh::updateMesh` now evaluates
  them on all processors right after every change.
- **Detecting a change.** A refinement and a coarsening in the same adaptation can leave
  the number of cells unchanged (1137 → 1146 → 1137 in the moving scalar test). A test on
  sizes then takes old data for new data: the swept-volume budget of a scalar jumped by
  5 % of the content. Changes are therefore flagged explicitly (`outOfDate()`), and
  per-cell data that must survive a change is kept in registered fields, which are mapped.
- **Moving bodies.** The adaptation runs before the bodies are moved to the new time, on
  the geometry of the old time, and `cutCellGeometry::update()` then moves them as
  before. The refined band must contain the body's motion over `refineInterval` steps
  (`nCutCellLayers` finest cells). The check after every change reports otherwise.

## 4. Code changes outside the new library

| file | change |
|---|---|
| `cutCellGeometry.H/C`, `cutCellGeometryMoving.C` | `updateMesh()`: resize the plain lists, reset the face geometry, `calcGeometry()`, old-time alpha of static bodies; `upToDate()`/`outOfDate()`, a flag set by the adaptive mesh before every change |
| `ibmLevelSetPatchDistMethod.H/C` | `movePoints()` returns `upToDate()`: the wall distance waits for the geometry |
| `ibmScalarTransport.H/C` | the fluid fractions of the previous step (swept-volume budget of moving bodies) are a registered field, mapped with the others |
| `ibmSimpleFoam.C`, `ibmPimpleFoam.C`, `createFields.H`, `Make/options` | `createDynamicFvMesh`, `mesh.update()` at the start of every step/iteration, blanking and reference cell after a change, link `libcutCellIBMAMR` |
| `Allwmake`, `Allwclean` | the new library |
| `cutCellGeometry.C` | the consistency check of the fluid fractions across coupled faces compares the geometric fractions (before sliver absorption). Before, it warned whenever an absorbed sliver cell lay next to a processor face, which is legitimate (already seen in parallel runs without balancing, frequent with it); the fractions themselves are unchanged |

## 5. Limitations

- Hexahedral background meshes only (`blockMesh`); a mesh refined by other tools without
  `cellLevel`/`pointLevel` is rejected at start. 2-D meshes must be one cell thick with
  `empty` front and back. Axisymmetric (wedge) meshes are not supported.
- Level-set bodies only: an `alphaField` body cannot be rebuilt on a changed mesh (fatal
  error at the first change).
- Mapping is first order: refined cells inherit the parent value, merged cells the
  volume average. This is harmless for the geometric band, which is built at the first
  step before anything has developed.
- `reconstructParMesh` (v2106) crashes on domains periodic in every solved direction, and
  `redistributePar -reconstruct` stops for lack of a non-coupled patch. Both happen on the
  unrefined mesh as well and are not caused by the AMR. Post-process such cases
  decomposed.
- Offline redistribution with `redistributePar` (tested with v2506) does not carry the
  refinement data of this AMR:
  - `-reconstruct` writes `cellLevel`, `pointLevel` and `level0Edge` but not the history
    (`refinementHistoryMultiDim`, a class OpenFOAM does not know). Without the history,
    no cell can be coarsened below its level at the time of reconstruction.
  - Parallel to parallel (4 → 3 ranks) deletes the refinement files and writes none. A
    restart then stops at the hexahedron check of `ibmRefineFvMesh`.
  - `fvMeshDistribute` (v2106 and v2506) needs a patch that is neither empty nor coupled.
    A domain periodic in every solved direction needs a dummy patch without faces
    (`dummy { type patch; faces (); }`; the load balancing in the solver adds one
    itself). Even with it, `-reconstruct` crashes on the periodic test case.

  The load balancing in the solver (section 3) replaces offline redistribution; a restart
  must use the number of processors the run was written with.
- Load balancing balances cells, not work, and pays only where the work per processor
  dominates. On the small 2-D Re 100 cylinder (≈ 9,000 cells on 6 or 8 processors) it
  makes the run 5–21 % slower although the imbalance drops from 0.45–0.73 to about 0.1:
  the time per step does not follow the largest subdomain there, and the scotch
  subdomains have more neighbours and processor faces than slabs (up to 20 % more
  pressure iterations). The redistributions themselves are cheap (one costs about an
  adaptation step). On the 3-D torus (205,000 cells on 8 processors) the run is 1.35
  times faster (section 7). Use it for large cases, not for small 2-D ones.
- The decomposition keeps whole refinement trees together, so its granularity is a
  cluster of up to `(2^D)^maxLevel` cells; the background mesh needs many cells per
  processor. Scotch partitions the cluster graph on the master processor (OpenFOAM's
  parallel scotch interface gathers it), which is cheap because there are only as many
  clusters as background cells.
- `fvMeshDistribute` (v2106 and v2506) can give the two sides of a new `processorCyclic`
  face different fluxes (domain periodic in every direction: both faces of such a patch,
  1.4e-5 against −3.6e-3, consistent before the redistribution). After every
  redistribution, the coupled faces whose two sides differ by more than 1e-8 (relative)
  get their flux rebuilt from the velocity, like the faces an adaptation changes, and the
  log reports their number. The check also catches small differences that the pressure
  flux of the previous step already has on a few processor faces (e.g. 2 faces with
  6e-7 in the 2-D cylinder, present before the redistribution), so the message is
  frequent in long runs; rebuilding a few faces is harmless.
- Observed while testing, not caused by the balancing:
  - Moving bodies in parallel: after some adaptations the two sides of a coupled face
    carry fluxes that differ by up to 1e-3 (also without balancing). The solver replaces
    the flux by `ibm.flux(U)` at the start of every step of a moving-body run, so it does
    not persist; its origin has not been traced.
  - The scalar source of `scalarTransportProperties` (`source { centre; radius; rate; }`)
    is selected once, by cell centre, when the solver starts. In an adaptive run it is
    selected on the background mesh and then inherited by the refined cells; a restart
    selects it again on the refined mesh, which covers a smaller area (the injection per
    step changed by a factor 1.75 in the oscillating scalar test).
- The flux rebuilt on the changed faces is not divergence-free. The next pressure
  solution corrects it, and the drag kick at an adaptation step is of the order of an
  ordinary time step (section 3).

## 6. OpenFOAM versions

Developed with v2106. All targets compile without errors with v2212, v2306, v2406, v2412
and v2506 (each in its own copy of the tree: the build directories do not carry the
version). Two changes were needed:
- `cutCellCorrectedSnGrad::corrected()` is `noexcept`, as the base class declares it
  since v2212.
- `ibmSpalartAllmaras` calls the model functions of the base class with whole fields
  from v2112 on (`SpalartAllmarasBase`, `#if OPENFOAM >= 2112`) and with internal fields
  before. The wall distance is used through a reference to its internal field, which
  is the type of `y_` in v2106. On v2106 an SA run (wall function, linking) is
  byte-identical to the base commit.

The solvers and AMR runs have been checked on v2506 (the 2-D AMR cylinder on 4 ranks,
the periodic scalar case: budget error 9e-12).

2026-11-10 (6692870): the wall-function code needs `#include "treeDataCell.H"` for
`cellTree()` (v2212 no longer includes it through `fvMatrices.H`); built against v2212
on istmdraco with `~/Codes/tools/buildtree_v2212.sh`, no warnings. Validation gates
G1–G9 and V1 on v2212 (`IBM_SUITE_FOAM` set to the v2212 environment, runs in
`~/runs/suite_v2212`): all pass, and against the v2106 baseline (e7d742a) every
regression metric is within its tolerance at zero distance; tutorial forces agree to
3e-17 (relative), the round-off Taylor-Couette forces to 2e-10. v2112 (HoreKa) is
not available locally: the same gates are to be run there before its first use.

Same on v2412 and v2506 (31d1827, `~/Codes/tools/buildtree_of.sh <tree> <version>`): no
errors, and no warnings after `hexRef` uses `Pstream::broadcast` from v2406 on (the
scatter calls are deprecated since 2024-01). G1–G9 and V1 pass on both, every regression
metric against the v2106 baseline at zero distance, tutorial forces to 1e-17, the
Taylor-Couette torque error 1.09e-4 at N = 320 as on v2106. The Re 100 AMR cylinder on 4
ranks (scotch, load balancing on, 100 steps, 10 adaptations) gives the same force history
on v2506 as on v2106 to every digit.

## 7. Tests

All runs are in `/home/ws/xt8786/runs/amr_study` (scripts `regression.sh`,
`amrtests.sh`, `shedding.py`).

| test | setup | result |
|---|---|---|
| regression, no `dynamicMeshDict` | 9 tutorials (steady, PIMPLE, STL, oscillating, rotating, eccentric rotation, two cylinders, two scalar cases), base commit 0d808c8 vs this branch | all written fields byte-identical |
| 2-D steady, body only | cylinder Re 1, periodic, background 21², `maxLevel 2` (= the 81² tutorial at the wall) | 441 → 861 cells in 2 passes, all cut cells at level 2; restart from the refined mesh changes nothing and continues smoothly |
| 2-D moving body | oscillating cylinder, background 21², `maxLevel 2`, `refineInterval 10` | 60 refinements and 29 unrefinements follow the body. Cut cells always at level 2. Force per unit depth within 1 % rms of the uniform 81² run; no force spikes at adaptation steps |
| 3-D steady, body + field | sphere Re 20, background 20×12×12, `maxLevel 3`, wake region \|U\| < 0.9, `refineInterval 100` | 2,880 → 10,552 cells at the first iteration, wake refined to 27,500; C_d = 2.75 (2 % blockage) |
| parallel 2-D / 3-D | Re 100 cylinder on 4 ranks (cut cells on processor boundaries), sphere on 4 ranks | refinement identical to serial cell for cell; unrefinement in parallel works after the `solutionD` fix |
| SA, linking | Re 1e4 cylinder, `ibmSpalartAllmaras`, `mergeAlpha 0.5` | wall distance rebuilt after the change, model runs on the refined mesh |
| STL vs analytic, 2-D | cylinder Re 1 periodic, the tutorial STL (28,800 facets) and the analytic cylinder, `maxLevel 2` | same refinement (441 → 861 cells), force differs by 2e-7 |
| STL vs analytic, 3-D | sphere Re 20 with \|U\| wake region, icosphere STL (20,480 facets) and analytic sphere | same refinement sequence (2,880 → 27,716 cells), force differs by 1.9e-4 (faceting) |
| non-convex STL, 3-D, 8 ranks | torus R = 1, r = 0.3 (6,144 facets), background 32×20×12, `maxLevel 3` | 7,680 → 205,000 cells, 15,248 cut cells, all cut cells and their neighbours at level 3 after every adaptation |
| scalars, static, 2-D | scalar tutorial `staticCylinder_fixedValue` on a 21² background, `maxLevel 3`, plume region T > 0.01, `refineInterval 5`, serial and 4 ranks | plume refined from 441 to 7,485 cells; budget error ≤ 3e-10 (content 26.3); serial and parallel content agree to 4e-9 |
| three scalars, STL, 2-D | `multiScalar` (fixed value, zero gradient, no wall) around the STL cylinder | budget errors ≤ 2e-10, 4e-12, 4e-12 |
| scalar, moving body, 2-D | `oscillating_fixedValue`, 21² background, `maxLevel 2`, 71 refinements, 29 unrefinements | budget error 1.3e-10 (content 2.14), swept-volume defect 3e-5 as without AMR |
| scalar, 3-D, STL, 1 and 4 ranks | STL sphere in a periodic box (12×8×8 background, `maxLevel 3`), upstream source, fixed-value wall, plume region | 9,600 → 46,800 cells; budget error ≤ 5e-11 (serial), 3.6e-10 (4 ranks); contents agree to 1.7e-7 |
| 2-D unsteady, wake | Re 100 cylinder, background 40×24 (h = D/2), `maxLevel 4` (D/32), vorticity region level 3, `refineInterval 10` | see below |

Re 100 cylinder (`tutorials/amr/cylinderRe100`, blockage D/H = 1/12). The adaptive run is
compared with a reference that uses the same rules but level 3 everywhere and
no adaptation (`regions near distance 100`, `refineInterval 0`). Both use D/32 at the wall
and D/16 in the wake, and are evaluated over t = 60–100 (Strouhal number from the peak of
the C_l spectrum):

| mesh | cells | C_d mean | C_l amplitude | St |
|---|---|---|---|---|
| fixed, level 3 everywhere | 61,968 | 1.4236 | 0.317 | 0.1765 |
| adaptive, vorticity > 0.5 | 8,700–9,000 (peak 10,700) | 1.4169 (−0.5 %) | 0.304 (−4 %) | 0.1764 |
| adaptive, flux recomputed on all faces (before the fix) | ≈ 9,000 | 1.4225 (t = 60–150) | 0.315 | 0.1764 |

The adaptive mesh has 7 times fewer cells than the reference. The lift amplitude is
still growing in both runs at t = 100 (adaptive, t = 60–150: 0.326). No cut cell left
level 4 in 1,500 adaptations.

### Load balancing

All in `/home/ws/xt8786/runs/amr_study/LB` (scripts `mk_re100.sh`, `mk_case.sh`, `runpar.sh`,
`launch_final.sh`, `timing.py`, `cellsproc.py`, `invariants.py`); the timings are from
one concurrent set on istmeridanus (`final_a`, 57 ranks on 64 cores).

| test | setup | result |
|---|---|---|
| no `dynamicMeshDict` | `regression.sh` (9 tutorials), `sa_reg/run.sh` (SA, wall function, linking) against the base commit | all written fields byte-identical |
| balancing off | `amrtests.sh` (steady cylinder, oscillating cylinder, 2-D and 3-D on 4 ranks) against 1d26377 | all written files byte-identical, processor directories included |
| exactness, first step | Re 100 cylinder, linear solvers to 1e-14: serial, 6 ranks `simple`, 6 ranks `scotch`, 6 ranks redistributed at the first step (imbalance 3.13 → 0.13) | drag equal to 12 digits, lift to 4e-12 |
| exactness, many steps | Re 100 cylinder to t = 3, `allowedImbalance 0.03`, debug sums | 9 redistributions, every field and old time level unchanged to the last digit |
| 2-D unsteady | Re 100 cylinder to t = 100, `simple` slabs (cylinder on one or two processors), 6 and 8 ranks | see below |
| 3-D steady | STL torus (t2b), 8 ranks, `simple` slabs along x, 300 iterations, wake region | see below |
| moving body, scalar, periodic | `oscillating_fixedValue` on a 21² background, 4 ranks, every direction periodic | `loadBalancing` patch added; 4 redistributions (0.49 → 0.003, 0.21 → 0.003, ...); 71 refinements and 29 coarsenings as in serial; budget error −1.1e-10 (serial −1.3e-10), content 2.13911914 (serial 2.13911915) |
| restart | oscillating cylinder, 4 ranks: 0 → 10 against 0 → 5 + restart 5 → 10 | same refinement sequence (60/29 against 36 + 24 / 17 + 12), same final mesh (927 → 906), drag equal to 6e-12; cells are coarsened after the restart (history read back) |
| v2506 | the moving scalar case and the Re 100 cylinder (6 ranks, on and off) | same numbers as v2106 to all printed digits |

Re 100 cylinder, t = 60–100 (imbalance: mean over the written times):

| run | imbalance | wall time | speed-up | C_d mean | C_l amplitude | St |
|---|---|---|---|---|---|---|
| serial | – | 2,276 s | 1 | 1.4169 | 0.323 | 0.1765 |
| 6 ranks | 0.45 | 1,614 s | 1.41 | 1.4169 | 0.333 | 0.1765 |
| 6 ranks, balanced (26 redistributions) | 0.10 | 1,700 s | 1.34 | 1.4171 | 0.326 | 0.1765 |
| 8 ranks | 0.73 | 1,483 s | 1.53 | 1.4169 | 0.319 | 0.1765 |
| 8 ranks, balanced (84 redistributions) | 0.09 | 1,795 s | 1.27 | 1.4170 | 0.335 | 0.1765 |

C_d and St do not depend on the balancing. The lift amplitude is still growing at t = 100
and scatters by ±2.5 % between all parallel runs, balanced or not. The balanced runs are
slower (section 5): the time per ordinary step rises from 124 to 158 ms on 8 ranks, with
135 instead of 112 pressure iterations per step, while all 77 redistributions after
t = 5 together cost 35 s.

STL torus, 8 ranks, 300 SIMPLE iterations, adaptations at iterations 1, 100, 200 and 300:

| run | imbalance after the adaptations | first iteration | ordinary iteration | total | drag force |
|---|---|---|---|---|---|
| slabs | up to 1.97 (first) and 2.70 (wake) | 716 s | 172 ms | 2,344 s | 2.3427761 |
| balanced (2 redistributions) | 0.008 and 0.010 | 352 s | 118 ms | 1,737 s | 2.3427763 |

Both reach the same mesh (204,968 cells) with all cut cells and their neighbours at the
finest level after every adaptation. The adaptations (level set of the STL surface on
every cell) dominate the run: the first iteration takes half the time balanced, the three
later adaptations 1,349 s instead of 1,576 s (their refinement passes run on the
distribution of the previous balancing, only the geometry rebuild on the new one), and
the whole run is 1.35 times faster.
