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
The criteria, the load balancing and the axisymmetric variant of multiDimAMR are
not used.

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
    (`dummy { type patch; faces (); }`). Even with it, `-reconstruct` crashes on the
    periodic test case.
- No load balancing: refinement around a body concentrates cells on the processors that
  hold it. For large 3-D cases, refine once (`refineInterval 0`), stop, reconstruct
  and redecompose (`reconstructParMesh`, `decomposePar`).
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
