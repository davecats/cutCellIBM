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
  including moving bodies and scalars; see section 6).

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
- **Moving bodies.** The adaptation runs before the bodies are moved to the new time, on
  the geometry of the old time, and `cutCellGeometry::update()` then moves them as
  before. The refined band must contain the body's motion over `refineInterval` steps
  (`nCutCellLayers` finest cells). The check after every change reports otherwise.

## 4. Code changes outside the new library

| file | change |
|---|---|
| `cutCellGeometry.H/C`, `cutCellGeometryMoving.C` | `updateMesh()`: resize the plain lists, reset the face geometry, `calcGeometry()`, old-time alpha of static bodies; `upToDate()` (mesh sizes of the last geometry) |
| `ibmLevelSetPatchDistMethod.H/C` | `movePoints()` returns `upToDate()`: the wall distance waits for the geometry |
| `ibmScalarTransport.C` | swept-volume diagnostic restarts after a mesh change |
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
- No load balancing: refinement around a body concentrates cells on the processors that
  hold it. For large 3-D cases, refine once (`refineInterval 0`), stop, reconstruct
  and redecompose (`reconstructParMesh`, `decomposePar`).
- The flux rebuilt on the changed faces is not divergence-free. The next pressure
  solution corrects it, and the drag kick at an adaptation step is of the order of an
  ordinary time step (section 3).

## 6. Tests

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
| 2-D unsteady, wake | Re 100 cylinder, background 40×24 (h = D/2), `maxLevel 4` (D/32), vorticity region level 3, `refineInterval 10` | see below |
