# cutCellIBM: cut-cell immersed boundary method for OpenFOAM-v2106

OpenFOAM implementation of the cut-cell immersed boundary method prototyped in
`.asset/fvm/SIMPLEC_faceArea_noSlip.ipynb` (theory: `.asset/fvm/doc/freeslip_ibm_part2.pdf`).
Every face carries its fluid area fraction θ, every cell its fluid volume
fraction α, and the wall segment of a cut cell is *defined* by geometric closure,
`Sw = -Σ_f θ_f Sf`. From that one identity the wall pressure force, no-penetration
and the no-slip wall traction follow without ghost cells, extrapolation or
velocity-component coupling. Free-slip and no-slip walls differ by one diagonal
term.

## Layout

```
src/cutCellIBM/                 library libcutCellIBM.so
  ibmBody/                        bodies: cylinder, sphere, box, plane, any
                                  searchableSurface (STL via triSurfaceMesh), alphaField
  cutCellGeometry/                θ, α, Sw, wall distance, wall/blank coefficients,
                                  cut-cell operators (grad, flux, ddtCorr), forces
  ibmMeanVelocityForce/           constant-flow-rate fvOption, α-weighted
applications/solvers/ibmSimpleFoam    steady SIMPLE/SIMPLEC (simpleFoam + IBM)
applications/solvers/ibmPimpleFoam    transient PIMPLE (pimpleFoam + IBM, static mesh)
applications/utilities/ibmSetGeometry compute and write the geometry fields only
tutorials/                      validation cases (see below)
validation/                     notebook reference script, field comparison, run_all.sh
etc/bashrc                      installs binaries into platforms/ of this tree
```

## Build and run

```bash
module load cae/openfoam+/v2106_prf && foamInit     # or: source /opt/OpenFOAM/OpenFOAM-v2106_prf/etc/bashrc
source etc/bashrc
./Allwmake
cd tutorials/cylinderRe1 && ./Allrun                 # blockMesh, ibmSetGeometry, ibmSimpleFoam
```

`etc/bashrc` redirects `FOAM_USER_APPBIN/LIBBIN` into `platforms/` of this
directory and prepends them to `PATH`/`LD_LIBRARY_PATH`. Parallel runs work with
the usual `decomposePar; mpirun -np N ibmSimpleFoam -parallel`.

## How the method enters the solvers

The mesh is not modified. θ and α enter the standard OpenFOAM operators
explicitly, which keeps every term identical to the notebook (full cell volume V
in `fvMatrix::A()`, so `rAU = V/a_P` is exactly the notebook's `rAtU·V`):

```c++
fvm::ddt(ibm.alpha(), U)                       // α V/Δt         (ibmPimpleFoam only)
+ fvm::div(phi, U)                             // phi = θ (U_f·Sf) = ibm.flux(U)
- fvm::laplacian(ibm.thetaInterpolate(nu), U)  // θ-weighted, no wall term
+ fvm::Sp(nu*ibm.noSlipCoeff(), U)             // wall traction ν A_wall/d_wall, d_wall = αV/(2A_wall)
+ fvm::Sp(nu*ibm.blankCoeff(), U)              // u = 0 in solid cells
== fvOptions(U) - ibm.grad(p)                  // Σ_f θ_f (p_f - p_P) Sf / V, wall force by closure
```

Pressure equation: `fvm::laplacian(ibm.thetaInterpolate(rAtU), p) == fvc::div(phiHbyA)`
with `phiHbyA = ibm.flux(HbyA)` (+ the SIMPLEC term with θ), solid rows made regular
by `ibm.constrainPressure(pEqn)`. `pEqn.flux()` is then automatically θ-weighted.
Everything else (`constrainHbyA`, `setReference`, `p.relax()`, fvOptions, residual
control, PIMPLE loops) is untouched simpleFoam/pimpleFoam.

Why this and not the alternatives considered:

* **Scaling `mesh.Sf()/V()` in place** would make all operators θ/α-aware
  implicitly, but `V=αV` changes the SIMPLEC mobility (αV/a_P instead of V/a_P),
  breaks `surfaceInterpolation` weights on closed faces (0/0) and requires floors
  for solid cells. Explicit operators are transparent and term-for-term the notebook.
* **`ibm_secondOrder`** (`~/simulation/motorBike/ibm_secondOrder`) computes α with
  isoAdvector's `cutCellIso` from a signed distance at mesh points, but its
  interface-reconstruction, delta-quotient and ghost-value machinery serves a
  different (second-order ghost-cell) method and provides no face fractions.
  Only its idea of a signed distance from a `searchableSurface`
  (`findNearest` + `getVolumeType`) is reused, for STL bodies.
* OpenFOAM's multiphase reconstruction schemes (PLIC/isoAlpha) need a flux and
  velocity and give per-cell interface planes, not consistent per-face fractions;
  the notebook's edge bisection + sub-sampling is simpler and exact for analytical
  shapes, so it was ported directly (generalised to polyhedral faces/cells).

### Geometry (`cutCellGeometry`)

1. Level set of the union of all level-set bodies (`min`), evaluated at mesh
   points, face centres and cell centres.
2. Every cut edge is bisected (`nBisect`, default 40) on the *exact* level set,
   in a canonical direction so both sides of processor/cyclic faces agree.
3. θ_f = area of the fluid polygon of the face / |Sf|.
4. α_P by sub-sampling `nSubSamples^d` points per candidate cell (`d` = number
   of solution directions, so 2-D cases sample a 2-D lattice), with a
   point-in-cell test for non-hex cells.
5. Slivers with α < `alphaMin` are absorbed and *all their faces closed* (C1 of
   the paper). θ is then synchronised across coupled faces with `min` (C2).
6. `Sw = -Σ θ Sf`, `A_wall = |Sw|`, `d_wall = αV/(2A_wall)`, no-slip coefficient
   `A_wall/(d_wall V)` on no-slip cells, blanking coefficient
   `Σ|Sf| δ_f/V` in solid cells. Fluid cells behind a fully closed face are wall
   cells too (staircase limit), as in the notebook.

Body types (`constant/ibmProperties`, several bodies allowed, each with its own
`wallType noSlip|freeSlip` and optional `invert true` for fluid inside):

```
bodies
{
    cyl   { type cylinder; centre (5 5 0); axis (0 0 1); radius 1; wallType noSlip; }
    ball  { type sphere;   centre (2 2 2); radius 0.5; }
    block { type box;      min (1 1 -1); max (2 2 1); }
    floor { type plane;    point (0 0.5 0); normal (0 1 0); }   // normal points into the fluid
    hull  { type triSurfaceMesh; file "hull.stl"; }             // any searchableSurface with inside/outside
    given { type alphaField; field alpha.body; solidFraction false; }
}
```

`alphaField` takes the fluid volume fraction from a field in the start time
directory; face fractions are reconstructed from the 0.5 iso-surface of the
point-interpolated field. Bodies of different kinds are unioned by taking the
minimum of α and θ (exact where they do not share cells).

Written fields (start time): `ibmAlpha`, `ibmTheta`, `ibmSw`, `ibmWallDistance`,
`ibmNoSlipCoeff`, `ibmBodyIndex`. Forces per body (pressure `Σ p_P Sw`, viscous
`Σ ν A_wall/d_wall u_P`, per unit density) are printed every step and written to
`postProcessing/ibmForces/<time>/forces.dat`.

### Constant flow rate: `ibmMeanVelocityForce`

Derived from `meanVelocityForce`: bulk velocity and force are weighted with αV,
so solid and cut cells receive no spurious momentum, and the controller gain is
the SIMPLEC mobility 1/(A−H1) instead of 1/A. For a uniform force the neighbour
couplings cancel, so 1/(A−H1) is the exact bulk response; 1/A underestimates it
by the diffusion number ν Δt/h² when viscosity is implicit, which makes the stock
controller overshoot by two orders of magnitude at the notebook's Δt. Increments
from the predictor and corrector calls are accumulated rather than overwritten.

### Case setup notes

* `fvSchemes`: `div(phi,U) Gauss linear`, `laplacianSchemes default Gauss linear corrected`,
  `ddtSchemes default Euler` (ibmPimpleFoam) or `steadyState`; interpolation `linear`.
* `fvSolution`: `SIMPLE`/`PIMPLE` `consistent yes` recommended. In ibmPimpleFoam
  keep `ddtCorr yes` (the OpenFOAM default) for production runs: it keeps the
  Rhie–Chow coupling independent of Δt, which matters at small time steps.
  Switch it off only to compare with the notebook, which has no such term
  (`tutorials/cylinderRe1_pimple` does this). The pressure
  reference cell must be a fluid cell (checked). `U` needs no relaxation in
  PIMPLE; `relaxationFactors { equations { U 0.9; } fields { p 1; } }` for SIMPLE.
* `transportProperties` Newtonian ν; `turbulenceProperties simulationType laminar`.
  Turbulence models are constructed for `nuEff()` only; their own transport
  equations are not IBM-aware, so only laminar is supported.
* The domain boundary should not cut through a body unless the body is defined
  consistently on both sides of a cyclic boundary (a warning reports mismatches).

## Validation against the notebook (`validation/`)

`reference_notebook.py` is the notebook's numerics as a script (cells 1–6
verbatim); `compare.py` maps an OpenFOAM case cell-by-cell and face-by-face onto
the notebook arrays; `run_all.sh` reproduces everything below. Case: cylinder
R = 1 in a 10×10 doubly periodic box, 81×81 cells, ν = 1, bulk velocity 0.5
(Re_D = 1), `alphaMin` 0.02. Drag per unit depth and density; the notebook gives
5.660008 (no slip) and 2.935529 (free slip).

| quantity | ibmSimpleFoam | ibmPimpleFoam, Δt = notebook Δt | STL body (SIMPLE) | alphaField body (SIMPLE) | free slip (PIMPLE) |
|---|---|---|---|---|---|
| α, cell by cell | exact (0) | exact | exact | exact (read back) | exact |
| θ, face by face, max diff | 7.6e-11 | 7.6e-11 | 6.4e-6 (3600-gon STL) | 0.11 (iso-surface reconstruction) | 7.6e-11 |
| wall distance, max diff | 1e-11 | 1e-11 | 2.8e-6 | 2.2e-2 | 1e-11 |
| u, max diff / max u | 2.2e-4 | 9.4e-6 | 2.2e-4 | 3.3e-3 | 1.7e-5 |
| p (fluid, mean removed), max diff / max p | 3.0e-2 | 9.6e-4 | 3.0e-2 | 9.8e-2 | 1.9e-4 |
| drag | 5.660383 (+6.6e-5) | 5.659997 (−1.9e-6) | 5.660383 | 5.651222 (−1.6e-3) | 2.935496 (−1.1e-5) |

* ibmPimpleFoam with the notebook's final pseudo-time step (Δt = 65h²/ν = 0.9907),
  Euler, one corrector, SIMPLEC, no ddtCorr is the same discrete fixed point as the
  notebook; the residual differences are linear-solver tolerances.
* ibmSimpleFoam converges to a slightly different fixed point because the
  Rhie–Chow term scales with the relaxed mobility instead of Δt; the drag differs
  by 7e-5 and the velocity by 2e-4. The pressure difference is concentrated in
  cut cells (rms 8e-4) where the Rhie–Chow dissipation acts.
* 4-rank parallel run of cylinderRe1: drag identical to the serial run to all
  12 printed digits (`tutorials/cylinderRe1_parallel`).
* `twoCylinders` (no-slip R = 1 and free-slip R = 0.7, mixed wall types): the
  forces on the two bodies sum to the driving force times the fluid volume,
  0.625915 + 0.202512 = 0.828427 = 0.0703985 · 11.76767, i.e. the discrete
  momentum balance closes; the free-slip body carries no viscous force.

## Limitations / next steps

* Laminar only (see above); a turbulence model would need θ-aware transport
  equations and wall treatment.
* Static bodies. Moving bodies need the wall velocity in the traction term,
  a wall flux in continuity and recomputation of the geometry per step.
* `ddtCorr` is implemented in Euler form only (`cutCellGeometry::ddtCorr`);
  with `backward` the Euler form is used and a warning is printed once.
* Bodies smaller than a cell, or entering a cell without crossing a vertex or
  face centre, are not resolved (as in the notebook).
* Small cells are handled by absorption below `alphaMin` (cell merging would
  avoid the O(alphaMin h) boundary shift).

## Possible improvements

* **Accurate alpha-field route.** The face fractions of an `alphaField` body
  are taken from the 0.5 iso-surface of the point-interpolated α, which smears
  the wall over one cell (θ errors up to 0.1 on single faces, rms 0.05 on cut
  faces, drag −0.16% on the cylinder). To make the α route as accurate as the
  level-set routes, make θ consistent with α cell by cell: search, per cell, the
  iso-value of the point field that reproduces the cell's α (what geometricVoF's
  `isoAlpha` does), cut the faces at that iso-value, and average the two values a
  shared face receives. That brings the error down to the reconstruction's O(h²)
  level. Estimated effort: a day, including checking on the cylinder.
* **ddtCorr with backward time integration.** Add the small extension mirroring
  `backwardDdtScheme::fvcDdtPhiCorr` (two old flux levels weighted with the
  backward coefficients) to `cutCellGeometry::ddtCorr`, so that the Rhie–Chow
  correction is Δt-independent also with `backward`.
