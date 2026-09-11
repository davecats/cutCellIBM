# cutCellIBM: cut-cell immersed boundary method for OpenFOAM-v2106

OpenFOAM implementation of the cut-cell immersed boundary (IBM) method prototyped in
`.asset/fvm/SIMPLEC_faceArea_noSlip.ipynb` (static bodies, theory in
`.asset/fvm/doc/freeslip_ibm_part2.pdf`) and `.asset/fvm/SIMPLEC_faceArea_moving.ipynb`
(moving bodies, theory in `.asset/fvm/doc/moving_ibm.pdf`).

Every face carries its fluid area fraction θ, every cell its fluid volume fraction α,
and the wall segment of a cut cell is *defined* by geometric closure,
`Sw = -Σ_f θ_f Sf`. From that one identity the wall pressure force, no-penetration
and the no-slip wall traction follow without ghost cells, extrapolation or
velocity-component coupling. Free-slip and no-slip walls differ by one diagonal term.
A body moving through the fixed mesh changes exactly one equation: continuity acquires
the wall flux `S = u_b · Sw` as a source, and that source *is* the moving-wall boundary
condition.

Contents

1. [Layout, build, run](#1-layout-build-and-run)
2. [Static bodies: discretisation and code](#2-static-bodies)
3. [Moving bodies: discretisation and code](#3-moving-bodies)
4. [Case setup reference](#4-case-setup-reference)
5. [Coupling to an external body dynamics or an α transport equation](#5-driving-the-bodies-from-outside)
6. [Validation against the notebooks](#6-validation-against-the-notebooks)
7. [Passive scalar transport](#7-passive-scalar-transport)
8. [Limitations and possible improvements](#8-limitations-and-possible-improvements)

## 1. Layout, build and run

```
src/cutCellIBM/                 library libcutCellIBM.so
  ibmBody/                        bodies (cylinder, sphere, box, plane, any searchableSurface
                                  such as an STL, alphaField) and their rigid motion (bodyMotion)
  cutCellGeometry/                θ, α, Sw, wall distance, wall/blank coefficients, cut-cell
                                  operators (grad, flux, ddtCorr), forces; cutCellGeometryMoving.C
                                  holds everything that only runs for moving bodies
  ibmMeanVelocityForce/           constant-flow-rate fvOption, α-weighted, zero target allowed
applications/solvers/ibmSimpleFoam    steady SIMPLE/SIMPLEC (simpleFoam + IBM), static bodies
applications/solvers/ibmPimpleFoam    transient PIMPLE (pimpleFoam + IBM), static or moving bodies
applications/utilities/ibmSetGeometry compute and write the geometry fields only
tutorials/                      validation cases (section 6)
validation/                     notebook reference scripts, field comparison, run_all.sh
etc/bashrc                      installs binaries into platforms/ of this tree
.asset/fvm                      the notebooks and their documentation
```

```bash
module load cae/openfoam+/v2106_prf && foamInit   # or: source /opt/OpenFOAM/OpenFOAM-v2106_prf/etc/bashrc
source etc/bashrc
./Allwmake
cd tutorials/oscillatingCylinder && ./Allrun       # blockMesh, ibmSetGeometry, ibmPimpleFoam
```

`etc/bashrc` redirects `FOAM_USER_APPBIN/LIBBIN` into `platforms/` of this directory and
prepends them to `PATH`/`LD_LIBRARY_PATH`. Parallel runs use the usual
`decomposePar; mpirun -np N ibmPimpleFoam -parallel`. After changing a class layout in the
library, rebuild the solvers as well (`./Allwmake` does): an application compiled against
an older header allocates the geometry object with the wrong size.

## 2. Static bodies

### 2.1 Discretisation

The mesh is not modified. θ and α enter the standard OpenFOAM operators explicitly, which
keeps every term identical to the notebook (the full cell volume V is kept in
`fvMatrix::A()`, so `rAU = V/a_P` is exactly the notebook's `rAtU·V`):

```c++
fvm::ddt(ibm.alpha(), U)                       // α V/Δt              (ibmPimpleFoam only)
+ fvm::div(phi, U)                             // phi = θ (U_f·Sf) = ibm.flux(U)
- fvm::laplacian(ibm.thetaInterpolate(nu), U)  // θ-weighted, no wall term
+ fvm::Sp(nu*ibm.noSlipCoeff(), U)             // wall traction ν A_wall/d_wall, d_wall = αV/(2A_wall)
+ fvm::Sp(nu*ibm.blankCoeff(), U)              // u = 0 in solid cells
== fvOptions(U) - ibm.grad(p)                  // Σ_f θ_f (p_f - p_P) Sf / V: wall force by closure
```

Pressure equation: `fvm::laplacian(ibm.thetaInterpolate(rAtU), p) == fvc::div(phiHbyA)`
with `phiHbyA = ibm.flux(HbyA)` (+ the SIMPLEC term with θ); the empty rows of solid
cells are made regular by `ibm.constrainPressure(pEqn)`, and `pEqn.flux()` is then
automatically θ-weighted. Everything else (`constrainHbyA`, `setReference`, `p.relax()`,
fvOptions, residual control, PIMPLE loops) is untouched simpleFoam/pimpleFoam.

Why this and not the alternatives:

* Scaling `mesh.Sf()/V()` in place would make all operators θ/α-aware implicitly, but
  `V=αV` changes the SIMPLEC mobility (αV/a_P instead of V/a_P), breaks the
  `surfaceInterpolation` weights on closed faces (0/0) and needs floors for solid cells.
* `ibm_secondOrder` (`~/simulation/motorBike/ibm_secondOrder`) computes α with
  isoAdvector's `cutCellIso` from a signed distance at mesh points, but its
  interface-reconstruction and ghost-value machinery serves a different (ghost-cell)
  method and provides no face fractions. Only its idea of a signed distance from a
  `searchableSurface` (`findNearest` + `getVolumeType`) is reused, for STL bodies.
* OpenFOAM's multiphase reconstruction schemes need a flux and a velocity and give
  per-cell interface planes, not consistent per-face fractions; the notebook's edge
  bisection plus sub-sampling is simpler and exact for analytical shapes, so it was
  ported directly and generalised to polyhedral faces and cells.

### 2.2 Geometry (`cutCellGeometry::calcGeometry`)

1. Level set of the union of all level-set bodies (`min`) at mesh points, face centres
   and cell centres.
2. Every cut edge is bisected (`nBisect`, default 40) on the *exact* level set, in a
   canonical direction so that both sides of processor/cyclic faces obtain the same
   crossing.
3. θ_f = area of the fluid polygon of the face / |Sf|.
4. α_P by sub-sampling `nSubSamples^d` points per candidate cell (`d` = number of
   solution directions, so 2-D cases sample a 2-D lattice), with a point-in-cell test for
   non-hex cells. A cell is a candidate if its vertices, face centres and centre do not
   all lie on the same side.
5. Slivers with α < `alphaMin` are absorbed and *all their faces closed* (condition C1
   of the paper). θ is then synchronised across coupled faces with `min` (C2), which
   also propagates the closing of absorbed cells across cyclic and processor boundaries.
6. `Sw = -Σ θ Sf`, `A_wall = |Sw|`, `d_wall = αV/(2A_wall)`, no-slip coefficient
   `A_wall/(d_wall V)` on no-slip wall cells, blanking coefficient `Σ|Sf| δ_f/V` in solid
   cells. A fluid cell behind a fully closed face is a wall cell too (staircase limit).
   Every wall cell is attributed to a body: the one owning it if it is cut, else the one
   across the closed face, else the nearest by level set.

Written fields (start time; every write time for moving bodies): `ibmAlpha`,
`ibmTheta`, `ibmSw`, `ibmWallDistance`, `ibmNoSlipCoeff`, `ibmBodyIndex`, and for
moving bodies `ibmUb`. Forces per body (pressure `Σ p_P Sw`, viscous
`Σ ν A_wall/d_wall (u_P - u_b)`, per unit density) are printed every step and written to
`postProcessing/ibmForces/<time>/forces.dat`.

## 3. Moving bodies

### 3.1 Discretisation (`.asset/fvm/doc/moving_ibm.pdf`)

The mesh never moves and no face carries a mesh velocity: this is not an ALE scheme.
Only θ, α and the wall segment depend on time. Reynolds transport over the fluid part of
a cut cell, bounded by the wet parts of the *fixed* faces and the *moving* wall segment,
gives

| term | body at rest | body moving |
|---|---|---|
| continuity | Σ_f θ_f φ_f = 0 | Σ_f θ_f φ_f = -S, with S = u_b · Sw (Sw into the body) |
| convective | no wall term | no wall term (relative flux through the wall is zero) |
| pressure | p_P Sw absorbed by closure | same |
| viscous, free slip | no wall term | no wall term |
| viscous, no slip | ν(0 - u_P) A_wall/d on the diagonal | ν(u_b - u_P) A_wall/d: same diagonal, new source |
| unsteady | αV(u^{n+1} - u^n)/Δt | (α^{n+1} u^{n+1} - α^n u^n) V/Δt |

The quantity taken as primary is the wall flux S (from the same θ as the flux operator),
not the geometric volume change: the old fluid volume is *integrated* from it,
`α^n V = α^{n+1} V - Δt S`. This makes the pressure right-hand side exactly compatible
(Σ_P S_P telescopes to round-off), makes the divergence of a uniform u = u_b equal -S
identically, and gives a discrete geometric conservation law: a fluid translating rigidly
with the body is an exact solution for any Δt, body path and mesh. Differentiating the
sub-sampled α instead amplifies its quadrature noise by 1/Δt and fails outright.

Two further points are essential:

* the carried flux is re-interpolated from u^n through the *new* θ before the outer
  loop (a face that has just opened has no history, one that has just closed carries a
  history that no longer exists); with it the free stream is preserved from the first
  outer iteration;
* a wall must not sweep through more than a wall-adjacent cell holds,
  `Δt |u_b| ≲ αV/A_wall = 2 d_wall`; the solver reports the smallest α^n and the number
  of cells violating it every step. A handful of violating cells is harmless.

### 3.2 Code

`ibmPimpleFoam` calls `ibm.update()` right after `++runTime`. For static bodies it
returns immediately and none of the moving-body fields exist; nothing else in the static
path changed (the static tutorials reproduce their previous results to all digits). For
moving bodies `update()`:

1. moves every body to the current time (`ibmBody::moveTo`), which sets its
   displacement, velocity and rotation angle from its `bodyMotion`;
2. stores the previous θ (for `ddtCorr`) and rebuilds the whole geometry with
   `calcGeometry()`; level-set queries transform the points into the body frame
   (periodic minimum image about the current centre, inverse translation and rotation), so
   the body definitions are the ones at rest and bodies may cross cyclic boundaries
   (`periodicLengths`);
3. forms the body velocity at the cell centres (`ibmUb`, u_b + ω × r for the cells a
   body owns), the mass source per unit volume `ibmMassSource = u_b·Sw/V`, and the
   swept old fraction `α^n = α^{n+1} - Δt S/V`, which is stored in **`alpha.oldTime()`**:
   `fvm::ddt(alpha, U)` then produces the two-level unsteady term of the table without
   any change to the solver.

The solver adds, behind `if (ibm.moving())`:

```c++
phi = ibm.flux(U);                               // re-interpolate the carried flux (after update)
pRefCell = ibm.pressureRefCell(pRefCellUser);    // the user's cell while fluid, else the farthest fluid cell
UEqn -= ibm.movingWallSource(nuEff);             // ν (noSlipCoeff + blankCoeff) u_b: wall traction source and blanking to u_b
pEqn -= ibm.massSource();                        // continuity: div(phi) = -S
```

and the continuity report becomes `div(phi) + S`. Forces use `u_P - u_b` in the viscous
part. With `ddtCorr yes`, the old flux is rescaled to the new wet area
(`θ_new/θ_old`, faces that just opened get the interpolated old velocity) before the Euler
correction is formed.

With moving bodies the initial velocity field is kept as given, also inside the solid:
the blanking source drives the solid interior to u_b from the first step on, and the
free-stream test is exact only if the solid starts at u_b (the notebook does the same).
For static bodies the solid velocity and pressure are zeroed at start-up as before.

`ibmSimpleFoam` refuses moving bodies.

### 3.3 Constant flow rate: `ibmMeanVelocityForce`

Stand-alone rewrite of `meanVelocityForce` (the stock class divides by |Ubar| and cannot
target a zero bulk velocity, the natural choice for a body stirring a quiescent fluid):

* bulk velocity and force are weighted with αV, so solid and cut cells receive no
  spurious momentum, and the fluid volume is re-summed every call (it changes in time);
* the controller gain is the SIMPLEC mobility 1/(A-H1) instead of 1/A: for a uniform
  force the neighbour couplings cancel, so this is the exact bulk response, whereas 1/A
  underestimates it by the diffusion number ν Δt/h² with implicit viscosity and the stock
  controller overshoots by two orders of magnitude at the notebook's Δt;
* the row sum A-H1 can vanish or turn negative in cells a moving body sweeps through by
  more than they hold, so it is floored at the transient diagonal α/Δt like the pressure
  equation (at 10⁻³ A for steady-state schemes);
* increments of the predictor and of every corrector call are accumulated within a step;
* `Ubar (0 0 0); flowDir (1 0 0);` is accepted.

The notebook pins the bulk velocity exactly by superposition of two solves; the controller
leaves a residual of about 10⁻⁶ in the bulk velocity before each correction, which is
what limits the agreement of the force histories to ~10⁻⁴ relative.

## 4. Case setup reference

### `constant/ibmProperties`

```
alphaMin        0.02;      // slivers with less fluid are absorbed into the body
nSubSamples     32;        // samples per direction for the cell fluid fraction (2-D: 32x32)
nBisect         40;        // bisection steps locating the wall on cut edges
writeFields     true;      // write the geometry fields at start-up
writeForces     true;      // postProcessing/ibmForces/<time>/forces.dat
periodicLengths (10 10 0); // moving bodies only: minimum image in the periodic directions (0: none)

bodies
{
    cyl   { type cylinder; centre (5 5 0); axis (0 0 1); radius 1; wallType noSlip;
            motion { type oscillating; amplitude (1 0 0); frequency 0.05; } }
    ball  { type sphere;   centre (2 2 2); radius 0.5;
            motion { type translating; velocity (0.5 0 0); angularVelocity (0 0 3); } }
    block { type box;      min (1 1 -1); max (2 2 1); }
    floor { type plane;    point (0 0.5 0); normal (0 1 0); }   // normal points into the fluid
    hull  { type triSurfaceMesh; file "hull.stl"; }             // any closed searchableSurface
    given { type alphaField; field alpha.body; solidFraction false;
            moving true; velocity (0.3 0 0); }                  // or velocityField Ub;
}
```

Common entries: `wallType noSlip|freeSlip` (default noSlip), `invert true` (fluid inside
the surface), `motion { ... }`. Motion types: `none` (default), `translating`
(`velocity`), `oscillating` (`amplitude`, `frequency`, optional `phase`; d = A sin(2πft+φ)),
`prescribed` (`displacement` and `velocity` as `Function1<vector>`, e.g. `table`; the
velocity must be the exact derivative of the displacement), `external` (state set from
code, section 5). Any type accepts `angularVelocity` (constant, about the body centre; the
level set is evaluated in the rotated frame, so it works for non-axisymmetric bodies and
STLs). `alphaField` bodies take the fluid fraction from a registered field or the time
directory, reconstruct θ from the 0.5 iso-surface of its point interpolation, and with
`moving true` re-fetch it every step (their velocity is uniform or a field).

### `constant/fvOptions`

```
momentumSource
{
    type            ibmMeanVelocityForce;
    selectionMode   all;
    fields          (U);
    Ubar            (0 0 0);     // bulk velocity of the fluid
    flowDir         (1 0 0);     // required when Ubar is zero
    relaxation      1.0;
}
```

### `system/fvSchemes`, `system/fvSolution`

* `ddtSchemes default Euler` (ibmPimpleFoam) or `steadyState` (ibmSimpleFoam);
  `div(phi,U) Gauss linear`; `laplacianSchemes default Gauss linear corrected`;
  interpolation `linear`.
* `PIMPLE { consistent yes; nOuterCorrectors 3; nCorrectors 1; ddtCorr yes; }` for moving
  bodies: with a single outer iteration the convective flux would be assembled from a
  geometry one displacement out of date. `ddtCorr yes` (the OpenFOAM default) keeps the
  Rhie-Chow coupling independent of Δt; switch it off only to compare with the notebooks.
  No under-relaxation of `U` in PIMPLE; `relaxationFactors { equations { U 0.9; } fields
  { p 1; } }` for SIMPLE.
* The pressure reference cell (`pRefCell`/`pRefPoint`) must be a fluid cell; with moving
  bodies it is replaced by the fluid cell farthest from the bodies whenever a body covers
  it.
* `transportProperties` Newtonian ν; `turbulenceProperties simulationType laminar`.
  Turbulence models are constructed for `nuEff()` only; their own transport equations are
  not IBM-aware, so only laminar is supported.
* Δt for moving bodies is a physical time step, not a relaxation parameter: keep
  `Δt |u_b|` below the smallest wet cell thickness `2 d_wall` (the log reports violations).

## 5. Driving the bodies from outside

The production use foreseen is a body whose position comes from a point-particle
equation integrated with the flow, or an α field solved by a separate scalar equation.
Both go through the same hook, `cutCellGeometry::update()`, which reads the bodies' state
at the start of every step:

* **Rigid body driven by the solver.** Give the body `motion { type external; }` and, in
  the solver, before `ibm.update()`:
  ```c++
  ibm.bodies()[i].motion().setState(displacement, velocity);   // displacement from the position at rest
  ```
  The hydrodynamic force needed by the particle equation is available from
  `ibm.forces(p, U, nu, Fp, Fv)` (per body, per unit density) at the end of the previous
  step; a light body would need an added-mass sub-iteration inside the outer loop.
* **α from a transport equation.** Register the solved field on the mesh under the name
  given in `field` and declare the body `type alphaField; moving true;` with its velocity
  (uniform `velocity`, or `velocityField` naming a `volVectorField` on the mesh).
  `update()` then re-fetches the field every step, reconstructs θ, and forms S and α^n
  from it exactly as for level-set bodies. Note that S is still taken from the wall flux
  u_b·Sw and α^n integrated from it, never from the difference of the supplied α fields
  (section 3.1), so the supplied field only needs to be right at the new time level.

## 6. Validation against the notebooks

`validation/reference_notebook.py` and `validation/reference_moving.py` are the two
notebooks' numerics as scripts (their cells copied verbatim, plotting removed, results saved);
`compare.py` / `compare_moving.py` map an OpenFOAM case cell by cell and face by face onto
the notebook arrays; `run_all.sh` reproduces everything below. Common setting: cylinder
R = 1 in a 10×10 doubly periodic box, 81×81 cells, ν = 1, `alphaMin` 0.02.

### 6.1 Static bodies (bulk velocity 0.5, Re_D = 1; notebook drag 5.660008 no slip, 2.935529 free slip)

| quantity | ibmSimpleFoam | ibmPimpleFoam, Δt = notebook Δt | STL body (SIMPLE) | alphaField body (SIMPLE) | free slip (PIMPLE) |
|---|---|---|---|---|---|
| α, cell by cell | exact (0) | exact | exact | exact (read back) | exact |
| θ, face by face, max diff | 7.6e-11 | 7.6e-11 | 6.4e-6 (3600-gon STL) | 0.11 (iso-surface reconstruction) | 7.6e-11 |
| u, max diff / max u | 2.2e-4 | 9.4e-6 | 2.2e-4 | 3.3e-3 | 1.7e-5 |
| drag | 5.660383 (+6.6e-5) | 5.659997 (−1.9e-6) | 5.660383 | 5.651222 (−1.6e-3) | 2.935496 (−1.1e-5) |

ibmPimpleFoam with the notebook's final pseudo-time step (Δt = 65h²/ν, Euler, one corrector,
SIMPLEC, no ddtCorr) is the same discrete fixed point as the notebook; ibmSimpleFoam differs
through the relaxation dependence of the Rhie–Chow term. A 4-rank run reproduces the serial
drag to all printed digits, and in `twoCylinders` (no-slip R = 1 plus free-slip R = 0.7) the
forces on the two bodies sum to the driving force times the fluid volume exactly.

### 6.2 Moving bodies (Δt = 0.025, 3 outer iterations, no slip, bulk velocity pinned to 0)

**Free-stream preservation** (`tutorials/freeStreamTranslating`): domain filled with
u = u_b, body translating at u_b, 5 steps. max |u − u_b| over the fluid:

| u_b | nOuter = 1 | nOuter = 3 | notebook |
|---|---|---|---|
| 0.37 | 0.0 (v: 7e-16) | 0.0 (v: 5e-16) | 5.6e-16 |
| 2.00 | 0.0 (v: 1e-15) | 0.0 (v: 2e-15) | 3.6e-15 |

**Oscillating cylinder** (`tutorials/oscillatingCylinder`): x = 5 + sin(2π·0.05 t),
t = 0…10 (400 steps). Force history: max |ΔF_x| = 4.4e-4 after the impulsive first 10 steps
(7e-5 of the peak force; the first step, F = −60.6, differs by 0.19%); driving force
3.86362e-2 vs 3.86361e-2; final fields at t = 10: u to 2.0e-7, p to 1.3e-5, α identical.
Mass source sum and continuity error at round-off every step.

**Galilean invariance** (`tutorials/galileanTranslating`): body translating at 0.5 through
a fluid of zero bulk velocity, t = 0…30, crossing the periodic seam twice. Mean force over
the last 400 steps: −5.69099 ± 0.0883 (OpenFOAM) vs −5.69040 ± 0.0884 (notebook as
published) and −5.69103 ± 0.0884 (notebook with the seam fix below); the fixed-body drag is
+5.66038. Final fields at t = 30 agree to 7.5e-5 in u. Away from the seam the histories agree
to ~1e-4 relative; while the body front or back crosses the seam (body centre 9.0–9.6,
10.4–11.0, 19.0–19.6) the published notebook deviates by up to 2.7% for a few steps,
because it copies θ across the periodic seam *before* absorbing slivers: an absorbed cell at
x = L closes its face on that side only, and the same face stays open for the cell at x = 0.
OpenFOAM synchronises θ across coupled faces with `min` *after* absorption, which keeps
closure exact (the mass-source sum is 1e-19 throughout). `reference_moving.py --seamFix`
applies the same rule to the notebook; with it the histories agree to 1.3e-3 at every step after the impulsive start (2.2e-4 of the mean force), mean −5.69103 ± 0.0884.

A 4-rank run of the oscillating cylinder reproduces the serial force history to 2e-12.

## 7. Passive scalar transport

`ibmScalarTransport` (library) solves one passive scalar T with the cut-cell operators;
`ibmScalarTransportList` reads `constant/scalarTransportProperties` (any number of scalars,
nothing if the file is absent) and both solvers call `scalars.solve()` once per time step
(ibmPimpleFoam, after the PIMPLE loop) or per iteration (ibmSimpleFoam). The immersed wall
is seen, per scalar and per body, either as a fixed-value wall (the analogue of no slip)
or as a zero-flux wall (the analogue of free slip):

```c++
fvm::ddt(ibm.alpha(), T) + fvm::div(phi, T) - fvm::laplacian(DT*ibm.theta(), T)
+ fvm::Sp(DT*wallCoeff, T) + fvm::Sp(DT*ibm.blankCoeff(), T)
== DT*(wallCoeff + ibm.blankCoeff())*Twall + q*ibm.alpha() + fvOptions(T)
```

* `wallCoeff` is the geometric Dirichlet coefficient A_wall/(d_wall V) (`ibm.wallCoeff()`,
  the same quantity that becomes `noSlipCoeff` for the velocity) masked to the wall cells
  of bodies seen as fixed-value walls; for zero-gradient walls there is no term at all,
  because the wall segment simply carries no flux. The masks of velocity and scalar are
  independent, so a free-slip body may hold a fixed scalar and vice versa, and each scalar
  chooses its own wall type per body.
* The blanking term drives the solid interior to the wall value (0 for zero-flux
  bodies), so a cell released by a moving body starts from a defined value.
* Cells containing the wall are fluid cells: T is transported into their fluid part
  through the wet fraction of their faces (the same `phi` the pressure equation returns),
  and never across the wall segment. This is what "the scalar is not lost into the solid"
  means discretely.
* For a moving body no mass-source term is added: Reynolds transport over the fluid part
  of a cell gives exactly this form, and a uniform T is preserved because
  (α^{n+1} − α^n)V/Δt = S and Σθφ = −S cancel.
* `q` is a prescribed source per unit fluid volume in a circular (2-D) or spherical region;
  it multiplies α, so it never injects into the solid.
* `fvOptions(T)` is the standard OpenFOAM mechanism: any source in `constant/fvOptions`
  whose `fields` list contains the scalar enters its equation (`scalarSemiImplicitSource`
  for explicit/implicit rates, `scalarCodedSource` for arbitrary expressions, ...), with the
  usual `constrain`/`correct` calls. Its contribution is counted in the budget. Two things
  to keep in mind: stock fvOptions act per cell volume V, not per fluid volume αV, and also
  in solid cells (where the blanking term wins), so a source that overlaps the body is
  first-order inconsistent there unless it is made α-aware like `ibmMeanVelocityForce`; and
  the coupling is one-way, the flow does not see the scalars.

Setup: one `0/<name>` field per scalar, and

```
// constant/scalarTransportProperties
scalars
{
    T
    {
        DT      0.1;
        source  { centre (5 8.5 0); radius 0.5; rate 1; }     // optional
        walls   { cylinder { type fixedValue; value 0; } }    // per body, optional
    }
    c   { DT 0.01; walls { ".*" { type zeroGradient; } } }
}
```

A body not matched in `walls` takes `scalarWall fixedValue|zeroGradient; scalarWallValue 0;`
from its own entry in `constant/ibmProperties`, and zero gradient if neither is given. The
single-scalar form `field T; DT 0.1; source {...};` at the top level of the file is still
accepted. Each scalar needs `div(phi,<name>) Gauss linear` in `fvSchemes` (a regex such as
`"div\(phi,(T|c)\)"` does), a solver entry in `fvSolution`, and `1.0` in the equation
relaxation factors for ibmSimpleFoam (the steady equation is linear at fixed flux;
under-relaxing it only slows the adjustment of the global level, which is set by the weak
wall sink). The scalars are independent of each other; sources coupling them (reactions)
go through fvOptions, or through a term added in `ibmScalarTransport::solve()`.

### Budget check

Every step the solver prints, per scalar, the fluid content Σ αVT, the injected amount
Σ qαV dt (plus what fvOptions added or removed), the
amount that left through fixed-value walls Σ D_T·wallCoeff·V·(T − T_wall) dt, the
swept-volume defect Σ (α_prev − α^n_swept) V Tⁿ, and their balance. The defect exists only
for moving bodies: the old fluid volume used by the momentum and scalar equations is
integrated from the wall flux (section 3.1) and is first-order consistent with the previous
geometry, and a cell absorbed into the body hands over no content. The balance must close
to round-off in all cases. `validation/scalar_budget.py` collects the last line of every
case in `tutorials/scalar` (Δt and durations as in the flow cases; source at (5, 8.5),
radius 0.5, rate 1, D_T = 0.1; T = 0 initially):

| case | steps | content | injected | through walls | swept defect | budget error | max |T − T_wall| in the solid |
|---|---|---|---|---|---|---|---|
| staticCylinder_zeroFlux (PIMPLE, Δt = 0.99) | 202 | 19.581366 | 19.581366 | 0 | 0 | −3.9e-12 | 0 |
| staticCylinder_fixedValue | 202 | 11.904949 | 19.581366 | 7.676417 | 0 | −5.7e-11 | 0 |
| steadyCylinder_fixedValue (SIMPLE) | 1834 | 16.323010 | rate 0.097847 | rate 0.097847 | – | −1.0e-11 | 0 |
| oscillating_zeroFlux | 400 | 0.978439 | 0.978472 | 0 | 3.3e-5 | −8.0e-11 | 1e-50 |
| oscillating_fixedValue | 400 | 0.968143 | 0.978472 | 0.010344 | −1.5e-5 | −6.1e-11 | 3e-56 |
| galilean_zeroFlux (seam crossings) | 1200 | 2.935650 | 2.935415 | 0 | −2.3e-4 | −7.5e-10 | 5e-13 |
| multiScalar, T (fixed value) | 202 | 11.904949 | 19.581366 | 7.676417 | 0 | −5.7e-11 | 0 |
| multiScalar, c (zero flux) | 202 | 19.581366 | 19.581366 | 0 | 0 | −3.9e-12 | 0 |
| multiScalar, d (zero flux, fvOptions decay 0.05/s) | 202 | 1.953492 | 1.953492 | 0 | 0 | −1.9e-12 | 0 |

With zero-flux walls and a static body the content equals the injected amount to 4e-12
after 200 steps, with the solid exactly at zero: no scalar is lost into the solid. With a
moving body the content differs from the injected amount by the swept-volume defect only
(3e-5 of the content for the oscillating cylinder, 8e-5 for the translating one), and the
solid stays at the wall value; the budget errors of 1e-10 are the accumulated linear-solver
tolerance. The flow is unaffected by the scalar (the forces of `oscillating_zeroFlux` equal
those of `oscillatingCylinder` to the last digit). `validation/fig_scalar.png` shows the
fields. `multiScalar` carries three scalars in one run: T and c reproduce the budgets of
`staticCylinder_fixedValue` and `staticCylinder_zeroFlux` to the last digit, and d adds a
first-order decay through `scalarSemiImplicitSource` (`fields (d)` only), whose sink is
accounted in the injected column so that the budget still closes.

## 8. Second-order wall treatment (`secondOrder true`)

`doc/secondOrderStudy.md` identifies what limits the order of the no-slip wall and tests the
remedy; this branch implements it. Switch it on per case with

```
// constant/ibmProperties
secondOrder   true;
// system/fvSchemes
laplacianSchemes { default Gauss linear cutCellCorrected; }
snGradSchemes    { default cutCellCorrected; }
```

With `secondOrder false` (the default) everything is bit-identical to before, and the
`cutCellCorrected` scheme then reduces to the plain corrected one.

What it does:

* the fluid centroid of every cut cell (from the sub-sampling) and the wet-face centroid
  of every cut face (from the polygon clipping) are stored (`ibmCentroidOffset`);
* the wall distance `d_wall` is the level set at the fluid centroid instead of
  `αV/(2A_wall)`;
* the `cutCellCorrected` snGrad scheme makes the diffusive fluxes across the faces of cut
  cells consistent with the centroids: implicit coefficient θA/(n·d) with d the
  centroid-to-centroid vector, explicit remainder `g_f·n − (g_f·d)/(n·d)` with least-squares
  gradients on the centroids. Through OpenFOAM's Gauss Laplacian this covers the momentum and
  pressure Laplacians, `pEqn.flux()` and `fvc::snGrad`, so the Rhie–Chow coupling stays
  consistent; the stencil does not grow and the implicit part keeps its diagonal dominance;
* the cut-cell gradient uses face values at the wet-face centroids and the wall pressure
  extrapolated to the wall segment, both from the least-squares gradient; the force uses the
  same wall pressure.

Verification (`tutorials/order`, `validation/richardson_drag.sh`, N = 20…160 for the operator
tests, α_min = 0.02):

| test | current | second order |
|---|---|---|
| Laplacian, Dirichlet (`ibmOrderTest LD`), L2 orders | 0.65, 1.22, 1.04 | 1.92, 2.06, 2.09 |
| convection–diffusion, Dirichlet (`CD`) | 0.55, 1.04, 1.21 | 1.76, 1.94, 1.99 |
| Laplacian, Neumann (`LN`) | 2.27, 2.43, 2.96 | 2.08, 2.33, 2.47 |
| convection–diffusion, Neumann (`CN`) | 1.94, 1.90, 1.82 | 1.91, 1.90, 1.83 |
| Taylor–Couette (exact NS), L2(U), N = 40/80/160 | 6.5e-3, 4.5e-3 (0.5), 2.1e-3 (1.1) | 1.2e-3, 3.4e-4 (1.8), 1.3e-4 (1.4) |
| Taylor–Couette, bulk cells, U rms, N = 80→160 | order 1.0 | order 1.7 |
| cylinder drag, no slip, Richardson N = 41/81/161 | 5.5431, 5.6604, 5.7577 (p = 0.3) | 5.7834, 5.8205, 5.8322 (p = 1.7) |
| cylinder drag, free slip | 2.8377, 2.9421, 3.0004 (p = 0.8) | 2.9161, 2.9848, 3.0237 (p = 0.8) |

The operators are second order for both wall conditions. In the coupled Navier–Stokes problem
the velocity is second order in the bulk and 5–16 times more accurate everywhere, but the
pressure of the cut cells does not converge (its mobility scales with αθ, so it is set by tiny
flux residuals: the small-cell problem of the pressure) and pollutes the near-wall pressure and
the forces at first order. A second-order force needs a cut-cell pressure treatment (cell
merging, or a wall pressure formed from the full cells); see the study document.

The Taylor–Couette test also exposed a defect of the moving-body path, fixed here: the body
velocity of a rotating wall was evaluated at the cell centre instead of the wall segment,
which created an O(h) spurious mass source alternating in sign around the body. Translating
bodies are unaffected (the moving-body tutorials reproduce their results to the last digit).

## 9. Limitations and possible improvements

* Laminar only; a turbulence model would need θ-aware transport equations and wall
  treatment.
* First order in time (implicit Euler, and the swept-volume integration inherits it) and,
  for no slip, first order in space; free slip is second order.
* Slivers are absorbed below `alphaMin` and lose their momentum when absorbed; mass
  conservation is unaffected. Cell merging would remove the O(alphaMin h) boundary shift
  and relax the step restriction.
* `ddtCorr` is implemented in Euler form only; with `backward` the Euler form is used and a
  warning printed once. The extension mirrors `backwardDdtScheme::fvcDdtPhiCorr`.
* The alpha-field route reconstructs θ from the 0.5 iso-surface of the point-interpolated
  α (θ errors up to 0.1 on single faces, drag −0.16% on the cylinder). Making θ consistent
  with α cell by cell (per-cell iso-value search as in geometricVoF's `isoAlpha`, faces cut
  at that value, shared faces averaged) would bring it to O(h²).
* Prescribed rigid motion; the `external` motion type and `forces()` are the hooks for a
  rigid-body ODE, and a light body would need an added-mass sub-iteration.
* One fluid region: a single pressure reference cell; a body splitting the fluid in two
  would need one reference per region.
* Geometry cost for moving bodies: the level set is evaluated at all points, face and cell
  centres every step and `nSubSamples^d` samples per cut cell; for 3-D STL bodies reduce
  `nSubSamples` (e.g. 8) or restrict the evaluation to a band around the previous
  position (not implemented).
