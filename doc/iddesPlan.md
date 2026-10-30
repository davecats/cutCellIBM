# Plan: Spalart–Allmaras IDDES with the cut-cell IBM (branch `feature_iddes_from_scalar_secondOrder`)

Goal: run `ibmPimpleFoam` (and `ibmSimpleFoam` for RANS) with Spalart–Allmaras (SA) RANS and
SA-IDDES, with immersed walls treated as consistently as body-fitted walls, including wall
functions; validate against body-fitted OpenFOAM runs of the same model on cheap cases.

Status of the code on the branch point (`scalar_secondOrder`, af89b23): the solvers construct
the turbulence model and call `turbulence->correct()`, but only `simulationType laminar` is
supported (warning at start-up): the model's own equation does not see the body, its wall
distance comes from the domain wall patches only, and the momentum equation uses the *cell*
`nuEff` in the no-slip wall traction and omits the transpose part of the viscous stress.

## 1. What OpenFOAM's SA models need from a wall

In v2106, `SpalartAllmarasDES::correct()` (inherited by DDES and IDDES; the RAS model
`SpalartAllmaras` is the same structure) assembles

    ddt(nuTilda) + div(phi, nuTilda) − laplacian(DnuTildaEff, nuTilda) − Cb2/σ |∇nuTilda|²
      = Cb1 S̃ nuTilda (1 − ft2) − Sp((Cw1 fw − Cb1/κ² ft2) nuTilda/d̃², nuTilda)

with standard operators, then `nut = nuTilda fv1`. A body-fitted wall enters through three
things only:

| ingredient | body-fitted OpenFOAM | what the IBM must supply |
|---|---|---|
| wall distance `y` (in d̃, fw, IDDES blending) | `wallDist` (meshWave from wall patches) | distance to the immersed surface |
| `nuTilda` at the wall | patch `fixedValue 0` | Dirichlet 0 on the cut-cell wall |
| wall shear in momentum | patch `nut` from a wall function (`nutUSpaldingWallFunction`, `nutLowReWallFunction`) | wall-function viscosity in the cut-cell wall traction |

The `alpha` in that equation is a phase fraction (`geometricOneField` for incompressible flow), so
the IBM cannot be injected through it: `correct()` has to be overridden. Everything else (`chi`,
`fv1`, `fv2`, `Stilda`, `fw`, the virtual `dTilda`, the IDDES blending and `correctNut`) is
protected and can be reused unchanged.

## 2. Strategy

### 2.1 Wall distance: a `patchDistMethod` (no change to the turbulence models)

`y_` of the SA models is a reference to `wallDist::New(mesh).y()`, and `wallDist` computes it with
a run-time selectable `patchDistMethod` (`fvSchemes: wallDist { method meshWave; }`). Add a method
`ibmLevelSet` to our library:

    y = min( meshWave distance to wall patches (if any),  φ_IBM(x) )

with φ_IBM the union level set of the bodies, which is a signed distance for every level-set body
type in the library (cylinder, sphere, box, plane, triangulated surface via `findNearest`). Solid
cells get a small positive floor (they are blanked anyway). Evaluate at the fluid centroid when
`secondOrder` is on. Moving bodies: `wallDist::movePoints()` is public; the IBM-aware model calls
it after `ibm.update()`. Alpha-field bodies have no distance; for them y is taken from a
meshWave pass seeded at the wall cells with the wall distance of the cut (first version: not
supported, fatal error with a clear message).

Alternative rejected: computing y inside the derived model and overriding every function that uses
`y_` (IDDES uses it in several blending functions; fragile against OpenFOAM versions).

### 2.2 The `nuTilda` equation: derived models that override `correct()` only

New classes in a separate library `libcutCellIBMTurbulence.so` (keeps the core library free of the
turbulence libraries; loaded with `libs (...)` in `controlDict`):

- `ibmSpalartAllmaras` (RAS) derived from `SpalartAllmaras`,
- `ibmSpalartAllmarasIDDES` (LES/DES) derived from `SpalartAllmarasIDDES` (DDES/DES come for free
  the same way if wanted).

Their `correct()` is the original one with the IBM operators, exactly as for the passive scalars
(`ibmScalarTransport`), which already handle Dirichlet walls, blanking and swept volumes:

    fvm::ddt(ibm.alpha(), nuTilda)                         // swept old volume for moving bodies
  + fvm::div(phi, nuTilda)                                 // phi = ibm.flux(U), already θ-weighted
  − fvm::laplacian(ibm.thetaInterpolate(DnuTildaEff), nuTilda)
  + fvm::Sp((nu/σ) ibm.wallCoeff(), nuTilda)               // Dirichlet nuTilda_w = 0
  + fvm::Sp(nu ibm.blankCoeff(), nuTilda)                  // nuTilda = 0 in solid cells
  − Cb2/σ alpha |∇_ibm nuTilda|²
  = alpha [ Cb1 S̃ nuTilda (1−ft2) ] − fvm::Sp(alpha (Cw1 fw − …) nuTilda/d̃², nuTilda)

Sources weighted with α (as the scalar sources), `nut = 0` in solid cells, and ∇ for the Cb2 term
and for S̃ (velocity gradient) taken with the cut-cell gradient using the wall value (0 for
nuTilda, u_b for U) instead of the zero-gradient wall of `fvc::grad`; this needs a small extension
`cutCellGeometry::grad(field, wallValue)`. The diffusivity at the wall is ν/σ because nuTilda_w = 0.
With `secondOrder true` the `cutCellCorrected` snGrad scheme and the centroid wall distance apply to
nuTilda automatically.

The IBM geometry is found through the registry (`mesh.lookupObject<cutCellGeometry>("ibmProperties")`);
the model refuses to run without it.

### 2.3 Wall functions: an IBM counterpart of `nutUSpaldingWallFunction`

Body-fitted SA practice is `nuTilda` fixedValue 0 plus a `nut` wall function on the patch; the wall
function only changes the shear on the wall face. The IBM analogue changes only the wall traction
of the cut cell, `ν Awall/d (u_P − u_b)`:

- per wall cell, solve Spalding's law for u_τ from the tangential slip |u_t,P − u_b,t| at distance d
  (the velocity sample point), as `nutUSpaldingWallFunction` does per face, and set the wall
  viscosity ν_w = u_τ² d/|u_t| (≥ ν); Spalding reduces to u⁺ = y⁺ for small y⁺, so the same law
  covers wall-resolved (y⁺ ≈ 1) and wall-modelled (y⁺ ≈ 30–300) meshes;
- momentum: `fvm::Sp(ν_w noSlipCoeff, U)` and the u_b source with ν_w (instead of the cell nuEff);
  the viscous force uses the same ν_w;
- per body in `ibmProperties`: `wallFunction spalding;` (default for turbulent runs) or `none`
  (ν_w = ν, wall-resolved low-Re).

Design choices to settle with the tests:
1. *Sample point.* First version: the cell value at distance d (consistent with the body-fitted wall
   function, which samples the wall cell). Known IBM weakness: tiny cut cells give small d and
   velocities pinned by the wall term, i.e. noisy τ_w. If the channel and airfoil tests show this,
   second version: an image point at a fixed distance d_s ≈ 1.5–2 h along n_w, velocity interpolated
   there, τ_w from (u_s, d_s), imposed on all wall cells (Capizzano 2011; Tamaki & Imamura 2017).
2. *Normal component.* The traction acts on the whole slip vector; with ν_w ≫ ν it also stiffens the
   normal component. Option: apply ν_w to the tangential part only (explicit correction). Decide
   from the channel test with a tilted wall.
3. *nuTilda in the first cell* for wall-modelled meshes: keep Dirichlet 0 at the wall (what the
   body-fitted reference does); the equilibrium value κu_τd as an alternative only if needed.

### 2.4 Solver changes (`ibmSimpleFoam`, `ibmPimpleFoam`)

- wall traction coefficient and u_b source from the IBM wall function when the model is IBM-aware;
- the missing transpose viscous term, `div(nuEff dev2(T(∇U)))`, assembled with θ-weighted faces and
  the cut-cell gradient (zero for constant ν, so laminar results stay bit-identical);
- the "laminar only" warning becomes an error for a non-IBM turbulence model, silent for the IBM ones;
- forces: viscous part with ν_w.

Laminar regression (all tutorials, order tests) must stay bit-identical.

## 3. Test cases (cheap first)

### T1. Turbulent plane channel, SA RANS (seconds per run) — isolates the wall treatment

Periodic channel at Re_τ = 550 and 2000, `ibmSimpleFoam` with `ibmMeanVelocityForce`, one cell in
x: a 1-D problem. Body-fitted: walls as patches, `nutUSpaldingWallFunction`. IBM: the same column
extended by a few cells, walls as two `plane` bodies placed *off* the mesh lines (offsets 0.1, 0.5,
0.9 of a cell, to expose the cut-cell sensitivity). Two resolutions: y⁺₁ ≈ 1 (wall-resolved,
`wallFunction none` vs `spalding`) and y⁺₁ ≈ 50 (wall-modelled). Compare u⁺(y⁺), ν_t, u_τ, bulk
velocity. Variant: channel tilted by 30° in a doubly periodic box (walls not aligned with the mesh).

### T2. NACA 0012, 2-D steady SA RANS — body-fitted vs IBM

Fully turbulent, inflow at angles of attack α = 0° and 10° (rotate the inflow, not the body, so
both meshes keep the same airfoil). Closed trailing edge as in the NASA Turbulence Modeling
Resource (TMR) "2DN00" case.

- Body-fitted reference: C-mesh (blockMesh or the TMR grid family converted with
  `plot3dToFoam` + extrusion), `simpleFoam` + `SpalartAllmaras`, both wall-resolved (y⁺ ≈ 1, the
  "truth") and with `nutUSpaldingWallFunction` at y⁺ ≈ 50 (the like-for-like comparison).
- IBM: Cartesian mesh graded towards the airfoil (blockMesh grading, or `refineMesh` boxes once the
  IBM is checked on split-hex cells), the airfoil as an extruded STL generated by a small Python
  script (`type triSurfaceMesh`), `ibmSimpleFoam` + `ibmSpalartAllmaras` + Spalding.
- Reynolds number: Re = 1·10⁶ as the primary case (y⁺ ≈ 50 needs h ≈ 2·10⁻³ c at the wall: about
  10⁴ cells in a 0.02 c band, 3–5·10⁴ cells in total, minutes per run); Re = 6·10⁶ second, for
  comparison with the TMR SA results and experiments (Ladson; Gregory & O'Reilly) at a coarser y⁺.
- Metrics: C_l, C_d (pressure/viscous split), C_p(x), C_f(x), velocity profiles at x/c = 0.5 and 0.9;
  three mesh levels for each approach.

OpenFOAM's `simpleFoam/airFoil2D` tutorial (SA + Spalding, Re ≈ 2.6·10⁶) is a ready body-fitted
fallback to get started.

### T3. IDDES

1. *Shielding check (cheap)*: the T2 airfoil at α = 10° with IDDES in `ibmPimpleFoam` on a thin
   spanwise slab: an attached boundary layer must stay in RANS mode (same C_l, C_d as SA RANS, no
   grid-induced separation). Exercises the IBM wall distance in the IDDES blending functions.
2. *Wall-modelled LES branch*: periodic channel at Re_τ ≈ 2000 with IDDES (the WMLES test of
   Shur et al. 2008), about 2·10⁵ cells, body-fitted vs IBM walls (plane bodies, off-grid).
   Compare mean u⁺, Reynolds stresses, u_τ.
3. *Optional demonstration*: NACA 0012 in deep stall (α ≈ 45°–60°, span 1 c periodic), the classic
   DES case; only after 1–2, as it costs 10⁶ cells and long averaging.

## 4. Work packages and order

| # | content | depends on |
|---|---|---|
| 1 | `ibmLevelSet` patchDistMethod; check y against meshWave on the body-fitted channel | — |
| 2 | `cutCellGeometry::grad(f, wallValue)`; library `cutCellIBMTurbulence` with `ibmSpalartAllmaras` | 1 |
| 3 | IBM wall function (Spalding, cell sample), solver changes 2.4, laminar regression | 2 |
| 4 | T1 channel (RANS): body-fitted vs IBM, y⁺ ≈ 1 and 50, wall offsets, tilted wall | 3 |
| 5 | image-point sampling and tangential-only ν_w, if T1 asks for it | 4 |
| 6 | T2 NACA 0012 RANS, Re = 10⁶ (then 6·10⁶) | 4/5 |
| 7 | `ibmSpalartAllmarasIDDES`; T3.1 shielding, T3.2 channel WMLES | 2–6 |
| 8 | moving bodies (nuTilda swept volume as for the scalars), docs (methods doc section, tutorials) | 7 |

## 5. Open questions

1. Target application: the wall-function emphasis (high Re, Cartesian meshes) versus wall-resolved
   IDDES (moderate Re) decides how much of WP5 is needed.
2. Reynolds number of the airfoil case: 10⁶ (cheap, body-fitted comparison only) or 6·10⁶ (TMR and
   experimental data, at least four times the cells)?
3. Moving bodies with turbulence in scope now, or static bodies first?
4. Mesh budget for the 3-D IDDES tests (cores, wall time).
