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

## 5. Decisions (2026-10-30)

1. General-purpose software: both wall functions (wall-modelled) and wall-resolved meshes must
   work, with RANS in the near-wall region in both cases (IDDES in its DDES/WMLES branches, never
   wall-resolved LES). Both `wallFunction spalding` and `none` are first-class; WP5 (image point,
   tangential-only ν_w) is in scope if T1/T2 show cut-cell sensitivity.
2. Airfoil at Re = 10⁶, compared with our own body-fitted runs only (no TMR/experiment needed).
3. Static bodies first, no scalars. Keep the design open for both: nuTilda assembled with
   `fvm::ddt(ibm.alpha(), …)` so that the swept volume works once bodies move, and the IBM pieces
   factored so that scalar transport can reuse the wall functions (wall heat/mass transfer) later.
4. Smallest possible cases first (2-D, ≤ 10⁵ cells, workstation); larger 3-D IDDES runs later on
   the workstations (istmio2, istmeridanus, istmhydra) or HoreKa (partition `cpuonly`, account
   `hk-project-exasim`).

## 6. Why the 2-D SA airfoil does not validate IDDES completely

SA-IDDES solves the same nuTilda equation as SA; it differs only in the length scale d̃, which
blends the wall distance y with the grid scale through the IDDES functions (f_d shielding, f_e,
f_B, with r_d and r_dt built from y and ∇U), and it runs transient. The IBM pieces (operators,
Dirichlet wall, wall function, wall distance) are shared, so a working SA airfoil validates most of
the machinery. Three things are not covered and are checked by T3.1 (cheap, 2-D/thin slab):

- y enters IDDES far from the wall as well (blending functions); the IBM distance must be right
  everywhere, not only in the wall cells, and so must h_max at cut cells;
- the grid scale: Cartesian IBM meshes have near-isotropic cells at the wall, whereas body-fitted
  boundary-layer meshes are stretched along the wall. With small isotropic cells the LES branch can
  switch on inside an attached boundary layer (modelled-stress depletion, grid-induced separation)
  where the body-fitted run stays RANS. Whether the f_d shielding keeps the near-wall region in RANS
  on IBM meshes is the IDDES-specific question; the `LESRegion` and blending fields are written and
  compared with the body-fitted run;
- the transient path (`ibmPimpleFoam`, ddtCorr) with a variable ν_t.

Real LES content (resolved turbulence) needs 3-D and is tested by T3.2 (channel, WMLES branch).

## 7. Status (2026-10-30)

Done (WP1–WP3): `wallDist { method ibmLevelSet; }` (distance at the fluid centroids, exact to
1e-12 in full cells for planar walls; the centre of a cut cell may lie in the solid, its centroid
not), `cutCellGeometry::fluidGrad` (gradient over the fluid with a wall value),
`libcutCellIBMTurbulence` with `ibmSpalartAllmaras`, per-body `wallFunction spalding | none` and
`cutCellGeometry::wallViscosity` (Spalding by bisection in u⁺, ν_w = ν Re_y/u⁺²), solver changes
(wall viscosity in the traction and the forces, transpose stress for turbulent runs). Laminar
tutorials (steady, transient, moving, rotating, second-order Taylor–Couette, scalar) are
byte-identical to `scalar_secondOrder`.

T1 channel, Re_τ ≈ 550 (`validation/data/channelSA/results.txt`), u_τ against body-fitted:

| mesh | IBM first order | IBM second order |
|---|---|---|
| wall-resolved, h₁ = 1.8e-3, f = 0.3 and 0.05 | −1.14 %, −1.09 % | −0.015 % |
| wall-resolved, h₁ = 0.9e-3, f = 0.3 | −0.57 % (first-order convergence) | +0.015 % (= body-fitted fine to 6 digits) |
| wall functions, y⁺₁ ≈ 47, f = 0.5 / 0.3 / 0.15 / 0.05 | −1.0 / −1.5 / +0.8 / +6.8 % | +1.8 / +2.8 / +7.5 / +24 % |

Wall-resolved: correct and robust to the cut position; `none` and `spalding` coincide at y⁺ < 1.
Wall functions: sensitive to the cut position, as anticipated in 2.3: a thin cut cell (y⁺ ≈ 2)
next to a full cell at y⁺ ≈ 60 carries the wall law at its own sample point, and the face shear
between the two cells is not the wall shear; the centroid-consistent fluxes of the second-order
treatment make it worse (they extrapolate a near-linear profile across a log layer).

WP5, option C (cell linking), done: `mergeAlpha a` in `ibmProperties` links every cut cell with
α < a to a fluid face neighbour (chains end at a cell with more fluid). Linking keeps the matrix
sparsity: a conductance c = mergePenalty·max(diagonals) (default 10⁴) on the link face of every
assembled equation (momentum, pressure, nuTilda), so that the sum of the two rows is the merged
cell's equation and the values coincide. A group (master and slaves) is then treated as one control
volume wherever the equations are nonlinear or geometric: wall distance (also for the turbulence
model) and wall coefficients from the merged centroid, wall function sampled at the master, SA
sources evaluated once at the master with the group's fluid volume (`mergedAlpha`), and the
gradients (∇U, ∇ν̃) averaged over the group. Linking the values alone was not enough (−2.5 % with
the sources at the slave's own distance and gradient; ν̃ 5× too small in the merged cell).

| u_τ vs body-fitted, f = 0.3 / 0.15 / 0.05 | first order | second order |
|---|---|---|
| wall functions, unlinked | −1.5 / +0.8 / +6.8 % | +2.8 / +7.5 / +24 % |
| wall functions, linked (mergeAlpha 0.5) | +0.08 / +0.06 / +0.17 % | +0.44 / +0.22 / +0.09 % |
| wall-resolved, linked | +0.06 / – / +0.01 % | +0.37 / – / +0.05 % |

Known issue: the wall-resolved first-order case with f = 0.05 does not converge fully (G wanders
by ±0.04 %, independent of the penalty and of the inner solver tolerances). Not yet linked: the
passive scalars, and links across processor boundaries (a cell whose only candidate master is on
another processor stays unlinked; counted in the log).

WP5, option A (wall-law face viscosity, no linking) for comparison: `wallFunctionAlpha a` gives a
wall cell with α < a a partner (its outward fluid neighbour, the same rule as a link master); u_τ
comes from Spalding at the partner's centroid, the thin cell's traction imposes u_τ², and the face
between the two carries ν_f = u_τ Δ/(u⁺(y⁺_partner) − u⁺(y⁺_cell)), blended by (n_f·n_w)². Only the
momentum equation changes; the SA equation of the thin cell is left as it is.

| u_τ vs body-fitted, f = 0.3 / 0.15 / 0.05 | option A | option C |
|---|---|---|
| wall functions, first order | −2.0 / −1.6 / −1.1 % | +0.08 / +0.06 / +0.17 % |
| wall functions, second order | −1.6 / −1.5 / −1.0 % | +0.44 / +0.22 / +0.09 % |
| wall-resolved, first order | −0.35 / – / −0.07 % (not converged) | +0.06 / – / +0.01 % |
| wall-resolved, second order | −0.01 / – / −0.02 % (not converged) | +0.37 / – / +0.05 % |

A removes most of the cut-position sensitivity of the unlinked scheme (spread 1 % instead of 8–22 %)
but keeps a bias of −1 to −2 % and a spread of 1 %, and its wall-resolved runs settle into limit
cycles (G spread up to 5.5e-4). C is about ten times more accurate and less sensitive. The likely
reason is that A leaves the thin cell's own SA equation (sources at its tiny distance and large
gradient) and the interpolated viscosities of its other faces untouched; fixing those would turn
A into C. Option A was removed again (commit after 3bc5a33); cell linking is the treatment.
