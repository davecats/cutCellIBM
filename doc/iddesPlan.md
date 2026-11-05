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

*Current state (2026-11-01, section 8): the SA analytic law of the wall instead of Spalding
(`wallFunction spalartAllmaras | none`), read from the cell `wallSampleDistance` = 1.5 h out along
the wall normal and blended in by h⁺, ν_w capped at ν + κu_τh, ν̃ = κu_τy in the modelled wall
cells; Spalding and the other variants were removed after the study. The original design follows.*

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
volume wherever the equations are nonlinear or geometric: every cell of a group gets the wall
distance of the merged centroid (wall coefficients, wall function, turbulence wall distance), and
the gradients (∇U, ∇ν̃) are averaged over the group. With its common distance, gradients and linked
values, every cell of a group has the source density of the merged cell, so the SA sources need no
special weighting. Linking the values alone was not enough (−2.5 % with the sources at the slave's
own distance and gradient; ν̃ 5× too small in the merged cell). (A first version evaluated the
sources only at the master with a group-volume weight; that weight only hid a bug in the slaves'
wall distance and was removed in 8a2d288.)

| u_τ vs body-fitted, f = 0.3 / 0.15 / 0.05 | first order | second order |
|---|---|---|
| wall functions, unlinked | −1.5 / +0.8 / +6.8 % | +2.8 / +7.5 / +24 % |
| wall functions, linked (mergeAlpha 0.5) | +0.08 / +0.06 / +0.17 % | +0.44 / +0.22 / +0.09 % |
| wall-resolved, linked (oscillation means) | +0.07 / – / +0.02 % | +0.38 / – / +0.06 % |

On the wall-resolved mesh (first cell 1.8e-3) every run, body-fitted `simpleFoam` included, ends in
a period-2 oscillation of SIMPLE (G amplitude 6.5e-5 body-fitted, up to 7e-4 IBM); wall-resolved
values are means over it (`results.txt`, correction section). Finer meshes and all wall-function
runs converge to round-off. Pressure linking had no effect in the channel (uniform pressure).

Pressure-velocity coupling of the merged cells (found on the airfoil): rAU = 1/A taken from the
linked momentum matrix holds the penalty, on the master as well as on the slaves, so rAU was ~10⁻⁴
of its value in every linked cell. Along a strip of linked cells (the surface nearly parallel to a
grid line, e.g. near the maximum thickness of the airfoil) the pressure Laplacian and the Rhie-Chow
damping between neighbouring groups vanished and a long-wave mode appeared (C_p ±1.5 over 20 cells,
reaching 10 cells into the flow; SIMPLE did not converge). Now the penalty enters only the momentum
predictor solve (on a copy of the matrix); pEqn takes A, H (and H1 for SIMPLEC) of the merged
control volume from the unlinked matrix, A_g = (Σ V_i A_i + Σ_l (upper_l + lower_l))/Σ V_i over the
cells and internal link faces of a group, H_g likewise with the coupling through the internal faces
removed (`mergeA`, `mergeH`), and the momentum corrector the V-weighted group mean of ibm.grad(p)
(`groupMean`). The channel is unchanged (u_τ to 2·10⁻⁴), laminar runs byte-identical, the airfoil
converges to round-off without the mode. Not yet linked: the passive scalars, and links across processor boundaries (a cell whose only candidate master is on
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
| wall-resolved, first order | −0.35 / – / −0.07 % | +0.07 / – / +0.02 % |
| wall-resolved, second order | −0.01 / – / −0.02 % | +0.38 / – / +0.06 % |

A removes most of the cut-position sensitivity of the unlinked scheme (spread 1 % instead of 8–22 %)
but keeps a bias of −1 to −2 % and a spread of 1 % with wall functions. C is about ten times more accurate and less sensitive. The likely
reason is that A leaves the thin cell's own SA equation (sources at its tiny distance and large
gradient) and the interpolated viscosities of its other faces untouched; fixing those would turn
A into C. Option A was removed again (commit after 3bc5a33); cell linking is the treatment.

T2 NACA 0012, Re = 10⁶ (`validation/data/naca0012/results.txt`). References: body-fitted
`simpleFoam` on the NASA TMR grid n0012_449-129 (far field 500 c), wall-resolved, and with wall
functions (wall grid lines removed up to a first cell of 2e-3). IBM: Cartesian core h = 2e-3
(y⁺ ≈ 50), R = 500, linking mergeAlpha 0.5; 10⁴ SIMPLE iterations, all converged to round-off
(dC_l/500 < 10⁻⁵).

| | C_l | C_d (pressure / viscous) |
|---|---|---|
| α = 0, body-fitted wall-resolved | 0 | 0.01105 (0.00209 / 0.00897) |
| α = 0, body-fitted wall functions | 0 | 0.01003 (0.00043 / 0.00959) |
| α = 0, IBM first / second order | 0 | 0.01274 / 0.01282 (0.0024 / 0.0104) |
| α = 10, body-fitted wall-resolved | 1.0477 | 0.01585 (0.00812 / 0.00773) |
| α = 10, body-fitted wall functions | 0.51 (stalled, unsteady) | 0.146 |
| α = 10, IBM first / second order | 0.9727 / 0.9839 (−7.2 / −6.1 %) | 0.02740 / 0.02603 (0.020 / 0.0074) |

- Linking on the airfoil exposed the pressure-velocity defect described under WP5 above (C_p ±1.5
  along strips of linked cells, no convergence); fixed in 07fca55.
- C_p of the IBM runs matches the body-fitted ones (rms 0.011 at α = 0; suction peak at α = 10
  −5.13 vs −5.39). IBM wall-function runs stay attached at α = 10, the body-fitted wall-function
  reference does not (also on the blockMesh C-mesh), so at α = 10 only the wall-resolved reference
  is usable.
- C_f per wall cell (the applied wall shear) has a sawtooth of about 10 % following the staircase of
  wall distances; the field is smooth: an image-point estimate (U, p at a fixed distance D along
  the normal, Spalding at D; `naca0012_surface.py --image D`) has 1 % noise. It is biased where the
  wall law does not hold to D (+9 % on the wall-resolved solution at D = 3e-3, invalid near the
  nose). Compared with that estimator, the IBM C_f at α = 0 lies on the wall-resolved one for
  x > 0.35 (0.00443 vs 0.00462 at mid-chord), the body-fitted wall-function one above (0.00518).
- Wall-function weakness shared by body-fitted and IBM runs: in a wall cell of a wall-function mesh
  the Gauss velocity gradient (U = 0 on the wall) is about 3× the log-law dU/dy at the cell centre,
  so the SA production and ν_t there are too large (ν_t = 2.4 κ u_τ y in the body-fitted
  wall-function cell at mid-chord, 0.85 wall-resolved); body-fitted wall functions give C_f +31 %
  at mid-chord and a laminar-like region near the nose, and stall at α = 10. At α = 10 the IBM
  upper-surface C_f is about half the wall-resolved one for x > 0.3 and the pressure drag 2.5×.

Open, candidates before IDDES: (i) log-law-consistent SA production in wall-function wall cells
(Ω or ν̃ of the wall cell from u_τ/(κ y)); (ii) wall-model input at a fixed matching distance
(2–3 h, image point) instead of the wall cell: a smooth applied wall shear independent of the cut
position, and the input WMLES needs for IDDES anyway (non-local stencil, processor halos, a
fallback where the boundary layer is thinner than the matching distance).

## 8. Wall-treatment study and clean-up plan (2026-11-01)

Question: is the way the IBM wall cells carry the SA wall law the best one? Tested on an inclined
channel (Re_τ ≈ 550, walls at tan α = 0, 1/8, 1/4, 1/2, 1, three lattice offsets each, uniform
cells at h⁺ ≈ 1, 10, 31, 93, initialised from the body-fitted profile) and on the NACA 0012 at
α = 3° (wall functions at h = 2e-3 and 1e-3; wall-resolved with `--levels 6`, h = 3.1e-5, 1.39M
cells). References: body-fitted wall-resolved, and body-fitted wall functions at the same h.

Variants of the wall cell's u_τ (all with linking unless marked NL):
`old` Spalding at the own fluid-centroid distance, no ν̃ wall law; `centroidSA` SA law
(Allmaras et al. 2012) at the centroid distance, ν̃ = κ u_τ y in the wall cells; `groupCentroid`
the same at the merged-group centroid; `halfCell` at h/2; `cellAverage` the law averaged over the
cut cell's fluid (16 wall-distance quantiles); `sample` the velocity of the cell 1.5 h out along
the normal, blended in with the weight w(h⁺) (smoothstep in log h⁺ between 1 and 10).

Results (u_τ error, means over offsets; NACA C_d against the wall-resolved body-fitted 0.01115):

| | channel h⁺≈10 aligned / inclined | h⁺≈31 | h⁺≈93 | NACA C_d h = 2e-3 / 1e-3 | C_f rms (1e-3) |
|---|---|---|---|---|---|
| old | +9.6 / +21..+29 | +8.4 / +10..+20 | +4.4 / −0.1..+10 | +21.0 / +20.5 % | 24 % |
| centroidSA | +8.5 / +21..+28 | +4.1 / +6..+17 | +0.9 / −2.5..+8.8 | +15.9 / +12.9 % | 18 % |
| groupCentroid | +8.5 / +21..+31 | +4.1 / +7..+18 | +0.9 / −0.7..+9.9 | +16.2 / +13.2 % | 18 % |
| halfCell | +12.3 / +26..+38 | +4.2 / +8..+23 | +0.4 / −1.0..+11 | +17.6 / +13.4 % | 18 % |
| cellAverage | +10.3 / +22..+32 | +11.2 / +15..+27 | +6.7 / +7..+17 | +26.5 / +26.1 % | 30 % |
| sample (capped) | +1.3 / +2.4..+8.6 | +0.4 / −1.0..+3.8 | +0.8 / −0.1..+3.2 | +13.1 / +7.9 % | 8 % |
| body-fitted wall functions | +9.1 | +8.4 | +4.3 | – / +8.6 % | 18–50 % |

- Wall-resolved (h⁺ ≈ 1): every linked variant within 0.8 % in the channel at all angles (the
  wall law has no weight there); unlinked about −1 %. NACA: linked +1.4 % C_d (second order
  +2.0 %), unlinked −1.6 %, C_l within 0.3 %, C_f rms 2 %. Pressure drag +4 % wall-resolved
  against +20 % with wall functions. (That this excess is the wall model's turned out wrong:
  it is the leading edge's resolution by the Cartesian core, section 9.)
- Unlinked runs with wall functions do not converge on inclined walls (G oscillating by a factor
  2–11) and diverge on the airfoil: linking stays.
- The cell-size-based distances (halfCell, cellAverage) do not beat the centroid; groupCentroid
  changes nothing, so the inclined bias of the centroid variants is not the linked-cell distance.
- Two defects found and fixed: (a) with a sample, ν_w = u_τ² d/|u_t| was unbounded as the wall
  cell's own slip vanished; the implicit wall drag then holds the cell at rest, a second steady
  state (+53..+180 % from an initialised field, +1..2 % from rest, crashes). Now capped at
  ν + κ u_τ h, above ν_w on the wall law (u⁺ ≥ ln(1 + κ d⁺)/κ); no effect on converged good
  states (NACA identical to 5 digits). (b) `polyMesh::findCell` is linear in the cells of an
  octree leaf on meshes refined around a small body (> 1 h setup); samples are now located by a
  walk through the wall cell's neighbours, identical results.

Clean-up plan (after committing the study state, so every variant stays reproducible from git):

1. Record the study in `validation/data/{channelSA,naca0012}/results.txt` and here; commit.
2. Remove `wallLawDistance` (centroid/groupCentroid/halfCell/cellAverage, the quantiles,
   `saUTauAverage`, `cellAverageUTau`), the `nuTildaWallLaw` switch and the Spalding wall
   function (`wallFunction spalartAllmaras | none`), and their `makeCase.py` options.
3. Keep: capped `sample` (`wallSampleDistance`, default 1.5) with the w(h⁺) blend, the ν̃ wall
   law with the same weight, linking, the `ibmWallShear` / `ibmWallLawWeight` output. The SA law
   at the wall cell's own centroid stays as the w → 0 end of the blend and as the internal
   fallback of cells whose sample is not available (on another processor, in the solid), not as
   an option. Wall functions require level-set bodies (fatal error otherwise).
4. No separate "none" treatment for wall-resolved walls: at h⁺ ≲ 1 the SA law gives
   ν_w = ν d⁺/u⁺(d⁺) ≈ ν, the law-consistent wall gradient, and switching to plain no-slip would
   need a new blend for 1 < h⁺ < 10 and new baselines. Per-body `wallFunction none` already
   exists for a plain no-slip comparison.
5. Checks: laminar byte-identical; wall-function `sample` channel and NACA h = 1e-3 identical to
   the study runs; wall-resolved channel and NACA unchanged.

Then: the ±1.5 % linked/unlinked spread on wall-resolved meshes, and the wall-function pressure
drag at the nose.

Clean-up done (2026-11-01, commit after bcd8d25): `wallFunction spalartAllmaras` (default) | `none`;
`wallSampleDistance` default 1.5 (0: the wall cell's own value); the wall-law samples are computed
on first use by the solver (laminar runs never compute them) and need level-set bodies (fatal
error otherwise). The cap is ν + κu_τ max(h, dWall): with h alone it cut into the own-cell law of
fallback cells whose fluid centroid lies further from an inclined wall than the cell's normal
extent h = 2V/Σ|S_f·n| (u_τ −0.75 instead of −0.73 % in one h⁺ 93 case). Checks: laminar
byte-identical; channel `sample` cases (h⁺ 1/10/31/93, aligned and inclined) byte-identical to
the study runs; NACA h = 1e-3 and the first 300 iterations of the wall-resolved run equal to
round-off (|ΔU| 3e-7, C_l and C_d to 5 digits; parallel runs).

## 9. Linked/unlinked spread and the wall-function nose drag (2026-11-02)

Runs and tables: `validation/data/naca0012/results.txt` (last block). NACA 0012, α = 3°.

**Linked/unlinked spread (wall-resolved) is a first-order discretisation error.** Refinement
series `--levels 5/6/7` (wall h 6.25e-5 / 3.1e-5 / 1.56e-5, 0.77M / 1.39M / 2.63M cells):

| wall h | linked C_d | unlinked C_d | spread |
|---|---|---|---|
| 6.25e-5 | +5.63 % | −1.69 % | 7.3 % |
| 3.1e-5 | +1.42 % | −1.60 % | 3.0 % |
| 1.56e-5 | +0.72 % | −0.67 % | 1.4 % |

(against the body-fitted TMR 449 grid, C_d 0.011147). The spread halves with h; linked converges
from above (viscous drag 9213 / 8903 / 8850 e-6 onto the body-fitted 8830–8846), unlinked from
below, non-monotone; the mean of the two at L7 is 0.011150. Both parts follow the wall shear:
linked over-predicts it over 0.2 < x < 0.8, unlinked under-predicts it over 0.01 < x < 0.5, and
the pressure part follows (thinner boundary layer, a little more circulation and nose suction);
surface C_p of the two agrees to 0.004. Unlinked second order is the most accurate at L6
(+0.75 %), as in the channel, but unlinked cells diverge with wall functions, so linking stays the
default. `wallFunction none` on L6 equals the SA law to 1e-4 % (w = 0 at h⁺ < 1). The finer
body-fitted grid (TMR 897, 150k iterations: C_d 0.011164, C_l 0.32328; L7 linked +0.57 %, unlinked −0.82 %)
shows the same +3.8 / −4.2 e-4 chordwise pattern (0.2–0.5 / 0.5–0.8) that every IBM run showed
against the 449 grid: that was the 449 reference's error, not the IBM's.

**The wall-function nose excess is the leading edge's resolution, not the wall model.**
- Second order (wall pressure extrapolated with lsqGrad(p)): C_d +8.06 / +13.04 % against
  +7.94 / +13.10 % first order (h = 1e-3 / 2e-3). Solver and CV force agree to 1e-6: the excess
  is in the solution, not in the force integral.
- Free slip, laminar, on the same Cartesian meshes (no boundary layer, no wall law): spurious
  pressure drag at x < 0.002 of +15.6 / +8.6 / +2.3 e-4 at h = 2e-3 / 1e-3 / 5e-4 against
  2.5e-4, total +10.7 / +3.8 / +0.6 e-4. The wall-function runs' excess over the wall-resolved
  IBM run is +10.2 / +5.0 e-4 at x < 0.002 and +8.7 / +3.6 e-4 in total, about the same.
  The NACA 0012 nose radius is 0.016 c: 8 cells at h = 2e-3.
- Two refinement levels at the leading edge only (x < 0.05, to 5e-4, +2.5k cells): the nose
  excess at x < 0.01 drops from +6.7 / +2.2 to +0.2 / +0.5 e-4; C_d +13.1 → +10.9 %
  (h = 2e-3), +7.9 → +7.3 % (h = 1e-3).
- What remains with wall functions is mostly viscous: +4.6 / +7.8 e-4 (+5 / +9 % of the
  viscous drag, 0.05 < x < 0.5), the equilibrium wall law in the thin front boundary layer (the
  body-fitted wall functions miss it the other way, −10 %), plus +2–3 e-4 pressure around the
  end of the nose refinement.

Consequences: wall-function meshes need the leading edge resolved for the inviscid flow
(h ≲ r_LE/30, here 5e-4), independent of the wall treatment; this is a meshing rule, not a code
change (a local level by x, or better a curvature-based refinement criterion, in the generator).
The wall-resolved spread needs no change: it converges at first order; quote ±0.7 % at L7 and
±1.5 % at L6.

## 10. SIMPLEC defaults and second-order linked cells (2026-11-02)

Runs and tables: `validation/data/naca0012/results.txt` (last block) and
`validation/data/secondOrder/taylorCouette_linked.txt`.

**SIMPLEC.** The slow convergence of the body-fitted references was SIMPLE with p relaxation 0.3,
not non-orthogonality (1 or 2 correctors change nothing). With SIMPLEC (`consistent yes`, p 1,
U/ν̃ 0.9) the TMR grids converge: 449 C_d 0.011137, 897 C_d 0.011134, C_l 0.32319 (100k
iterations, drift −5e-8 per 1000); the SIMPLE 897 value 0.011164 was 0.27 % high. 897 SIMPLEC is
the reference from here on. IBM: wall functions converge with U/ν̃ 0.9 (1k iterations instead of
4k); the linked wall-resolved run ends in a limit cycle with 0.9 (C_d ±0.15 %) and converges with
U/ν̃ 0.7 to the SIMPLE value (0.011305) in a third of the iterations. The NACA generator now writes
SIMPLEC with p 1, U/ν̃ 0.7.

**Why linked cells were first order.** The group shares one value, but the outer-face fluxes
used each cell's own centroid as the point that value belongs to. Fix, active only with
`secondOrder` and `mergeAlpha > 0` (`cutCellGeometry::valueOffset_`):
- the value point of every group cell is the group centroid (α V-weighted fluid centroids);
  outer-face distances, face-value offsets and the wall-pressure moments use it; faces inside a
  group keep their own centroids (their fluxes cancel in the merged equation);
- one least-squares gradient per group over its outer faces (per-cell fits without the
  intra-group faces were near-singular: a slave's remaining neighbours lie along the wall);
- a master that is itself cut also takes the group wall distance for its wall coefficients
  (only slaves did before);
- the moving-wall velocity is taken at the foot of the normal from the group centroid.
Unlinked second order and linked first order are byte-identical to before.

Taylor–Couette, second order, bulk L2 velocity error (N = 40 / 80 / 160 / 320):

| | N 40 | N 80 | N 160 | N 320 |
|---|---|---|---|---|
| unlinked | 1.10e-3 | 2.76e-4 (2.00) | 7.25e-5 (1.93) | 1.78e-5 (2.03) |
| linked, before | 3.89e-3 | 2.16e-3 (0.85) | 1.30e-3 (0.73) | 6.01e-4 (1.11) |
| linked, now | 1.20e-3 | 3.21e-4 (1.90) | 7.43e-5 (2.11) | 3.25e-5 (1.19) |

Open: at N 320 the linked error near the rotating inner wall (5.9e-5 at r 0.65 against 2.3e-5
unlinked) decays into the gap; mid-gap both are 1.2e-5. Something in the moving-wall treatment
of groups is still first order; it does not affect static walls. (Resolved in §11: it was the
continuity flux, and it is not specific to moving walls.)

NACA 0012, α = 3°, wall-resolved, C_d against 897 SIMPLEC:

| level (wall h) | linked 1st | unlinked 1st | linked 2nd | unlinked 2nd |
|---|---|---|---|---|
| L5 (6.25e-5) | +5.75 % | −1.58 % | +5.50 % | +2.88 % |
| L6 (3.1e-5) | +1.53 % | −1.49 % | +1.35 % | +0.87 % |
| L7 (1.56e-5) | +0.83 % | −0.56 % | +0.74 % | +0.65 % |

At second order the linked/unlinked spread drops by about 5 per level (2.6 / 0.48 / 0.09 %), and
both series approach about +0.65 %; C_l is +0.25 %. That remainder is not the wall discretisation:
unlinked pressure drag is the same at L6 and L7 (+1.0 %), viscous drag +0.4–0.6 % changes little
from L6 to L7. Its source (outer Cartesian resolution, far field, wall distance in SA) is not
identified; the 449/897 body-fitted pair differs by 0.03 % in C_d and 0.36 % in C_l. On these grids
the linking change itself moves C_d by −0.25 / −0.2 / −0.1 points: the linked first-order error on
the airfoil was mostly resolution, the gain is the vanishing linked/unlinked spread.

## 11. Second-order continuity flux (2026-11-02)

Details and tables: `doc/secondOrderStudy.md`, "Fix of 2026-11-02".

**Cause of the remaining linked error.** `ibm.flux()` (continuity flux of `phiHbyA`, initial φ,
`ddtCorr`) was the linear interpolate at m_f, without the g_f·(x'_f − m_f) correction that every
other face value has. A first-order flux on every cut face; for linked cells, whose value point
is far from the group's faces, large and systematic. Diagnosis: with the outer cylinder rotating
the excess appeared at both walls (not the moving wall); at Ω = 10⁻³ all errors scaled by 10⁻³
(a Stokes error, not convection); scalar Dirichlet tests with linking were second order
(`ibmOrderTest` now links). Ruled out first, without effect: limiting the moving-wall point to the
group box, the level-set distance of the group centroid, a quadratic group fit.

Taylor–Couette at N 320: bulk 1.94e-5 (unlinked) / 1.80e-5 (linked), both of order 2.1; cut cells
1.7e-5 / 2.7e-5 (were 9.6e-5 / 1.5e-4); the wall gradient and the cut-cell pressure now converge
(pressure L2 3.5e-3 → 5.1e-4 unlinked, 7.1e-3 → 4.5e-4 linked). Parallel = serial to 1.4e-8.
All second-order results change; first order is unchanged.

**Still first order: the cut-cell pressure, hence the forces.** Cut-cell pressures scatter by
about h·∂p/∂r (a viscous error: it scales with Ω, not Ω²), only their mean converges, at first
order. Cylinder drag: pressure part orders 0.8–0.95, free slip 1.0–1.1 once sliver absorption is
switched off (α_min = 0.02 absorbs eight cells at N = 161 and 321 only, which moves the free-slip
drag by up to 0.9 %). Remedies: an O(h)-truncation (quadratic) viscous flux in the cut cells
(deferred 2026-11-03, too complex for now), or a wall pressure fitted on full cells (pursued).

NACA 0012 wall-resolved, second order, C_d against 897 SIMPLEC (wr5, SIMPLE as wr4):

| level | linked, before → now | unlinked, before → now |
|---|---|---|
| L5 | +5.50 → +5.28 % | +2.88 → +2.76 % |
| L6 | +1.35 → +1.26 % | +0.87 → +0.83 % |
| L7 | +0.74 → +0.70 % | +0.65 → +0.64 % |

The flux fix is a small gain on the airfoil: the L5 spread is wall resolution, not the flux. At
L7 the linked/unlinked spread is 0.06 %, the viscous drag is within 0.03 % of 897 (was +0.4–0.6 %),
and the remaining +0.64 % is pressure drag (+2.9–3.5 %, up from +1.0–2.0 %): the first-order
cut-cell pressure above, now no longer partly cancelled by the viscous error. C_l +0.25 %.

## 12. Linking by default, SA-neg, parallel (2026-11-03)

Linking is the default and independent of the decomposition (`doc/linkingPlan.md`, Status;
`cutCellIBM.tex`, "Cell linking"). Consequences for the turbulence cases:

**Flow-rate forcing with linked cells** (8aebbce). `ibmMeanVelocityForce` corrected every cell
by rAU·α·ΔG. That correction is not uniform within a group, so it broke the group-uniform
velocity that the pressure corrector produces. The 45° wall-function channel then had a U
residual floor of 1e-9 and GAMG hit maxIter every iteration. The response is now the group mean
(floor 2e-12; cylinderRe1 linked: same drag, 1565 instead of 2145 iterations).

**Inclined wall-function channel, rerun** (tables: `runs/ibmSA_study/linkC/inc/table.txt`).
Unlinked runs are unchanged. Linked runs have errors of the same size as before (n56,
tan α = 1/4: u_τ +7.06 % against +6.26 %). The conclusions of §8 stand.

**SA-neg** (e99057b; Allmaras, Johnson & Spalart, ICCFD7-1902, 2012), switch
`ibmSpalartAllmarasCoeffs { negativeNuTilda true; }` (default).
- Problem: in the NACA 0012 wall-function case, ν̃ undershot in about 95 full cells at the edge
  of the far wake (x = 3–33 c) every iteration. bound() clipped them, so the ν̃ residual
  stalled at 2e-5 and the run never reached a fixed point. Serial and parallel could only
  agree to the drift (1.4e-4 in drag).
- Where ν̃ < 0:
  - production C_b1 (1 − C_t3) Ω ν̃;
  - destruction −C_w1 (ν̃/d)²;
  - diffusivity (ν + f_n ν̃)/σ, with f_n = (C_n1 + χ³)/(C_n1 − χ³);
  - ν_t = 0, and no clipping.
- Where ν̃ ≥ 0, every term is the standard one to the bit.
- Effect on existing cases: a channel without undershoot is byte-identical with the switch on or
  off; a channel that clipped once differs by 1e-11.
- NACA h2 wall functions, linked, 175k cells (`runs/ibmSA_study/linkC/nacaneg`): converges by
  residual control (p 1e-7) in about 3600 iterations; ν̃ residual 4e-10.
  - Serial / scotch 8 / scotch 16: F_y 1.617871e-3 / 1.617842e-3 / 1.617849e-3; F_x
    −2.16692e-5 / −2.16622e-5 / −2.16680e-5.
  - Restarted with residual control 1e-11/1e-12 and run to 20000 iterations, the gap shrinks
    (F_y 2.8e-7 / 1.1e-6, F_x 4e-6 relative) and is still closing in the parallel runs.
  - Why the gap: one iteration from the same state, with all equations solved to round-off
    (pressure by PCG to 1e-17), gives the same p to 4e-11 and U to 1e-10. So the
    discretisation is independent of the decomposition, and the gap is convergence. A link over
    a processor face is a stiff mode that parallel GAMG reduces slowly, and with an absolute
    p tolerance (1e-10) the pressure solve stops before it is resolved.
  - A real decomposition dependence found on the way and fixed (b81d779): the solid-row
    diagonal of the pressure matrix (only the residual normalisation) lacked the coupled-face
    coefficients.

**NACA 0012 wall-resolved with linking by default** (second order, C_d against 897 SIMPLEC; runs
`naca3/wrC`, `naca3/wrN`; data in `validation/data/naca0012/results.txt`, last block):

| level | linked before (wr5) | linked default | unlinked (wr5) |
|---|---|---|---|
| L5 | +5.28 % | +5.41 % | +2.76 % |
| L6 | +1.26 % | +1.27 % | +0.83 % |
| L7 | +0.70 % | +0.70 % | +0.64 % |

- The new linking leaves L6 and L7 unchanged; L5 is 0.13 points worse. C_l is unchanged.
- All wall-resolved runs, here and in §9–11, clipped ν̃ in about 97 % of their iterations, and
  their ν̃ residual stalled at 1e-5. The convergence criterion used there (ΔC_d over 500
  iterations) did not show this.
- Restarted with SA-neg, L5–L7 meet residual control (p 1e-7, U/ν̃ 1e-8) after 65–78
  iterations, with the same coefficients to the digits shown. The clipping blocked the fixed
  point but did not bias the forces, so the earlier tables stand.

## 13. The linked wall-resolved excess: convective face values (2026-11-04)

Details: `doc/secondOrderStudy.md`, "Convective face values of cut and linked cells";
data in `validation/data/naca0012/results.txt` (last block).

- The 0.13 points that L5 lost with linking by default (+5.28 → +5.41 %) come from no longer
  absorbing slivers (alphaMin 0.02 restores +5.28 %).
- The wall law is not the cause (switched off: +5.84 / +2.83 %).
- The linked excess is viscous drag at x < 0.2 plus pressure drag. The two-point wall formula
  is exact on the body-fitted solution, so the error is in the solution.
- Cause: the convective face value (linearUpwind, U and ν̃) is extrapolated from the cell
  centre, but a slave carries the group value from about 0.5–0.7 h off the wall.
- `ibmConvection true` (the default since 773df1b) extrapolates from the value point to the wet-face
  centroid. C_d linked +1.13 / +0.47 / +0.61 %, unlinked +2.16 / +0.73 / +0.63 % at L5/L6/L7
  (default: +5.41 / +1.27 / +0.70 and +2.76 / +0.83 / +0.64 %).
- The remainder at L6/L7 is pressure drag (+3 %), independent of the level and of the linking.
- The residual-control convergence in §12 is weak: restarts stop after a few iterations.
