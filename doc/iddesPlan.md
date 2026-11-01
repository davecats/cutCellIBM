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
  against +20 % with wall functions: the wall-function excess is the wall model's, not the
  nose geometry's.
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
