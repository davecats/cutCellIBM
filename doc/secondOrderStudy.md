# Towards a second-order cut-cell wall: an order-of-accuracy study

Exploration, without changes to the OpenFOAM code, of what limits the order of the
cut-cell immersed boundary for no-slip (Dirichlet) and free-slip (Neumann) walls, and of
the smallest set of ingredients that makes both second order. Everything is measured on
the notebook's discretisation with `validation/order_study.py` (2-D Cartesian, periodic
box, level-set body, sub-sampled fractions, sliver threshold 0.02), on manufactured
solutions with the fluid being the interior of a disk of radius 1:

| test | equation | wall | exact solution |
|---|---|---|---|
| LD | ∇²w = 1 | w = 0 | w = −(R² − r²)/4 |
| LN | ∇²w = 16r² − 8R² | ∂w/∂n = 0 | w = r⁴ − 2R²r² |
| CD | u·∇T − ∇²T = q, u = Ω(−y, x) | T = 0 | T = (R² − r²) x y |
| CN | same, Ω = 10 | ∂T/∂n = 0 | T = (R² − r²)² x y |

Solid-body rotation is tangential to the wall and divergence-free, so the convection
tests exercise the wall treatment with convection without a flow solve; the cell Péclet
number at N = 20 is 1.5. Errors are α-weighted L2 norms over the fluid cells against the
exact solution at the fluid centroid of each cell; Linf is reported in the script output.

Variants (all deferred corrections are iterated to convergence here; in the solver they
would be one correction per outer iteration):

- wall term, Dirichlet: **V0** current, d = αV/(2A_wall); **V1** d = distance of the fluid
  centroid to the wall (from the level set, or |x_c − x_wall|); **V2** V1 plus a curvature
  correction of the two-point gradient taken from the cell's own Laplacian (turns into a
  scaling of the diffusive row); **V3** V1 plus a three-point quadratic along the wall
  normal with one bilinearly interpolated probe at d + h;
- diffusive face flux: **F0** current, θA(T_N − T_P)/h; **F1** consistent with the fluid
  centroids: implicit part θA(T_N − T_P)/(n·d) with d the centroid-to-centroid vector, explicit
  remainder g_f·n − (g_f·d)/(n·d) with least-squares gradients on the centroids (the
  over-relaxed non-orthogonal form; exact for a linear field);
- convective face value: **C0** current, ½(T_P + T_N); **C1** value reconstructed at the
  wet-face centroid from the two cells' linear reconstructions (deferred).

## Results

L2 error and observed order, N = 20 … 320 (h = 3/N):

**LD, Laplacian with Dirichlet wall**

| variant | 20 | 40 | 80 | 160 | 320 |
|---|---|---|---|---|---|
| V0/F0 (current) | 1.79e-2 | 1.14e-2 (0.65) | 4.92e-3 (1.22) | 2.40e-3 (1.04) | 1.23e-3 (0.96) |
| V1/F0 | 1.65e-2 | 9.45e-3 (0.80) | 4.09e-3 (1.21) | 2.14e-3 (0.93) | |
| **V1/F1** | 1.55e-3 | 4.09e-4 (1.92) | 9.81e-5 (2.06) | 2.31e-5 (2.09) | 5.52e-6 (2.07) |
| V2/F1 | 2.49e-3 | 6.58e-4 (1.92) | 1.63e-4 (2.01) | 5.18e-5 (1.66) | |
| V3/F1 | 2.32e-3 | 6.00e-4 (1.95) | 1.52e-4 (1.98) | 3.48e-5 (2.13) | |

**CD, convection–diffusion with Dirichlet wall**

| variant | 20 | 40 | 80 | 160 | 320 |
|---|---|---|---|---|---|
| V0/C0/F0 (current) | 1.13e-2 | 7.71e-3 (0.55) | 3.76e-3 (1.04) | 1.62e-3 (1.21) | 9.45e-4 (0.78) |
| **V1/C1/F1** | 3.46e-3 | 1.06e-3 (1.70) | 2.78e-4 (1.94) | 6.97e-5 (1.99) | 1.77e-5 (1.98) |
| V3/C0/F1 | 5.86e-3 | 1.79e-3 (1.71) | 4.77e-4 (1.91) | 1.17e-4 (2.02) | |
| V3/C1/F1 | 5.45e-3 | 1.71e-3 (1.68) | 4.55e-4 (1.91) | 1.12e-4 (2.03) | |

**LN, Laplacian with Neumann wall** (source projected onto the compatible space, see below)

| variant | 20 | 40 | 80 | 160 | 320 |
|---|---|---|---|---|---|
| F0 (current) | 4.09e-3 | 8.45e-4 (2.27) | 1.56e-4 (2.43) | 2.01e-5 (2.96) | 5.57e-6 (1.85) |
| F1 | 3.14e-3 | 7.41e-4 (2.08) | 1.48e-4 (2.33) | 2.66e-5 (2.47) | 6.44e-6 (2.05) |

**CN, convection–diffusion with Neumann wall**

| variant | 20 | 40 | 80 | 160 | 320 |
|---|---|---|---|---|---|
| C0/F0 (current) | 3.28e-3 | 8.53e-4 (1.94) | 2.29e-4 (1.90) | 6.45e-5 (1.82) | 1.75e-5 (1.88) |
| C1/F1 | 2.95e-3 | 8.16e-4 (1.85) | 2.22e-4 (1.88) | 6.36e-5 (1.80) | 1.73e-5 (1.88) |

## What the numbers say

1. **The no-slip wall is first order because of the faces, not the wall term.** Improving
   the wall distance alone (V1/F0) or the wall gradient alone (V3 on top of F0, not shown,
   identical to V1/F0) leaves the order at ≈1. The dominant error is the diffusive flux
   through the faces of a cut cell: (T_N − T_P)/h assumes the two cell values sit h apart
   along the face normal, whereas they represent the fluid parts, whose centroids are
   neither h apart nor aligned with the normal. Making the face flux consistent with the
   centroids (F1) together with the centroid wall distance (V1) gives clean second order,
   with errors 200 times smaller than the current scheme at N = 320. The same holds with
   convection (CD).

2. **The wall gradient needs no probe.** With F1 in place, the two-point wall flux
   ν(T_P − T_w)/d_c with the centroid distance is enough for second order; the three-point
   quadratic (V3) is not better. The reason is that the two-point difference with the
   fluid-average value is a second-order estimate of the flux through the wall segment of
   a cell whose interior is resolved to second order by its other fluxes; the first-order
   argument (the gradient "sits at d/2") applies to point values, not to cell averages
   coupled to consistent neighbours. This is the important simplification: no stencil
   beyond the face neighbours, no interpolation into the fluid, no ghost values.

3. **The row-scaling correction (V2) is inconsistent** (order drops to 1.7 and Linf to 1.1):
   the cell Laplacian contains tangential curvature and is not the wall-normal second
   derivative. Do not use it.

4. **The free-slip wall is already second order.** On the Neumann tests the current
   scheme gives orders 1.9–2.4 in L2 (1.7–1.9 in Linf), consistent with the paper, and F1
   changes nothing. The apparent orders of 1.5–1.8 obtained before projecting the source
   were an artefact of the test, not of the wall: absorbing slivers removes fluid, the
   discrete source of the pure-Neumann problem stops being compatible, and the doubled
   diagonal of the reference cell deposits the mismatch as a point source with a
   logarithmic error field. Raising the sliver threshold makes this worse (order 1.4 at
   α_min = 0.1, 1.2 at 0.3). The pressure equation of the flow solver is compatible by
   construction, so this affects only Neumann *tests*; it does show that the boundary
   displacement by absorption, O(α_min h), is what remains at first order for Neumann
   walls, and it argues for keeping α_min small or for cell merging.

5. **Convection needs no special wall treatment.** The wet-face centroid value (C1) is a
   small refinement (10–30 % on the coarsest meshes, nothing asymptotically); the
   ingredients that matter are the same as without convection.

## Recommended path to second order

Minimal set, in order of importance:

1. **Fluid centroids.** Store the centroid of the fluid part of every cell (the
   sub-sampling that gives α gives it for free) and treat every cell value as living
   there: sources, initial conditions, error norms, the wall pressure force
   (p_wall = p_c + ∇p·(x_wall − x_c), untested here but the same argument).
2. **Centroid-consistent diffusive fluxes** across the faces of cut cells (F1): implicit
   coefficient θA/(n·d) and an explicit remainder with least-squares gradients on the
   centroids, applied identically to the momentum Laplacian, the pressure Laplacian and
   the Rhie–Chow face gradient so that flux and velocity corrections stay consistent. In
   OpenFOAM this is exactly the `corrected` Gauss Laplacian evaluated with cut-cell
   `deltaCoeffs` (1/(n·d) on the centroids) and correction vectors, i.e. a small custom
   `laplacianScheme` in the library plus a least-squares gradient on the centroids; the
   stencil does not grow. Uncut cells are untouched (d = h n).
3. **Wall distance from the centroid** (V1): d_c = φ(x_c) for a signed-distance level set
   or the distance to the wall segment; one line.

Optional: wet-face centroid values in the convective flux (C1), as a deferred correction.

Not needed: three-point normal probes, ghost cells, interpolation into the fluid, any
change to closure, blanking, the mass source or the flow-rate forcing.

Stability. The implicit part of F1 is still an M-matrix with a diagonal at least as large
as now (n·d ≤ h in a cut cell is floored at 0.1h); the explicit remainder vanishes for
uniform and linear fields, so free-stream preservation and the closure properties are
retained, and it converged in 11–26 fixed-point iterations in all tests here, which in the
solver becomes one evaluation per outer iteration. Corrections are antisymmetric across
faces, so conservation is exact. For moving bodies the centroids and n·d are recomputed
with the geometry every step.

Cost. One least-squares gradient per corrected field per outer iteration, the centroids,
and one scalar per face.

## What this study does not cover

- The coupled Navier–Stokes problem (pressure equation, Rhie–Chow flux, wall pressure
  force). The Dirichlet/Neumann Laplacian and convection–diffusion results transfer to the
  momentum equation and the scalar directly; the pressure needs the same F1 and the
  extrapolated wall pressure, and should be verified on a manufactured Stokes or
  Navier–Stokes solution (e.g. Wannier or a rotating-cylinder Couette flow).
- Time accuracy (implicit Euler, swept-volume integration): separate, first order.
- Three dimensions: the ingredients are dimension-independent.

A direct verification route in OpenFOAM once implemented: the passive scalar of
`ibmScalarTransport` with a prescribed solid-body-rotation flux and the manufactured
sources above reproduces the CD and CN tests without any flow solve.

Reproduce: `python validation/order_study.py --N 20 40 80 160 320` (about two minutes;
`--noProject` shows the incompatibility artefact, `--alphaMin` its dependence on the
sliver threshold, `--omega` the Péclet number).

## Implementation in OpenFOAM (branch `scalar_secondOrder`) and results

Implemented as recommended, behind `secondOrder true` in `ibmProperties` (default off,
bit-identical results otherwise):

* `cutCellGeometry`: fluid centroids from the sub-sampling, wet-face centroids from the
  polygon clipping, per-face centroid vector and cut delta coefficient 1/(n·d), wall distance
  from the level set at the centroid, a least-squares gradient on the centroids
  (`lsqGrad`), the explicit face-gradient remainder (`snGradCorrection`), wet-centroid face
  values and an extrapolated wall pressure in `grad()` and `forces()`.
* one `snGradScheme`, `cutCellCorrected`, supplying those delta coefficients and the
  remainder; selected for `laplacianSchemes` and `snGradSchemes`, it makes the momentum
  Laplacian, the pressure Laplacian, `pEqn.flux()` and `fvc::snGrad` consistent without any
  change to the solvers.
* `ibmOrderTest`: the manufactured tests of this study on the library operators.

### Operators (`tutorials/order/laplaceDisk/Allrun`)

L2 errors and orders, N = 20/40/80/160 (raw lines in `validation/data/secondOrder/orderTest.txt`),
identical to the Python study to the printed digits for the variants both compute (the
second-order CD/CN rows here are V1/C0/F1: the C1 convective value is not implemented in
OpenFOAM, whereas the CD row of the Python table above is V1/C1/F1):

| test | secondOrder false | secondOrder true |
|---|---|---|
| LD | 1.79e-2, 1.14e-2 (0.65), 4.92e-3 (1.22), 2.40e-3 (1.04) | 1.55e-3, 4.09e-4 (1.92), 9.81e-5 (2.06), 2.31e-5 (2.09) |
| LN | 4.09e-3, 8.45e-4 (2.27), 1.56e-4 (2.43), 2.01e-5 (2.96) | 3.14e-3, 7.41e-4 (2.08), 1.48e-4 (2.33), 2.66e-5 (2.47) |
| CD | 1.13e-2, 7.71e-3 (0.55), 3.76e-3 (1.04), 1.62e-3 (1.21) | 4.06e-3, 1.19e-3 (1.76), 3.11e-4 (1.94), 7.84e-5 (1.99) |
| CN | 3.28e-3, 8.53e-4 (1.94), 2.29e-4 (1.90), 6.45e-5 (1.82) | 3.33e-3, 8.85e-4 (1.91), 2.37e-4 (1.90), 6.67e-5 (1.83) |

### Coupled Navier–Stokes with an exact solution (`tutorials/order/taylorCouette`)

Taylor–Couette flow between an inner cylinder of radius 0.5 rotating at Ω = 1 (no slip) and
an outer cylinder of radius 1.5 at rest (no slip, fluid inside), both immersed, ν = 1,
transient to steady state (`ibmPimpleFoam`, Δt = 0.02, t = 6). Exact solution
u_θ = A r + B/r, p = A²r²/2 + 2AB ln r − B²/(2r²). Errors from
`validation/taylorCouette_error.py` (output in `validation/data/secondOrder/taylorCouette.txt`),
at the fluid centroids for the second-order runs and at the cell centres for the first-order
runs (which do not write `ibmCentroidOffset`). The bulk and cut-cell columns come from an
ad-hoc analysis that was not saved; the script's own subset (cells with α ≥ 0.5) gives
L2(U) orders 1.80 and 1.40 and L2(p) orders 0.35 and 0.39 (second order):

| N | L2(U), first order | L2(U), second order | bulk U rms, 1st / 2nd | bulk p rms, 1st / 2nd | cut-cell p rms, 1st / 2nd |
|---|---|---|---|---|---|
| 40 | 6.48e-3 | 1.16e-3 | | | |
| 80 | 4.48e-3 (0.53) | 3.35e-4 (1.79) | 4.0e-3 / 2.7e-4 | 3.5e-3 / 1.9e-3 | 0.26 / 0.13 |
| 160 | 2.10e-3 (1.09) | 1.27e-4 (1.40) | 2.0e-3 / 8.4e-5 (1.0 / 1.7) | 4.6e-3 / 1.0e-3 (– / 0.9) | 0.20 / 0.11 |

("bulk" = full cells farther than 2h from a wall; the pressure range is 0.08.) The velocity
is second order in the bulk and 5–16 times more accurate everywhere; the overall L2 order
decays to 1.4 at N = 160 because the near-wall full cells stagnate at about 4e-4. The
pressure is second order in the bulk but does not converge in the cut cells (rms 0.11–0.26,
larger than the whole pressure range) and pollutes the near-wall full cells at first order.
This is the small-cell problem of the pressure: the mobility of a cut cell, r_At θ, scales
with αθ because the wall term dominates its momentum diagonal, so its pressure is set by
tiny flux residuals. The velocity is barely affected (the wall term pins it), the forces are.

This test also exposed a defect of the moving-body path: the body velocity of a rotating
wall was evaluated at the cell centre, so u_b·Sw did not vanish for a tangential motion and
an O(h) spurious mass source, alternating in sign around the body, checkerboarded the
cut-cell pressures (errors ±0.4). It was then evaluated at x_w = x_c + d_wall n_w, which only
reduces the defect (no single point can remove it); the source is now computed from the faces,
see "Fixes of 2026-09-29" below, and x_w is used for the wall traction only.

### Forces by Richardson extrapolation (`validation/richardson_drag.sh`)

Periodic cylinder at Re_D = 1 (`tutorials/cylinderRe1`), drag per unit depth from the
solver's force, `ibmSimpleFoam` (at most 8000 iterations), N = 41/81/161/321; the N = 321 runs
and the free-slip N = 161 runs stopped at the iteration limit with residuals of 2e-9 to 2e-8
(final log lines in `validation/data/secondOrder/richardson.txt`; the N = 321 values below
are the final ones and replace earlier snapshots taken during the runs):

| wall | treatment | N = 41 | 81 | 161 | 321 | order (41/81/161) | order (81/161/321) |
|---|---|---|---|---|---|---|---|
| no slip | first | 5.5431 | 5.6604 | 5.7577 | 5.7952 | 0.27 | 1.37 |
| no slip | second | 5.7834 | 5.8205 | 5.8322 | 5.8354 | 1.66 | 1.89 |
| free slip | first | 2.8377 | 2.9421 | 3.0004 | 3.0354 | 0.84 | 0.73 |
| free slip | second | 2.9161 | 2.9848 | 3.0237 | 3.0511 | 0.82 | 0.51 |

The second-order treatment brings the no-slip drag much closer to its limit (the first-order
sequence is still 1.5 % away at N = 161) and its increments shrink faster, but the free-slip
drag, a pure pressure force, converges at first order or below with both treatments. The
no-slip total converges at an apparent order of 1.7–1.9, but largely by cancellation: with
the second-order treatment the pressure drag (2.802, 2.860, 2.890, 2.907) increases at order
about one and the viscous drag (2.982, 2.961, 2.942, 2.929) decreases. Consistently with the Taylor–Couette
diagnosis, the forces are limited by the cut-cell pressure and by the O(α_min h) boundary
displacement of sliver absorption, not by the operators.

### Conclusion and next step

The recommended ingredients deliver second order for the Laplacian and the
convection–diffusion operators with both wall conditions, and for the velocity of the coupled
problem, at the cost of one snGrad scheme and a least-squares gradient per corrected field.
Second-order pressure and forces need one more ingredient that this study did not cover: a
cut-cell pressure that is not determined by the cut cell alone. Candidates, in order of
simplicity: (i) evaluate the wall pressure force from a least-squares fit of the full cells
next to the wall instead of the cut cell's own pressure (post-processing only, no change to
the solution); (ii) cell merging or linking of cut cells with α below ~0.5 to a full
neighbour, which removes the small-cell problem for pressure and momentum alike and would
also relax the moving-wall step restriction.

## Fixes of 2026-09-29: rotating-wall source and linearly exact pressure gradient

The two open issues of `doc/openIssues/openIssues.pdf`, implemented and measured. Raw records in
`validation/data/fixes`.

### 1. Mass source of a rotating wall

`S_P = u_b(x_e)·Sw` with one evaluation point per wall cell cannot give the flux of a rotation
through the wall: `u_b(x_e)·Sw = ω·((x_e − x_m)×Sw)` depends on the offset of x_e across the
normal, so the cell centre, `x_P + d n_w` and the centroid all leave O(|ω| h A_wall) sources.
The source is now

    S_P = U·Sw + S_rot,P,    S_rot,P = −Σ_{f open} θ_f (ω×(x'_f − c))·S_f,

with x'_f the centroid of the geometric wet face (before absorption) and the minimum image
about the current centre c on periodic faces: by the divergence theorem the exact flux of the
rotation through the wall of the post-processed cell (the chord, plus the geometric wet parts of
faces closed by sliver absorption, which move with the body like the rest of the wall). Σ S = 0
by telescoping for any rigid motion and mesh; translating bodies are bit-identical.

The handout proposed the purely geometric form (all geometrically wet faces, the flux of removed
cells handed to their fluid neighbours ∝ θ^geo|S_f|). It is exact for a body turning in place,
but a cell next to a removed one keeps the flux through its closed face, O(|ω| R θ h), which its
continuity cannot carry, and an area-weighted share of the removed cell's flux does not cancel it
face by face. Measured on the Taylor–Couette case: cut-cell pressure errors at the rotating wall three to six
times larger,
and with the second-order treatment the velocity lost its second order (L2(U) 1.3e-4 → 8.9e-4 at
N = 160). An offline check with the exact TC field (face-centre fluxes, as the solver's
interpolation gives for a linear field) located the whole difference in the cells next to
removed ones; elsewhere the geometric and the open-face source are equally consistent. Two other
variants were tried and rejected: a flux correction θ_f (ω×(x'_f − x_f))·S_f in `ibm.flux()`
(worse), and a traction target whose normal component is taken from S (better than the
handout form, worse than the old source).

The price of the open-face form: a cylinder turning in place has S ≠ 0 in the cells next to a
removed cell (5 of 63 wall cells on an 81² mesh, |α^{n+1} − α^n| ≤ 0.018), against all 63 with
the point evaluation. A uniform field is still preserved exactly.

| test | before | after |
|---|---|---|
| cylinder turning in place, 81², Σ S | −1.9e-4 (one cell with d_wall ≈ 5 wrapped by the minimum image) | −1.2e-18 |
| cylinder turning in place, cells with S ≠ 0 | 63 of 63 | 5 of 63 |
| eccentric rotation, max over steps of \|Σ S\| | — | 8.1e-17 (4 ranks: same max \|S\| every step) |
| free stream, oscillating, Galilean cases | | byte-identical (movingBody, scalar, scalar_secondOrder) |

Taylor–Couette, L2 errors at N = 40/80/160:

| treatment | source | L2(U) | L2(p) |
|---|---|---|---|
| first order (movingBody) | point, x_P + d n_w | 1.11e-2, 7.61e-3, 3.83e-3 | 4.35e-2, 3.48e-2, 2.68e-2 |
| first order (movingBody) | open faces | 1.13e-2, 7.65e-3, 3.89e-3 | 3.28e-2, 3.12e-2, 2.21e-2 |
| first order (movingBody) | geometric + handover | 1.08e-2, 7.53e-3, 3.53e-3 | 1.22e-1, 5.22e-2, 1.48e-1 |
| second order | point, x_c + d_c n_w | 1.16e-3, 3.35e-4, 1.27e-4 | 2.17e-2, 1.70e-2, 1.30e-2 |
| second order | open faces | 1.17e-3, 3.04e-4, 8.09e-5 | 1.20e-2, 8.97e-3, 4.94e-3 |
| second order | geometric + handover | 2.82e-3, 6.32e-4, 8.93e-4 | 9.87e-2, 4.02e-2, 1.24e-1 |

The wall velocity of the no-slip traction is still a point value, x_P + d n_w or x_c + d_c n_w.
For a nearly full cell with a tiny wall segment d = αV/(2A_wall) can be several units, and the point
lay far outside the cell, wrapped by the minimum image (u_b = 4.5 instead of about 1 in one cell of
the rotating cylinder). The distance is now limited to where n_w leaves the cell's bounding box
(the traction coefficient keeps d): |u_b| = 1.01 in that cell. On movingBody the Taylor–Couette
errors become L2(U) = 1.13e-2, 7.30e-3, 3.76e-3 and L2(p) = 3.28e-2, 2.89e-2, 2.14e-2 (16–124
cells limited); translating cases are bit-identical.

### 2. Second-order pressure gradient

With `secondOrder true` the gradient is now

    ∇_V p|_P = Σ_f θ_f [p̄_f + g_f·(x'_f − m_f)] S_f + Σ_k [p_P + g_P·(x_k − x_c)] S_k,

m_f = w_f x_c,O + (1 − w_f) x_c,N the point the linear interpolate of centroid values belongs
to, applied on every open face of a cell with a centroid offset and on every partially wet face;
the wall pieces k are the wall polygon(s) of the geometric cut (the wall traces on the cut faces
chained into loops, fan-triangulated; the chord times the depth in 2-D) and the geometric wet
parts of faces closed by the post-processing. Σ_k S_k = Sw is checked (≤ 2e-11); only
Σ_k (x_k − x_c) S_k is stored. `forces()` uses the same pieces. Cells with too few fluid
neighbours for the least-squares fit (3-D slivers) take the mean gradient of their neighbours.

`ibmOrderTest` problems GL (p = G·x) and GQ (quadratic p), relative error of ∇_V p/(αV) in the
wall cells (L2), N = 20/40/80/160:

| problem | previous form | new form | new form, against V_geo |
|---|---|---|---|
| GL | 0.79, 0.76, 0.74, 0.98 | 3.4e-2, 2.5e-2, 1.1e-2, 9.1e-3 | 5.6e-15, 1.6e-14, 2.2e-14, 6.8e-14 |
| GQ | 0.29, 0.41, 0.34, 0.48 | 3.0e-2, 1.6e-2, 7.8e-3, 4.7e-3 | 2.2e-2, 1.2e-2, 6.0e-3, 3.7e-3 |

V_geo is the volume of the geometric polyhedron, obtained from the operator itself as
[∇_V x]_x = [∇_V y]_y. A linear pressure is exact to round-off, also in 3-D (sphere, both
sides, 12³–36³ cells, ≤ 1.2e-13); the error against αV is the sub-sampling error of α, large in
relative terms in the smallest cells. For a quadratic pressure the cut-cell gradient converges
at first order: the least-squares gradient of a one-sided stencil is first-order accurate, the
face values second order, and the cut cell divides them by its small volume. The handout
expected second order here; that holds for the volume-integrated gradient relative to h^(d−1),
not for the cell average. The Laplace/convection problems LD, LN, CD, CN are byte-identical.

Taylor–Couette with each fix alone and with both (second order, L2 at N = 40/80/160):

| | L2(U) | L2(p) |
|---|---|---|
| before | 1.16e-3, 3.35e-4, 1.27e-4 | 2.17e-2, 1.70e-2, 1.30e-2 |
| source only | 1.17e-3, 3.04e-4, 8.09e-5 | 1.20e-2, 8.97e-3, 4.94e-3 |
| gradient only | 1.15e-3, 3.60e-4, 1.43e-4 | 2.22e-2, 1.72e-2, 1.30e-2 |
| both | 1.17e-3, 3.04e-4, 8.21e-5 | 1.23e-2, 9.58e-3, 5.13e-3 |
| both, evaluation point limited to the cell | 1.17e-3, 3.04e-4, 8.12e-5 | 1.23e-2, 9.58e-3, 5.13e-3 |

With both, the velocity converges at order 1.9–2.0 in L2 and 1.8 in Linf (4.4e-4 at N = 160,
against 1.7e-3 before), and the pressure error halves. All of it comes from the source; the
linearly exact gradient changes the errors by at most 7 %. This answers the question of the
handout: the non-convergence of the cut-cell pressure (order 0.4 then 0.9, largest error 0.19 at
N = 160) is not a consistency defect of the gradient but the conditioning of the cut-cell
pressure rows (small-cell problem). Cell merging remains the remedy.

Hydrostatic balance (`tutorials/order/hydrostatic`): fluid at rest inside an immersed container
(an inverted cylinder of radius 1.1) around an immersed obstacle (radius 0.4, off the mesh
lines), g = (0, −1, 0) applied to the fluid volume αV (a coded fvOption), doubly periodic box
whose boundaries all lie in the solid, steady state (t = 10). Exact: U = 0, p = −y + const.
(A first version in a closed box with walls was dominated by the box walls: the same spurious
velocity without any body.)

| N | first order: max \|U\|, rms \|U\| | second order, before: max, rms | second order, after: max, rms | second order, cut-cell p error before / after |
|---|---|---|---|---|
| 40 | 4.0e-4, 5.8e-5 | 2.3e-4, 5.0e-5 | 1.1e-4, 2.8e-5 | 4.6e-3 / 6.9e-3 |
| 80 | 9.9e-5, 1.2e-5 | 6.0e-5, 9.1e-6 | 2.6e-5, 5.2e-6 | 2.6e-3 / 3.7e-3 |
| 160 | 2.5e-5, 2.5e-6 | 1.5e-5, 1.8e-6 | 6.4e-6, 9.3e-7 | 1.6e-3 / 2.2e-3 |

The linearly exact gradient halves the spurious velocity (second order in all cases); the
first-order treatment is unchanged. What remains is the imbalance between the pressure force of a
cut cell, now V_geo g, and the body force αV g: the sub-sampling error of α, relatively large in
small cells, which also makes the cut-cell pressure error slightly larger than with the previous
gradient. Weighting the body force with V_geo (available from the operator, [∇_V x]_x V), or
computing α from the polyhedron (exact in 2-D), would make the balance exact.

Drag of the periodic cylinder (`validation/richardson_drag.sh`, second order, both fixes; the
first-order runs are unchanged), per unit depth, N = 41/81/161/321 (the N = 321 and free-slip
N = 161 runs stopped at 8000 iterations with residuals ≤ 3e-8):

| wall | quantity | N = 41 | 81 | 161 | 321 | before the fixes (41 … 321) |
|---|---|---|---|---|---|---|
| no slip | total | 5.7847 | 5.8208 | 5.8324 | 5.8355 | 5.7834, 5.8205, 5.8322, 5.8354 |
| no slip | pressure | 2.7830 | 2.8461 | 2.8816 | 2.9025 | 2.8017, 2.8599, 2.8903, 2.9065 |
| no slip | viscous | 3.0018 | 2.9747 | 2.9507 | 2.9330 | 2.9817, 2.9606, 2.9418, 2.9288 |
| free slip | total | 2.8981 | 2.9728 | 3.0191 | 3.0523 | 2.9161, 2.9848, 3.0237, 3.0511 |

The gradient moves the components by up to 0.7 % and the no-slip total by less than 2e-4; the
orders are unchanged (no-slip total 1.6 then 1.9, pressure drag 0.8, free slip 0.7 then 0.5).
The forces, like the Taylor–Couette pressure, are limited by the conditioning of the cut-cell
pressure and by the O(α_min h) displacement of sliver absorption, not by the consistency of the
operators.
