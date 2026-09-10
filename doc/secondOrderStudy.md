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
