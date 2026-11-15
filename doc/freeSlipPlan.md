# Handout: zero-shear free slip on curved walls (written 2026-11-03; implemented 2026-11-14)

## Status (2026-11-14): implemented behind `zeroShear` (default false)

Commits f4f2366 (switch), 361c23c (mergeH, finding 12 of `doc/validationSuite.md`), and the
normal viscous flux (commit after 361c23c). Runs and scripts: `~/runs/zeroshear` (`scripts/`:
`v2zs.py`, `tc_zs.py`, `sph_zs.py`, `mms_zs.py`). G1–G9 bit-identical with the switch off.

What the switch does (`cutCellGeometry::calcZeroShear`, `zeroShearSource`; user decisions
2026-11-14 in brackets):
- [global key `zeroShear` in ibmProperties] every freeSlip wall cell of a body with a smooth
  level set gets the wall flux of the Laplacian form ν (A_w/V)(H_w·w_t,w + ω × n + (∂u_n/∂n) n),
  w = u − u_b, n into the fluid (fixes 1, 2, 5, 6, 8, 10).
- [curvature from the level set] H_w: finite-difference Hessian of the body's level set at the
  fluid centroid (step 1e-3 of the cell size), projected on the tangent plane, moved to the wall
  per principal direction κ_w = κ/(1 − φκ); limited where the radius of curvature is under twice
  the distance (the friction then stays below the no-slip coefficient). triSurfaceMesh and
  alpha-field bodies keep ∂u/∂n = 0 with a warning (fix 7; normal fit for them: open).
- [Robin factor] tangential wall value from the value point, w_t,w = (I + d H_w)⁻¹ w_t,P (exact
  for a rigid rotation), plus the tangential offset of the wall-piece centroid u_P + δ·∇u at
  second order (`zeroShearTangential`, default true; changes V2 by < 1e-3) (fix 3).
- normal part ∂u_n/∂n = −tr(P ∇u P) from continuity and the least-squares gradient, explicit.
  n·∇u·n, equal for a solenoidal field, depends on the cell's own normal velocity with a
  negative self-coefficient and made SIMPLE without linking and PIMPLE with it diverge.
- implicit part: the largest |eigenvalue| of the friction tensor (concave walls: negative
  friction, +|f| implicit, −2|f| explicit).
- forces: the flux in the viscous force (force = driving force × fluid volume in V2); moment:
  the transpose part ν Σ S_w × u_w replaces the rigid-rotation term for these cells (fix 4).

Results (zeroShear true, second order):
- Exact rigid rotations (Taylor–Couette geometry, ibmPimpleFoam to steady state; 3-D spherical
  shell R = 0.5 in 1, ibmSimpleFoam), u L2 at N = 40, 80, 160 (d/h = 8, 16, 32), torques → 0 at
  order 2: concave free-slip outer wall, inner turning (A) 9.1e-4, 1.5e-4, 2.7e-5; stress-free
  inner cylinder turning, fluid at rest (Bi) 9.3e-5, 2.6e-5, 8.8e-6; static stress-free inner,
  outer turning (C) 2.3e-4, 6.5e-5, 1.2e-5; sphere turning in fluid at rest (SB) 4.0e-4,
  8.7e-5, 2.3e-5; static sphere, outer shell turning (SC) 1.6e-3, 6.1e-4, 2.1e-4 (order 1.4–1.5).
  Old free slip: A 0.81, C 0.15, SC 0.12, not converging.
- Manufactured solution (`mms_zs.py`: TC geometry, inner wall stress-free, outer no slip,
  ψ = (r − R1)(r − R2)²(1 − 6r) sin 2θ, slip 2 sin 2θ and ∂u_r/∂r = −8 cos 2θ at the wall,
  p = 0): u L2 0.0222, 0.0055, 0.0014, 3.4e-4 (N = 40 … 320, order 2.0), p L2 0.10, 0.050, 0.017. Without the
  normal flux u L2 0.051, 0.018, 0.0092, 0.0040 (order ≤ 1.2) and the wall-cell pressure an O(1)
  layer that grows (L∞ 5.1 … 6.9): the "normal viscous stress carried by the cut-cell pressure"
  of the diagnosis is not consistent, it costs one order. No-slip control: 3.5e-3, 9.0e-4,
  2.3e-4.
- V2 (periodic cylinder Re_D = 1, drag from the driving force, = the reported force), linked:
  3.9641, 3.9895, 3.9974, 3.9990 at N = 41, 81, 161, 321 (orders 1.7, 2.4; limit ≈ 3.9993);
  body-fitted zero shear 4.0234, 4.0132, 4.0068, 4.0032 at M = 32 … 256 (order ≈ 0.8, limit
  ≈ 3.998–4.000). Without the normal flux the series converged at order ≈ 0.9 (3.8933 …
  3.9932 at N = 641); without the Robin factor (the experiment) more regular in V2 (3.9622 …
  3.9980 at N = 641) but first order in the exact tests (error cancellation).
- Unlinked (alphaMin 0.02): 3.9747, 3.9936, 4.0032, 4.0119, above the reference and rising:
  the absorbed slivers (α < 0.02, made solid) leave steps of closed faces across which the slip
  velocity cannot pass (finding 13). Same with the old free slip; the concave Taylor–Couette
  case A unlinked: u L2 0.082 (N = 40) with alphaMin 0.02, 4.3e-4 with alphaMin 1e-9; V2
  unlinked with alphaMin 1e-9: 3.9748, 3.9936, 3.9980, 3.9991, as the linked series.

Open:
- Default on: changes every free-slip result (G7 free-slip tutorials, V2); baselines need the
  user. Then fix 9 (free-slip pimple tutorial, notebook comparison, `cutCellIBM.tex`).
- Unlinked free slip with alphaMin > 0 (finding 13).
- triSurfaceMesh and alpha-field bodies (curvature from the wall normals).
- V11 BUB (3-D, needs the V10 SPH setup); MIG-G2 stress-free sphere (V6).
- Turbulent free slip: the force uses the solver's wall viscosity (molecular for the
  turbulent solvers), the momentum nuEff.
- First order (secondOrder false): no normal flux and no tangential correction.

Branch `feature_iddes_from_scalar_secondOrder`. Diagnosis and data:
`doc/secondOrderStudy.md`, "Free slip on curved walls: the missing curvature shear";
`validation/data/secondOrder/freeSlip_zeroShear.txt`. Postponed in favour of the first-order
Taylor–Couette problem with linking.

## Diagnosis

- The free-slip wall has no flux in the Laplacian form ∇·(ν∇u), so it imposes ∂u/∂n = 0 for
  every Cartesian component. Zero tangential stress on a wall with curvature tensor H (Hessian of
  the signed distance, positive in the fluid) is P ∂u/∂n = H·u_t (cylinder: ∂u_θ/∂r = u_θ/R).
  The IBM wall therefore pushes the fluid along with a traction ν H·u_t: a spurious thrust.
  OpenFOAM's `slip` has the same defect.
- Periodic cylinder, Re_D = 1: body-fitted slip 3.058, zero shear 4.007 (M = 128); the old IBM
  free slip converges to the slip value. The curvature term is about 24 % of the drag.
- It is also the source of the 1–3 % linked/unlinked free-slip spread: with the zero-shear term
  the spread is 0.06–0.3 % at N = 321 and all variants converge to about 4.00.
- The normal part of the same inconsistency (∂u_n/∂n = 0 imposed, continuity needs
  −∇_t·u_t) stays: its normal viscous stress is carried by the cut-cell pressure (one-cell layer,
  first order in the pressure split, consistent in the total force).

## The experiment (reference implementation)

Worktree `runs/ibmSA_study/freeslip_fix` (detached at 8de56d7, uncommitted), files
`applications/solvers/ibmSimpleFoam/createFields.H` (end) and `UEqn.H`: switch `zeroShearSlip`
in `ibmProperties`; κ from the Laplacian of the level set at the cell centre by central
differences (ε = 1e-4), moved to the wall with κ_w = κ/(1 − φκ); coefficient κ_w A_w/V on
free-slip wall cells; `fvm::Sp(ν c, U) − ν c (U·n_w) n_w`. Runs: `runs/ibmSA_study/freeslip_ibm`
(`run.sh`, `table.py`), body-fitted `runs/ibmSA_study/freeslip_bf` (`validation/cylinder_bodyFitted.py`).

## Fixes to do

1. **In `cutCellGeometry`, not the solver.** A `zeroShearCoeff_` (or a tensor coefficient) built
   with the wall coefficients (`cutCellGeometry.C`, wall-coefficient loop after "Dirichlet wall
   coefficient"), used by ibmSimpleFoam, ibmPimpleFoam and the turbulent solvers.
2. **Curvature tensor in 3-D.** H = Hessian of the level set at the wall point (shape operator in
   the tangent plane), friction ν A_w H·u_t/V. Implicit part with the largest positive eigenvalue
   (or the trace/2), the rest explicit. Exact for the analytic bodies; triSurface bodies need a
   smooth enough level set (check the Hessian noise, maybe fit over the wall pieces).
3. **Wall point and second order.** Evaluate H at the wall point, not the cell centre; with
   `secondOrder`, use the wall value u_t,w = u_t,P + P(∇u·(x_w − x_v)) (x_v: value point,
   `valueOffset_`), as for the no-slip traction.
4. **Forces.** Add ν H·u_t A_w to the viscous force of free-slip bodies (`calcForces`, next to the
   noSlipCoeff term); check force = driving force × fluid volume in the periodic cylinder.
5. **Concave walls** (κ < 0, e.g. the outer cylinder of Taylor–Couette with free slip): negative
   friction, explicit only.
6. **Moving bodies.** u·n = u_b·n on the wall adds the tangential gradient of u_b·n; for a rigid
   body P ∂u/∂n = H·(u − u_b)_t + (rigid-rotation part). Derive and test with a translating and a
   rotating free-slip cylinder.
7. **Alpha-field bodies** have no level set: no curvature (keep ∂u_t/∂n = 0, warn).
8. **Linking.** Per-cell coefficients are fine (the group equation sums them); check the slaves
   use their own A_w, not the group distance.
9. **Default.** Behind a switch first; it changes every free-slip result by about 30 %. Then
   rerun the free-slip Richardson series, the free-slip pimple tutorial and the notebook comparison
   (`validation/fig_pimple_freeSlip_vs_notebook.png`, the notebook has the same ∂u_t/∂n = 0
   condition), and update `cutCellIBM.tex` (free-slip paragraph claims zero tangential traction).
10. **Optional:** the normal viscous stress 2ν∂u_n/∂n at the wall (full traction form), which
    would remove the one-cell pressure layer of free-slip walls.
