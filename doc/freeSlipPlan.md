# Handout: zero-shear free slip on curved walls (written 2026-11-03, not started)

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
