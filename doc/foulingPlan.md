# Plan: reacting scalars and fouling growth on the cut-cell IBM (branch `fouling`)

Reference implementation: `foulingmodel_ibmdavide` (A. Krimmel, KIT GitLab, head
4c3e252, branches `main` / `intoduceEps`), an ibmSimpleFoam/ibmPimpleFoam built on the
static (`main`) cut-cell IBM with a temperature equation, three protein species, a
deposition law and a deposit field fed back to the IBM as `alpha.body`.

Base of the new branch: `scalar_secondOrder`, which contains `movingBody` and `scalar`
as ancestors. One branch with switches (`secondOrder true|false`, growth through the
moving-body path or quasi-static) covers the three requested variants; three parallel
fouling branches would triple the maintenance without a benefit. This is a decision the
user can overrule.

## 1. What the reference solver does

Fields: `U`, `p`, `T` (temperature), `folded`, `unfolded`, `aggregated` (protein
concentrations), `alpha` (deposit mass per cell volume, capped at ρ_dep·ε),
`foulingLayerIndicator` (1 once a cell is full: it joins the solid), `steelIndicator`
(pipe wall region), `alpha.body` = (1 − alpha/(ρ_dep ε))(1 − indicator), the fluid
fraction the IBM reads (`type alphaField`, iso-value `isoValue` in `ibmProperties`).

Equations (steady inner iterations, θ/α do not appear explicitly):

* Temperature, whole domain including deposit and steel (conjugate):
  `div(φ ρc_p, T) − laplacian(k, T) = 0`, with `k` and `ρc_p` blended per cell from
  fluid / deposit / steel values through `alpha.body` and `steelIndicator`.
* Species, fluid only: `div(φ, c) − laplacian(Γ, c) = reactions`, with the face
  diffusivity Γ set to zero on faces between cells of different `indicator` (so no
  diffusive flux into the deposit; the convective flux is blocked by φ = θ(u·S) = 0).
  Reactions (Arrhenius, `k exp(−E_a/(R T)) c^n`):
  unfolding `folded → unfolded` (order n_u, sink implicit in the folded equation),
  aggregation `unfolded → aggregated` (order n_a, sink implicit on `intoduceEps`).
* Deposition, explicit in the fouling time: in "surface cells" (cells with a face to a
  full cell, to a wall patch with alpha = 1, or across a processor patch)
  `r = (k_d − s_τ ū(z)) exp(−E_a,d/(R T)) unfolded^n_d · Σ_faces 1/h`,
  `alpha += r Δt`, capped; ū(z) is the cross-section mean velocity of the pipe slice
  at axial coordinate z from the inlet velocity and the free cross-section (a 1-D
  surrogate for wall shear). A cell reaching ρ_dep ε gets `indicator = 1`, its species
  are zeroed, and the time not used up (`leftoverDT`) is passed to a neighbour cell
  (`alphaLeftoverEqn.H`). The deposited protein is **not** removed from `unfolded`.
* Time decoupling: the fouling time is the OpenFOAM clock. Per fouling step the flow and
  T are re-converged (SIMPLE, up to `maxSubIter`, until the T residual is below
  `convergenceToleranceT`), then the species (until their residuals are below
  `convergenceToleranceSpecies`), then the deposit is advanced; a new
  `cutCellGeometry` is constructed from the updated `alpha.body`. Flow and species are
  re-solved only when a cell changed state (`indicatorUpdated`).

Points that do not carry over (design of the reference, not needed with the cut-cell
operators): the indicator/leftover logic, the surface-cell search (hex-only: 6-face
arrays, `(0 1 0)` face test, patch check on the first face only), the axial-slice ū(z)
with hard-coded z axis and inlet data, the `40000` area factor of the PIMPLE variant,
the zero-Γ interface faces, the pseudo-time `setTime/setDeltaT` manipulation, and the
per-step reconstruction of `cutCellGeometry` (re-registers the dictionary and fields
every step). The commit note "convection and diffusion not blocked at exactly the same
position" is resolved by construction here: both fluxes carry the same θ_f.

## 2. What the current branches already provide

| need | status on `scalar_secondOrder` |
|---|---|
| species transported in the fluid part only, no flux through the wall | `ibmScalarTransport`, zero-gradient wall: θ-weighted convection and diffusion, blanked solid (budget closes to 1e-11) |
| several scalars, per-scalar walls, fvOptions sources | `ibmScalarTransportList` (2026-09-11) |
| deposit as a body given by a fluid-fraction field, updated in time | `type alphaField; moving true;` re-fetched by `ibm.update()`, θ from the 0.5 iso-surface of the point-interpolated field |
| geometry-consistent growth (mass source, swept volume) | `movingBody` path: S = u_b·Sw, α^n from the wall flux |
| local wall area, wall distance, wall shear, wall temperature per cell | `Awall()`, `dWall()`, `Sw()`, `wallCoeff()`, `cellBody()`; viscous wall force per cell ν A_wall u_P/d_wall is already what `forces()` sums |
| mixed bodies (steel pipe as level set + deposit as alpha field) | union by min α, `cellBody` gives the owner |
| second-order operators | level-set bodies only: alpha-field bodies get no fluid centroid, wall distance stays the α V/(2 A_wall) estimate |
| constant flow rate in a periodic section | `ibmMeanVelocityForce` (α-weighted) |

## 3. Design

### 3.1 Reactions between scalars (`ibmReaction`, library)

Generic Arrhenius reaction between registered scalars of the list, read from
`scalarTransportProperties`:

```
reactions
{
    unfolding   { reactant folded;   product unfolded;   k0 1e3; Ea 5e4; order 1; }
    aggregation { reactant unfolded; product aggregated; k0 1e2; Ea 4e4; order 2; }
}
temperature T;          // scalar giving T in exp(-Ea/(R T)); absent: isothermal, exp term = 1
```

Rate per unit fluid volume `r = k0 exp(−Ea/(R T)) c_r^n`, entered as `Sp(−α r/c_r)` in
the reactant equation (implicit, as on `intoduceEps`) and `+α r` in the product
equation (explicit, evaluated with the latest reactant). Multiplying by α keeps the
solid free of reactions and makes the budget exact: production and consumption
`Σ α V r` go into the per-scalar budget line, so the sum over the chain closes. Stoichiometry
factor optional (`yield`). The scalars are solved sequentially in dictionary order and
iterated in the outer loop; nothing else changes in `ibmScalarTransport::solve()` beyond
one added source term. Effort: half a day plus the 0-D test.

### 3.2 Conjugate scalar mode (temperature)

Per scalar `solid conjugate;` replaces the wall treatment: no wall coefficient, no
blanking, no `theta` in the Laplacian. The equation becomes

```
ddt(C, T) + div(φ c_f, T) − laplacian(k, T) == fvOptions(T)
k = α k_f + Σ_bodies (1−α)_b k_b,   C = α C_f + Σ_bodies (1−α)_b C_b
```

with `k_b`, `C_b` per body (`solidDiffusivity`, `solidCapacity` in the body entry:
deposit and steel differ) using `cellBody` for the ownership of the solid part. φ is the
cut-cell flux, so convection stops at the wall on its own; conduction continues through
the deposit into the steel. This is the blended (grayscale) treatment of the reference
and is first order at the wall. The budget for such a scalar reports the boundary fluxes
instead of the wall flux. Effort: half a day. A second-order conjugate wall (two-sided
flux continuity at the wall segment) is a separate, later item.

### 3.3 Deposition as a surface reaction (`ibmDeposition`, library)

Operates on the wall cells of the IBM, which are exactly the cells that touch the
deposit or the pipe: no surface-cell search, no indicator, polyhedral and parallel for
free (all quantities are per cell).

```
// constant/foulingProperties
species        unfolded;          // deposited scalar
temperature    T;
depositDensity 1200;  packing 0.8;         // rho_dep, epsilon
rate  { k0 ...; Ea ...; order 1; shear { type linear; factor s_tau; }  }
body           deposit;           // alphaField body that receives the growth
growth         quasiStatic;       // or: moving (u_b = r/(rho eps) n_w through the moving-body path)
maxAlphaChange 0.2;               // adaptive fouling time step
consumeSpecies true;              // remove the deposited mass from the fluid (reference: false)
```

Per wall cell: `r_w = (k0 − s_τ τ_w) exp(−Ea/(R T_P)) c_P^n` [kg m⁻² s⁻¹] with the local
wall shear stress τ_w = ν |u_P,t| / d_wall (the same traction the no-slip term applies),
which replaces the 1-D ū(z) surrogate; the reference law can be reproduced by a
`shear { type meanVelocity; ... }` variant if wanted. Then

* species sink `Sp(−r_w A_wall /(c_P V))` in the deposited scalar (a Robin-type wall
  flux; switchable to reproduce the reference where nothing is removed);
* deposit mass `m += r_w A_wall Δt_f` per cell, solid fraction `s = m/(ρ_dep ε V)`,
  body fluid fraction `alpha.body = 1 − s`, written and registered so the alphaField body
  picks it up at the next `ibm.update()`;
* overflow: if `s` would exceed 1 the excess mass goes to the face-neighbour wall cells
  with room (mass conservative), and `Δt_f` is chosen so that `max Δs ≤ maxAlphaChange`
  (adaptive fouling step) so overflow stays rare;
* budget: `Σ m V` against `∫ Σ r_w A_wall dt` and the species budget, so that
  inflow − outflow − reacted − deposited − Δcontent closes.

Cells below `alphaMin` are absorbed by the geometry as usual: a filled cell becomes
solid, its neighbours become wall cells, growth continues. Effort: one to two days.

### 3.4 Solver `ibmFoulingFoam` (quasi-steady outer loop)

New application rather than a modified ibmSimpleFoam. The OpenFOAM clock is the fouling
time; per step:

1. `ibm.update()` if the deposit changed by more than `geometryUpdateTolerance`
   (alphaField body, `moving true`, zero velocity: geometry re-fetched, S = 0, the
   quasi-static case; `growth moving` supplies u_b instead);
2. inner SIMPLE iterations for U, p and the conjugate scalars until their
   `residualControl` is met (the fluid and T are steady on the fouling time scale);
3. sequential species iterations with the reactions until their residuals are met
   (`steadyState` ddt, relaxation; no pseudo-time);
4. deposition step, adaptive Δt_f, write.

The transient alternative (ibmPimpleFoam + scalars each step + deposition each step)
needs no new code beyond 3.1 to 3.3, but its time step follows the flow, which is what
the time decoupling avoids; it stays as a check case. Effort: one day.

### 3.5 Second order (`secondOrder true`) with the deposit body

The operators of `scalar_secondOrder` (centroid-consistent fluxes, cutDeltaCoeffs) need
the fluid centroid of cut cells, which the alpha-field route does not compute. Add it from
the reconstructed geometry: wet-face centroids are already available from the face
polygons; the cell centroid follows from the divergence theorem over the wet faces and the
wall segment (whose centroid is the mean of the cut-edge points). The centroid wall
distance uses the signed distance to the reconstructed segment. Effort: one day. Until
then `secondOrder true` runs but degrades to first order in the deposit cells.

## 4. Missing features, ranked

1. **Accurate θ from an α field** (README "possible improvements"): the 0.5 iso-surface of
   the point-interpolated field gives θ errors up to 0.1 and a wall position error of
   about h/10. With the deposit defined only by α this is the accuracy floor of the whole
   model. Fix: per cell, the iso-value that reproduces the cell's α (as geometricVoF's
   isoAlpha), faces cut at that value, shared faces averaged. One day; do first.
2. Reactions between scalars (3.1).
3. Conjugate scalar mode with per-body solid properties (3.2).
4. Surface reaction / deposition class with local shear and overflow handling (3.3).
5. Quasi-steady solver with geometry update on demand (3.4).
6. Pipe wall as a domain patch: wall patches are not IBM walls, so A_wall = 0 there and no
   deposition starts on them. Recommended: model the pipe as an IBM body (plane or
   inverted cylinder; the reference already treats the steel as a region), which also gives
   the conjugate conduction into the steel. Alternative, if the mesh must keep its wall
   patches: count wall-patch faces as wall area in `ibmDeposition` (half a day).
7. Fluid centroids for alpha-field bodies (3.5) so `secondOrder true` is second order in
   the deposit cells too.
8. Growth through the moving-body path (`growth moving`): u_b = r_w/(ρ_dep ε) n_w in the
   wall cells provides S and the swept α consistent with the deposit; only relevant if
   the transient solver is used, since S ~ 1e-9 m³/s is negligible in the quasi-steady
   loop.
9. Second-order conjugate wall (flux continuity across the segment): later.
10. Not planned: turbulence (laminar only, as before).

## 5. Validation

* Reactions: 0-D batch (uniform fields, no flow, no body) against the analytic solution
  of the two-step chain; 1-D plug flow against the analytic profile.
* Conjugate T: 1-D three-layer slab (still fluid, deposit, steel) against the piecewise
  linear profile; annulus conduction around the cylinder, order of convergence.
* Deposition: channel with a uniform `unfolded` and constant rate: the deposit thickness
  grows as r t/(ρ_dep ε); check the wall distance and the fluid volume against it, and
  that the flow rate held by `ibmMeanVelocityForce` raises the wall shear as the
  section narrows. Then the shear-limited case (k0 − s_τ τ_w → 0) reaches a steady
  thickness.
* Budgets: species chain and deposit mass close to solver tolerance over a full fouling
  run; the wall-flux, reaction and deposition lines are printed per scalar as now.
* Comparison with the reference: needs a case from the fouling repository (none is
  included); the deposit thickness history and the outlet concentrations are the
  quantities to compare, with `consumeSpecies false` and the mean-velocity shear model
  to match its assumptions.

## 6. Inputs needed

* A representative case of the reference solver (mesh, `transportProperties` with the
  constants, initial fields), so the tutorial reproduces a run the user knows.
* Units and orders of the three reactions and of the deposition law (the reference uses
  dimensionless "dimCorrector" factors); ε is the packing fraction of the deposit?
* Whether the deposited protein should be removed from the fluid (`consumeSpecies`).
  Physically yes; the reference does not.
* Local wall shear (recommended) or the cross-section mean velocity model for the
  shear-dependent rate; both can be offered.
* Pipe wall as an IBM body (recommended) or as a mesh patch (item 6).
