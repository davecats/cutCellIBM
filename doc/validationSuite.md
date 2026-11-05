# Validation suite: design

Started 2026-11-05. Principles: `doc/developerNotes.md`. Each case states which components it
**stresses** and which it **relies on**. A **gate** is a cheap check with a binary tolerance that
runs on every change; it must pass before the cases that rely on it are run or interpreted. A
**validation case** compares with an independent reference (exact solution, body-fitted solution
of the same problem, literature) with a stated tolerance and an error bar from mesh refinement.
Every check reports one of PASS, FAIL, or REPORT (no tolerance yet, or a known open defect whose
fix the case will accept).

## 1. Components

| id | component | where |
|---|---|---|
| K1 | geometry: cut fractions, wall area and normal, sub-sampling, closure Σ S_k = S_w; 2-D, 3-D; analytic, STL, alpha-field bodies; invert; several bodies | `cutCellGeometry.C`, `ibmBody` |
| K2 | linking: groups, shared value point, mergeAlpha | `cutCellGeometry.C` |
| K3 | no-slip wall traction ν(u − u_b)/d, second-order wall gradient | wall coefficients |
| K4 | free-slip wall (zero shear still missing, `doc/freeSlipPlan.md`) | wall coefficients |
| K5 | pressure gradient in cut cells (linearly exact form) | `grad()` |
| K6 | continuity: cut-face fluxes, second-order continuity flux, p–U coupling (SIMPLEC, PIMPLE) | solvers, `cutCellGeometry` |
| K7 | convection at the wall (`ibmConvection`) | `cutCellGeometryTemplates.C` |
| K8 | wall force (and torque, C1) integration | `forces()` |
| K9 | moving bodies: wall-flux source, swept, absorbed and uncovered cells, rotation, minimum image | `cutCellGeometryMoving.C`, `bodyMotion` |
| K10 | flow-rate forcing | `ibmMeanVelocityForce` |
| K11 | passive scalars: fixedValue and zeroFlux walls, budget | `ibmScalarTransport` |
| K12 | turbulence near IBM walls: SA with the IBM wall distance | `ibmSpalartAllmaras`, `ibmWallDist` |
| K13 | wall law: SA law, sample at wallSampleDistance, blending in h+ | `cutCellGeometryWallFunction.C` |
| K14 | parallel: decomposition independence, groups and samples across processors | everywhere |
| K15 | AMR (after the merge of `feature_amr`): refinement, unrefinement, mapping, load balancing | `ibmRefineFvMesh` |

## 2. Inventory of what exists

| existing | stresses | reference today | becomes |
|---|---|---|---|
| `tutorials/order/laplaceDisk` (`ibmOrderTest` LD LN CD CN GL GQ, N = 20…320, + 3-D sphere GL) | K1 K3 K5 (K2) | manufactured solutions; `data/secondOrder/orderTest*.txt` | gate G2 |
| `tutorials/order/hydrostatic` | K5 K6 | rest state; `data/fixes/hydrostatic.txt` | gate G3 |
| `tutorials/order/taylorCouette` | K3 K5 K6 K2 (K8 torque) | exact; `taylorCouette_*.py` | validation V1 |
| `cylinderRe1` (+ `_stl`, `_alphaField`, `_parallel`), `richardson_drag.sh` | K1 K3 K5 K8 K14 | notebook (`reference_notebook.py`), Richardson | gate G7 (regression), validation V2 |
| `cylinderRe1_pimple`, `_pimple_freeSlip`, `twoCylinders` | K4 K6, several bodies | notebook | G7, V2 |
| `freeStreamTranslating` | K9 | exact (uniform flow preserved) | gate G4 |
| `rotatingCylinder`, `eccentricRotation` | K9 (source) | Σ S = 0, max \|S\| | gate G5 |
| `galileanTranslating`, `oscillatingCylinder` | K9 | notebook (`reference_moving.py`), Galilean equivalence with cylinderRe1 | validation V3 |
| `tutorials/scalar/*` | K11 | budget (`scalar_budget.py`) | gate G6, V4 |
| `tutorials/turbulence/channelSA` | K12 K13 K10 | body-fitted SA channel (same model) | validation V7 |
| `tutorials/turbulence/naca0012` | all of K1–K13 together | body-fitted TMR grid (same model); `data/naca0012` | validation V8 (integral, last) |
| `runs/ibmSA_study/freeslip_bf` (running) | — | body-fitted zero-shear periodic cylinder | reference of V2 free slip |

Gaps: no laminar case at moderate Re (convection at the wall, K7, is tested only inside the
turbulent NACA, where it cannot be separated from K12/K13); no 3-D flow case; no torque output; no
decomposition gate beyond `cylinderRe1_parallel`; regressions are checked by hand
(`data/secondOrder/tutorials_secondOrderDefault.txt`).

## 3. Gates (minutes each, all on every change)

| id | check | tolerance |
|---|---|---|
| G1 | geometry: Σ αV against the exact area/volume (cylinder 2-D, sphere 3-D, inverted), closure Σ S_k = S_w, linear pressure exact | volume error converging at ≥ 1.8; closure and linear p ≤ 10⁻¹⁰ |
| G2 | operator orders, `ibmOrderTest` LD LN CD CN GL GQ, N = 20…160, mergeAlpha 0 and 0.5 | order ≥ 1.8 (LD LN CD CN GL), ≥ 0.9 (GQ); errors within 1 % of the stored baseline |
| G3 | hydrostatic container | max \|U\| ≤ stored baseline × 1.05 and converging at ≥ 1.8 |
| G4 | free stream around a translating body | max \|U − u_b\| ≤ 10⁻¹⁰ |
| G5 | rotation sources (`rotatingCylinder`, `eccentricRotation`) | \|Σ S\| ≤ 10⁻¹² ; max \|S\| ≤ stored × 1.05 |
| G6 | scalar budgets (`tutorials/scalar`) | closure ≤ 10⁻¹⁰ relative |
| G7 | regression: forces of every tutorial against the stored baseline | ≤ 10⁻⁶ relative (larger = intended change, baseline updated in the same commit with a reason) |
| G8 | decomposition: cylinderRe1 and channelSA (linked, wall law) serial vs 4 and 6 ranks, converged | ≤ 10⁻⁸ relative in the forces |

## 4. Validation cases

| id | case | stresses | relies on | reference | tier |
|---|---|---|---|---|---|
| V1 | Taylor–Couette, N = 40…320 | K3 K5 K6 K2, torque | G1 G2 | exact | 1 |
| V2 | periodic cylinder Re_D = 1, no slip and free slip, Richardson | K3 K4 K8 | G1–G3 | body-fitted Richardson (no slip: notebook; free slip: `freeslip_bf`) | 1 |
| V3 | translating and oscillating cylinder | K9 | G4 G5 | Galilean equivalence with V2, notebook | 1 |
| V4 | scalars, static and moving | K11 | G6 | body-fitted steady cylinder scalar | 1 |
| V5 | **new:** laminar cylinder in a channel, Re_D = 40 steady and 100 periodic, three meshes | K7 K3 K5 K8 | V2 | body-fitted OpenFOAM in the same domain (separation length, C_d, C_p, St); literature only as plausibility | 1 |
| V6 | **new:** MIG, sphere in a pipe (below) | 3-D K1 K3 K4 K7 K8 K10 | MIG-G1, MIG-G2, C1 C2 | Rivero-Rodriguez & Scheid 2018 | 2 |
| V7 | SA channel: cut fraction f sweep, tilt tan α, wall-resolved and wall-function, linked and unlinked | K12 K13 | V1 V2 V5 | body-fitted SA channel | 1 (wall-resolved), 2 (sweeps) |
| V8 | NACA 0012, α = 0, 10, wall-resolved and wall-function | integral | everything | body-fitted TMR grid | 2 |
| V9 | AMR (after the merge): V5 and V3 with refinement against the uniform fine mesh | K15 | V3 V5 | the same solver on the uniform mesh | 2 |
| V10 | **new:** SPH, rigid sphere in uniform flow, Re 20–300 (§10) | 3-D K1 K3 K7 K8, separation | G1–G3, V2 | Johnson & Patel 1999; Magnaudet et al. 1995 | 2 |
| V11 | **new:** BUB, stress-free sphere in uniform flow, Re 50–300 (§10) | K4 in 3-D | V10 | Blanco & Magnaudet 1995; Moore; exact surface identity | 2 |
| V12 | **new:** ROT, cylinder rotating in a uniform stream, Re 200, α 0–4 (§10) | K3 with wall velocity, K8 torque, K7 | V2 V5, C1 (C2 for α ≥ 2 steady) | Mittal & Kumar 2003 | 1 (α = 0, 1, 2), 2 |
| V13 | **new:** OSC, cylinder oscillating in-line in fluid at rest, Re 100, KC 5 (§10) | K9, K8 (added mass) | V3 | Dütsch et al. 1998 | 1 |

Tier 1: hours in total on hydra, the baseline before the AMR merge. Tier 2: days.

## 5. Options

Every case runs with the defaults (secondOrder, ibmConvection, mergeAlpha 0.5, wall law SA,
wallSampleDistance 1.5, ibmConvectionWallLaw off). Off-default, one at a time, only where the
option acts: mergeAlpha 0 (G2, V1, V2, V5, V7), ibmConvection off (V5, V7, V8), wall law none
(V7, V8), ibmConvectionWallLaw on (V7, V8). Full product only for interacting options:
ibmConvection × mergeAlpha in V5 and V7 (the NACA study showed they interact). secondOrder false
only in G2 and V1 until first order is removed. nSubSamples (8, 32) only in G1.

## 6. Driver and report

- `validation/suite/`: one module per gate or case (setup from the tutorial or a generator,
  commands, extraction to a dict of metrics, reference, tolerance, components stressed and relied
  on), and `suite.py` (`run <ids> --tier --host hydra`, `report`). Case directories outside the
  repository (`runs/suite/<commit>/<id>/`); per-case results as small JSON in
  `validation/suite/results/<commit>/`, the accepted baseline in `validation/suite/baseline/`.
- Report: text summary (one line per check with value, reference, tolerance, status), a
  component × case matrix coloured by status, plots per case (convergence plots with the
  observed order; profiles and curves against the reference), as markdown + PNG in the results
  directory and one HTML page published as a private artifact.

## 7. Order of work

1. Driver and gates G1–G8, tier 1 cases (V1–V5, V7 wall-resolved, V12 α = 0, 1, V13) on the
   current branch: the baseline.
2. Merge `feature_amr`; accepted when the gates pass and tier 1 agrees with the baseline.
3. C1, C2 and the free-slip zero-shear fix, each accepted by its gates; then tier 2 (V10 SPH and
   V11 BUB before V6 MIG, V7 sweeps, V8, V9, V12 complete).

## 8. Case MIG (V6): inertial migration of a sphere in a pipe (Rivero-Rodriguez & Scheid 2018)

Reference: J. Rivero-Rodriguez & B. Scheid, *Bubble dynamics in microchannels: inertial and
capillary migration forces*, J. Fluid Mech. 842 (2018) 215–247, doi:10.1017/jfm.2018.78 (steady
FEM, COMSOL, mesh-converged; the PDF is not in the repository).

### Problem

A periodic train of spheres of diameter d in a circular pipe of diameter d_h = 1, mean flow
J = 1, centre at transverse position ε (eps_* = (1 − d)/2 is contact), in the frame of the sphere:
pipe wall moving at −V e_x, mean flow J − V relative to the sphere, period L. Re = J d_h/ν, Ca = 0
(undeformed sphere). Two interfaces:

- **rigid** (infinite Marangoni, = solid sphere): no slip, the sphere rotates at Ω e_z,
  torque-free;
- **stress-free** (clean bubble, = inviscid drop): zero tangential stress, no penetration.

V (and Ω) follow from the force (torque) balance along x (about z). Output: the balanced transverse
force per unit sphere volume, f = L_y/V_B (f/Re in units ρ = J = d_h = 1 is L_y/V_B), V,
Ω, and the pressure correction factor β, Δ_s p_T/∂_x p_P = 1 + β α_G, α_G = V_B/(L π/4).

### Reference data (in the repository)

- `validation/reference_riveroScheid2018.py`: the polynomial fits of appendix C, pure linear
  inertial regime (Re → 0): f1 = f/Re, V0/V_P, β0, Ω0/8ε against (d, ε), 0.1 ≤ d ≤ 0.9,
  ε ≤ 0.95 eps_*, < 1 % error (paper). The captions of tables 3 and 4 are swapped (table 3 is the
  rigid interface); `check()` reproduces the text: centred position unstable for d < 0.727 rigid
  (text 0.73) and d < 0.829 stress-free (text ≈ 0.85), d = 0.4 rigid maximum f/Re 0.297 at
  0.44 eps_* and neutral position 0.745 eps_* (figure 4a).
- `validation/data/riveroScheid2018/` (from the vector graphics of the PDF,
  `validation/riveroScheid2018_extract.py`, ~0.3 % of the axis range):
  `fig9_rigid_d0.4.dat`, f/Re against ε/eps_*, d = 0.4 rigid, Re = 8, 16, 32, 64, 128 and the
  linear curve; `fig10b_rigid_neutral.dat`, `fig10d_stressFree_neutral.dat`, neutral equilibrium
  position against d for Re = 4…128.
- Cross-checks: the linear curve of figure 9 agrees with the fit to 0.007 up to 0.8 eps_*, but at
  0.95 eps_* the fit gives −0.81 against −0.96 (the fit is poor at the wall: use the figure
  there); figure 10b at Re = 4 agrees with the neutral positions of the fit to three digits
  (d = 0.2, 0.4, 0.6: 0.272, 0.224, 0.158).

Scope (decided 2026-11-05): the zero-capillary part only (§3.1, figures 4–10, tables 3–4).
Everything with Ca > 0 (§3.2–3.3, figures 11–19, table 5, eq. C2) needs a deformable interface
with surface tension. The setup below (IBM pipe, Re = 2 for the linear regime, gates, cases,
tolerances) was accepted the same day.

### Components

Stressed: three-dimensional cut cells on a doubly curved convex wall (first 3-D flow case; 3-D
so far only in the operator tests); no-slip traction with the wall velocity of a rotating body
(rigid); free slip on a doubly curved wall (stress-free: the zero-shear term of
`doc/freeSlipPlan.md`, fix 2, is missing, so a failure is expected until it is implemented, and
this is its 3-D acceptance case); pressure and viscous force integration in 3-D, where the lift is
small against the stresses (fore–aft cancellation); torque; convection at the body at Re = 2–128
(`ibmConvection`, `secondOrder`); linking in 3-D; `ibmMeanVelocityForce` with a body in the
domain.
Relies on: the pipe wall (gate MIG-G1), periodicity, symmetry plane.
Not involved: turbulence, wall functions, moving geometry, scalars.

### Setup

- Mesh: Cartesian, x periodic (L = 3, checked with L = 2, 4 at the worst case Re = 128), half
  domain z ≥ 0 with a symmetry plane through the sphere centre (as in the paper), mirror-symmetric
  in x about the sphere centre (so the creeping-flow lift vanishes discretely too). Base
  h = 1/32 (refined near the pipe wall if G1 asks for it), static refinement (refineMesh, hanging
  nodes) to h = d/32, d/48, d/64 in a box around the sphere containing the sphere and its gap to
  the wall. The gap must hold ≥ 3 fine cells: cells cut by both bodies attribute their force to
  one body only, so ε ≤ 0.9 eps_* at d = 0.4 on the finest mesh.
- Pipe: IBM cylinder with `invert true`, sliding at −V along its axis. Alternative considered:
  body-fitted O-grid pipe (no pipe error, wall velocity a plain patch value) — rejected as the
  primary because the sphere then sits in skewed cells near the wall, and IBM on non-Cartesian
  cells is untested; kept as a cross-check at one point if wanted.
- Sphere: `sphere` body, rigid: `noSlip`, rotation Ω e_z about its centre; stress-free:
  `freeSlip`.
- Flow rate: `ibmMeanVelocityForce` with Ubar = (J − V) A_f L/(A_f L − V_B), the fluid-volume
  average that corresponds to a flux A_f (J − V) through every cross-section (A_f: discrete pipe
  section from ibmAlpha in a slice away from the sphere). Δ_s p_T from its pressure gradient G;
  Δp = (G − G_0) L with G_0 of the empty pipe on the same mesh, so the pipe discretisation error
  cancels in β.
- ν = 1/Re, ρ = J = d_h = 1, steady SIMPLEC, ibmSimpleFoam.
- Force balance: x-force on the sphere = periodic-pressure + viscous force + G V_B (the mean
  pressure gradient acting on the sphere volume). Outer Newton iteration on (V, Ω) (rigid) or V
  (stress-free) until |F_x|/V_B and |T_z|/V_B are below 10⁻³ |f|_max of the configuration; the
  Jacobian from two perturbed runs at the first point, then Broyden updates, restarts from the
  previous solution. Linearity in Stokes flow makes this two to three outer steps.
- Linear regime: Re = 2 (f/Re = f1 + O(Re²); figure 9 shows < 0.003 between Re = 8 and the
  linear curve at mid ε, so Re = 2 is within 2·10⁻⁴); one check at Re = 1.

### Prerequisite code changes (general, not case-specific)

- **C1 torque output:** torque of pressure and viscous wall forces about the centre of rotation
  of each body, in `forces.dat` and the log.
- **C2 tangential wall velocity on a static geometry:** a body whose motion leaves it invariant
  (sphere rotating about a diameter, body of revolution about its axis, cylinder sliding along its
  axis) keeps its geometry fixed, has no wall flux, and only supplies a wall velocity to the
  no-slip traction. Today such a body counts as moving and ibmSimpleFoam refuses it. The
  invariance is checked at start-up (level set at the cut points after a virtual small motion),
  fatal otherwise. Used here for the rotating sphere and the sliding pipe.
- **C3 driver:** `validation/suite/mig/` setup, outer (V, Ω) iteration, extraction, report.

### Gates

- **MIG-G1 pipe:** empty pipe, IBM wall sliding at −V: exact u = −V + 2J(1 − 4r²),
  ∂_x p = 32/Re. Error of G and of the profile on the three base meshes, order ≥ 1.5, G within
  0.2 % on the base mesh used (otherwise refine at the pipe wall). Stresses the concave 3-D IBM
  wall and C2 (sliding).
- **MIG-G2 sphere, Stokes:** sphere in a large periodic box (Re = 0.01): drag of the translating
  rigid sphere 6πμRU, of the stress-free sphere 4πμRU (Hadamard–Rybczynski), torque of the
  rotating sphere 8πμR³Ω, two box sizes with extrapolation in R/L_box; d/h = 16, 32, 48. Rigid and
  torque within 1 % at d/h = 32 and converging; stress-free reported (expected to fail before the
  free-slip fix). Stresses force and torque integration (C1) and C2 (rotation) without the pipe.

### Cases

| id | interface | d | Re | ε/eps_* | reference | quantities |
|---|---|---|---|---|---|---|
| MIG-R1 | rigid | 0.4 | 2 (1) | 0, 0.2, 0.4, 0.6, 0.75, 0.9 | fit (table 3), fig. 4 | f/Re, V, Ω, β |
| MIG-R2 | rigid | 0.4 | 8, 32, 128 | 0.2, 0.4, 0.6, 0.75, 0.9 | fig. 9 | f/Re |
| MIG-R3 | rigid | 0.2, 0.6 | 2, 128 | neutral position (secant on ε) | fit, fig. 10b | ε_eq, V |
| MIG-R4 | rigid | 0.6, 0.8 | 2 | 0.1 | fit (slope at the centre) | sign and value of ∂_ε f |
| MIG-S1 | stress-free | 0.4 | 2 | 0, 0.2, 0.4, 0.6, 0.75, 0.9 | fit (table 4) | f/Re, V, β |
| MIG-S2 | stress-free | 0.2 | 2, 32 | neutral position | fit, fig. 10d (→ √2/4 for small d) | ε_eq, V |

Convergence: MIG-R1 at ε/eps_* = 0.4 and 0.9 and MIG-R2 at Re = 128, ε/eps_* = 0.4 on all three
meshes (Richardson estimate gives the error bar of the other points). Options: the default
(secondOrder, ibmConvection, linking) everywhere; ibmConvection off and mergeAlpha 0 at MIG-R1
ε/eps_* = 0.4, 0.9 and MIG-R2 Re = 128.

Tolerances on the finest mesh (proposed): f/Re within 3 % of max |f/Re| of the curve (0.01 at
d = 0.4), V within 0.5 %, Ω within 1 %, ε_eq within 0.01. β is a small signal: at d = 0.4 rigid,
β α_G = G/G_P − 1 is 3·10⁻⁵ at the centre and 4·10⁻⁴ at 0.9 eps_*, so G and G_0 must be converged
to 10⁻⁷ relative; β within 5 % where |β α_G| > 10⁻⁴, reported elsewhere.

Report: tables of the points against the reference with the error and the Richardson bar; plots
in the form of figures 4a–d (f/Re, V, β, Ω against ε/eps_*, d = 0.4, both interfaces), figure 9
(f/Re at Re = 8–128) and figures 10b/d (neutral positions) with the reference curves.

Cost (estimate): finest mesh ~0.6 M cells half domain; about 40 converged points × 2–3 outer steps
plus the convergence and option runs: a few days on hydra, the linear and Re = 128 sets first.

## 9. Candidate cases for the gaps (proposed 2026-11-05, to be selected)

Literature search for independent references, each case chosen to stress one gap. "Exact" and
"spectral" references allow an order-of-accuracy statement; correlations only plausibility.
The user supplied the papers of the sphere, bubble, rotating and oscillating cylinder on
2026-11-06: these four are specified in §10 (V10–V13); the others are still candidates.

| gap | candidate | stresses | reference | remark |
|---|---|---|---|---|
| K7 at a flat wall | **Hiemenz stagnation flow** onto a plane IBM wall, tilted tan α = 0, 1/3, 1; 3-D: Homann flow | K7 K3 (no curvature, no separation) | exact NS solution: f''' + f f'' − f'² + 1 = 0, f''(0) = 1.232588 (Homann 1.311938); Dirichlet data from it on the outer boundary | 2-D, minutes: a gate (order of u, p, wall shear) |
| K7 at a curved wall, laminar | **DFG 2D-1** (Re 20, steady) and **2D-2** (Re 100, periodic), Schäfer & Turek 1996 | K7 K3 K5 K8, separation, small lift | 2D-1: C_D 5.57953523384, C_L 0.010618948146, Δp 0.11752016697 (Nabh 1998, spectral); 2D-2: official reference file of the FeatFlow page (C_D, C_L max, St ≈ 0.3) | replaces the body-fitted reference of V5 |
| 3-D, single curvature | **DFG 3D-1Z** (Re 20, cylinder in a square channel), 3D-2Z (Re 100) | 3-D K1 K3 K8 | 3D-1Z: c_D 6.18533, c_L 0.009401 (Braack & Richter 2006); 3D-2Z: c_D max 3.2978 (FeatFlow finest level) | channel walls body-fitted, only the cylinder immersed |
| 3-D, double curvature | **sphere in uniform flow**, Re 20–300, Johnson & Patel 1999 (JFM 378) | 3-D K1 K3 K7 K8 | C_D, wake length, separation angle (steady axisymmetric up to Re ≈ 210), steady oblique wake 210–270, shedding at 300 | the tables need the paper |
| K4 free slip, finite Re | **spherical bubble** (stress-free sphere) in uniform flow, Re 1–300 | K4 in 3-D, no separation at all | Hadamard–Rybczynski C_D = 16/Re (Re → 0), Moore C_D = 48/Re (1 − 2.21/√Re) (Re ≫ 1), Mei & Klausner 1992 correlation, Blanco & Magnaudet 1995 (numerical) | acceptance case of the zero-shear fix with MIG-S; a spurious wake is the signature of the missing curvature shear |
| K8 torque, rotation with convection | **rotating cylinder in uniform flow**, Re 200, rotation rate 0–2 (Mittal & Kumar 2003, JFM 476) | K3 with wall velocity, K8 lift and torque, K7 | lift, drag, torque (power), shedding suppressed above α ≈ 1.9 | uses C1, C2 |
| K9 against an independent reference | **in-line oscillating cylinder** in fluid at rest, KC = 5, Re = 100 (Dütsch et al. 1998, JFM 360) | K9 | phase-averaged velocity profiles (LDA and computation), force; C_d 2.09, C_i 1.45 as quoted by later papers | today V3 is checked against the notebook, i.e. the same method |
| K12 K13 between channel and NACA | **NASA TMR 2-D flat plate (2DZP)** and **bump-in-channel (2DBC)**, SA | K13 in a developing layer (flat plate), K12 K13 with mild curvature and pressure gradient (bump) | grid-converged SA results of CFL3D and FUN3D on the TMR grid families | the same model as ours: compares wall treatment, not model |

Not needed from the literature: decomposition (G8) and regression (G7) are self-references.

## 10. Cases SPH, BUB, ROT, OSC (V10–V13, specified 2026-11-06)

Reference data in the repository: `validation/reference_sphere.py`, `reference_bubble.py`,
`reference_rotatingCylinder.py`, `reference_oscillatingCylinder.py` (tables of the papers, the
definitions, `check()`), and the figures that exist only as bitmaps digitised by
`validation/bitmapFigures_extract.py` into `validation/data/johnsonPatel1999/` and
`validation/data/mittalKumar2003/` (only the authors' computed markers; ~0.3 % of the axis range;
checks against the text: C_D(α = 0) 1.318 vs 1.316, St(0) 0.1932 vs 0.1934, |C_L|(Re 212) 0.0235
vs 0.024). The PDFs are not in the repository.

Cross-checks of independent references: sphere C_D at Re 50 and 100, Johnson & Patel (digitised)
1.583, 1.100 against Magnaudet et al. 1.581, 1.092 and Masliyah & Epstein 1.59, 1.08; Moore's
formula reproduces the Moore row of Blanco & Magnaudet table I to three digits and their numerical
bubble drag is within 1.5 % of it for Re ≥ 100; Mittal & Kumar C_M/α = (C_P − C_D)/α² is smooth
(0.238, 0.234, 0.234, 0.234, 0.233 for α = 1…3), so the two digitised curves are consistent.

Common setup: Cartesian background mesh with static refinement (refineMesh, hanging nodes) in
nested boxes around the body, D = 1, U = 1, ν = 1/Re; convergence on three levels of the finest
box, D/h given per case; outer-boundary distance checked once at the most sensitive point;
defaults (secondOrder, ibmConvection, linking), ibmConvection off and mergeAlpha 0 once per case.

### SPH (V10): rigid sphere in uniform flow

Reference: Johnson & Patel 1999 (JFM 378, 19–70): body-fitted FV, outer radius 15 D, order ≈ 1.7;
figures 4, 5, 10 digitised; Re = 300 values from the text. Magnaudet, Rivero & Fabre 1995 and
Masliyah & Epstein 1970 for C_D at Re = 1–100 (table II of Blanco & Magnaudet).

- Domain x ∈ [−10, 20], half domain z ≥ 0 with a symmetry plane through the centre (for Re ≤ 200
  a quarter domain with two planes is equivalent; the half domain keeps the planar symmetry of
  the Re = 250 and 300 wakes, which is the (x, y) plane in the paper), lateral boundaries at 10 D
  slip (checked at 15 D for Re = 50); inflow U e_x, outflow zero gradient / p = 0.
- D/h = 24, 32, 48 (Re ≤ 100), 32, 48, 64 (Re 150–300) in a box [−1, 3] × [0, 1.2]², coarsening
  by 2 per nested box.
- Steady ibmSimpleFoam for Re = 20, 50, 100, 150, 200, 250 (Re 250 started from Re 200 with a
  small lateral perturbation); ibmPimpleFoam for Re = 300 (Δt = 0.01 D/U, 30 shedding periods,
  statistics over the last 15).
- Quantities: C_D; separation angle θ_s and bubble length x_s from the wall shear and the axis
  velocity; vortex centre (x_c, y_c); |C_L| at Re = 250; C_D, C_L mean and amplitude, St at 300.
- Tolerances (finest mesh): C_D within 1.5 % (J&P and Magnaudet et al. differ by 0.8 % at
  Re = 100); θ_s within 1.5°, x_s within 0.03 D, (x_c, y_c) within 0.02 D; |C_L|(250) within
  5 %; Re 300 C_D within 2 %, St within 2 %, |C_L| within 10 %. The onset Re of the non-axisymmetric
  wake (211, Natarajan & Acrivos 210) reported, not judged.
- Stresses: 3-D cut cells on a doubly curved convex wall, separation, convection at a curved wall,
  3-D force integration; the first 3-D flow case, and the rigid part of MIG without the pipe and
  without the small lift. Relies on G1–G3, V2.
- Cost: ~1–3 M cells half domain at the finest level; a day on hydra for the steady set.

### BUB (V11): stress-free sphere (clean bubble of fixed shape) in uniform flow

Reference: Blanco & Magnaudet 1995 (Phys. Fluids 7, 1265), table I (sphere, Re = 50–1000, drag to
< 0.5 % by their grid study), Moore's 48/Re (1 − 2.21/√Re); Hadamard–Rybczynski 16/Re for Re → 0
(in MIG-G2). Exact identity on a stress-free sphere at every Re: surface vorticity ω = 2 u_θ/a
(zero shear ∂u_θ/∂r = u_θ/r, u_r = 0); free slip ∂u_t/∂n = 0 gives u_θ/a.

- Mesh, domain and solver of SPH, `freeSlip` sphere; Re = 50, 100, 200, 300 (steady and
  axisymmetric: a clean spherical bubble has no separation at any Re).
- Quantities: C_D; tangential shear stress on the surface, max |μ(∂u_θ/∂r − u_θ/r)| over the
  sampled wall points, in units of μU/a; reversed flow on the axis behind the bubble.
- Tolerances: C_D within 2 % of table I; wall shear → 0 at the observed order; any recirculation
  is a FAIL. Expected to FAIL until the zero-shear term (`doc/freeSlipPlan.md`, fix 2) is
  implemented: this is its 3-D acceptance case at finite Re (with MIG-S at Re → 0).
- Optional after an ellipsoid body (STL or analytic): oblate spheroids χ = 1.5–1.95, table III
  (curvature-dependent shear, standing eddy for χ ≥ 1.65 in a range of Re).
- Stresses K4 in 3-D at a doubly curved wall, without separation, so the error is the boundary
  condition alone. Relies on SPH.

### ROT (V12): cylinder rotating in a uniform stream

Reference: Mittal & Kumar 2003 (JFM 476, 303–334): stabilised FEM, outer boundary 100 D, 160
points on the circumference, figure 11 digitised (time-averaged C_L, C_D, power coefficient C_P,
r.m.s. of C_L and C_D, St, for α = a ω/U = 0…4.3); text: St(0) = 0.1934, C_D(0) = 1.316, shedding
for α < 1.91, unsteady again for 4.34 < α < 4.70. Torque from C_M = (C_P − C_D)/α.

- 2-D, Re = 2aU/ν = 200, cylinder at the origin, square domain of half width L = 100 D (their
  L; checked with 50 D at α = 2); inflow U, outflow zero stress, slip at the sides, as in the
  paper. D/h = 40, 80, 160 in the finest box.
- α = 0, 0.5, 1, 1.5, 1.9 unsteady (ibmPimpleFoam, rotating cylinder as today, Δt such that
  CFL ≤ 0.5, 40 periods, statistics over the last 20); α = 2, 2.5, 3, 4 steady (ibmSimpleFoam
  with C2, or ibmPimpleFoam to the steady state until C2 exists).
- Quantities: time-averaged C_L, C_D, C_M (needs C1), r.m.s. of C_L and C_D, St; the shedding
  threshold (C_L,rms at α = 1.9 vs 2).
- Tolerances: C_L and C_D means within 2 % or 0.02 absolute; St within 1 %; r.m.s. within 5 %;
  C_M within 5 % for α ≥ 1.5 (the digitised C_P − C_D is too uncertain below); shedding at 1.9
  (C_L,rms 0.069) and steady at 2.
- Stresses the no-slip traction with a tangential wall velocity up to 4U, the torque, the lift
  under convection; the steady part is the 2-D acceptance case of C2. Relies on V2, V5.
- Cost: 2-D, ~10⁵ cells; tier 1 runs α = 0, 1, 2 (hours), the rest in tier 2.

### OSC (V13): cylinder oscillating in-line in fluid at rest

Reference: Dütsch, Durst, Becker & Lienhart 1998 (JFM 360, 249–271): FV in the cylinder frame,
outer boundary at 60 D, grid and time-step study (table 1) and measured LDA profiles. x(t) =
−A sin(2πft), Re = U_max D/ν = 100, KC = 2πA/D = 5; Morison fit of the in-line force F(t) =
−½ρDc_d ẋ|ẋ| − ¼πρD²c_i ẍ gives c_d = 2.09, c_i = 1.45 on the finest set (C); 2.10/1.45 and
2.15/1.42 on the coarser sets. β = Re/KC = 35, KC = 1–4 (table 2, symmetric attached regimes A*,
A): c_d 5.39, 3.72, 2.98, 2.28, 1.97; c_i 1.36, 1.35, 1.34, 1.31, 1.30.

- 2-D, moving body in the inertial frame (ibmPimpleFoam, `oscillating` motion of V3), fluid at
  rest, domain half width 30 D (checked with 60 D), refined box covering the stroke ±A + 1 D;
  D/h = 40, 80, 160; Δt = T/720 (and T/360 once).
- Quantities: c_d, c_i from the least-squares Morison fit
  (`reference_oscillatingCylinder.morison_fit`) over cycles 5–10; force history against figure 10
  as a plot only. The force here is in the inertial frame and compares with c_i (c_m = c_i + 1
  belongs to the accelerated frame).
- Tolerances: c_d and c_i within 2 % (sets B1 and C of the paper differ by 0.5 %); converging with
  the mesh; KC = 1–4 at β = 35 within 3 %.
- Optional: the phase-averaged velocity profiles of figure 8 (x = −0.6, 0, 0.6, 1.2 D, phases
  180°, 210°, 330°) are bitmaps with overlapping markers and not digitised; the paper offers the
  LDA data on request.
- Stresses moving geometry (cells changing status, fresh cells, the force on an accelerating
  body: added mass from the pressure) against a reference computed with a different method.
  Relies on V3 (Galilean consistency).
