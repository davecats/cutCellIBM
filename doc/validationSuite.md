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
| K4 | free-slip wall; zero shear (curvature friction, normal viscous flux), `zeroShear` default on since 2190505 (`doc/freeSlipPlan.md`) | wall coefficients, `zeroShearSource` |
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
(`data/secondOrder/tutorials_secondOrderDefault.txt`). The cases of §9–§11 (adopted 2026-11-06)
fill these gaps.

## 3. Gates (minutes each, all on every change)

| id | check | tolerance |
|---|---|---|
| G1 | geometry: Σ αV against the exact area/volume (cylinder 2-D, sphere 3-D, cylinder through the walls z = 0, L in 3-D straight and tilted (V14), inverted), closure Σ S_k = S_w, linear pressure exact | volume error converging at ≥ 1.8; closure and linear p ≤ 10⁻¹⁰ |
| G2 | operator orders, `ibmOrderTest` LD LN CD CN GL GQ, N = 20…160, mergeAlpha 0 and 0.5 | order ≥ 1.8 (LD LN CD CN GL), ≥ 0.9 (GQ); errors within 1 % of the stored baseline |
| G3 | hydrostatic container | max \|U\| ≤ stored baseline × 1.05 and converging at ≥ 1.8 |
| G4 | free stream around a translating body | max \|U − u_b\| ≤ 10⁻¹⁰ |
| G5 | rotation sources (`rotatingCylinder`, `eccentricRotation`) | \|Σ S\| ≤ 10⁻¹² ; max \|S\| ≤ stored × 1.05 |
| G6 | scalar budgets (`tutorials/scalar`) | closure ≤ 10⁻¹⁰ relative |
| G7 | regression: forces and moments of every tutorial against the stored baseline | ≤ 10⁻⁶ of the body's load scale (larger = intended change, baseline updated in the same commit with a reason) |
| G8 | decomposition: cylinderRe1 and channelSA (linked, wall law) serial vs 4 and 6 ranks, cylinderRe1 on 5 slabs that split linked groups, converged | ≤ 10⁻⁸ relative in the forces |
| G9 | stagnation flow onto a plane wall, exact NS solution (Hiemenz 2-D: tan α = 0 at cut fractions 0.1, 0.5, 0.9 and tan α = 1/3, 1; Homann 3-D: tan α = 0, 1/3), §11 | order ≥ 1.8 (L2 of u and p), ≥ 1.5 (tangential wall shear); errors within 1 % of the stored baseline |

## 4. Validation cases

| id | case | stresses | relies on | reference | tier |
|---|---|---|---|---|---|
| V1 | Taylor–Couette, N = 40…320 | K3 K5 K6 K2, torque | G1 G2 | exact | 1 |
| V2 | periodic cylinder Re_D = 1, no slip and free slip, Richardson | K3 K4 K8 | G1–G3 | body-fitted Richardson (no slip: notebook; free slip: `freeslip_bf`) | 1 |
| V3 | translating and oscillating cylinder | K9 | G4 G5 | Galilean equivalence with V2, notebook | 1 |
| V4 | scalars, static and moving | K11 | G6 | body-fitted steady cylinder scalar | 1 |
| V5 | **new:** DFG 2D-1 (Re 20, steady) and 2D-2 (Re 100, periodic), cylinder in a channel, D/h = 20…160 (§11) | K7 K3 K5 K8, separation, small lift | V2, G9 | Nabh 1998 (spectral); FeatFlow 2D-2 reference files; body-fitted OpenFOAM as a diagnostic | 1 (D/h ≤ 80), 2 |
| V6 | **new:** MIG, sphere in a pipe (below) | 3-D K1 K3 K4 K7 K8 K10 | MIG-G1, MIG-G2, C1 C2 | Rivero-Rodriguez & Scheid 2018 | 2 |
| V7 | SA channel: cut fraction f sweep, tilt tan α, wall-resolved and wall-function, linked and unlinked | K12 K13 | V1 V2 V5 | body-fitted SA channel | 1 (wall-resolved), 2 (sweeps) |
| V8 | NACA 0012, α = 0, 10, wall-resolved and wall-function | integral | everything | body-fitted TMR grid | 2 |
| V9 | AMR (after the merge): V5 and V3 with refinement against the uniform fine mesh | K15 | V3 V5 | the same solver on the uniform mesh | 2 |
| V10 | **new:** SPH, rigid sphere in uniform flow, Re 20–300 (§10) | 3-D K1 K3 K7 K8, separation | G1–G3, V2 | Johnson & Patel 1999; Magnaudet et al. 1995 | 2 |
| V11 | **new:** BUB, stress-free sphere in uniform flow, Re 50–300 (§10) | K4 in 3-D | V10 | Blanco & Magnaudet 1995; Moore; exact surface identity | 2 |
| V12 | **new:** ROT, cylinder rotating in a uniform stream, Re 200, α 0–4 (§10) | K3 with wall velocity, K8 torque, K7 | V2 V5, C1 (C2 for α ≥ 2 steady) | Mittal & Kumar 2003 | 1 (α = 0, 1, 2), 2 |
| V13 | **new:** OSC, cylinder oscillating in-line in fluid at rest, Re 100, KC 5 (§10) | K9, K8 (added mass) | V3 | Dütsch et al. 1998 | 1 |
| V14 | **new:** DFG 3D-1Z (Re 20) and 3D-3Z (Re(t) ≤ 100), cylinder across a square channel (§11) | 3-D K1 K3 K8, immersed wall meeting a body-fitted wall | V5 | Braack & Richter 2006; FeatFlow level 4 | 2 |
| V15 | **new:** TMR flat plate 2DZP, SA, Re 5·10⁶, wall-resolved (cut fraction sweep) and wall function (§11) | K12 K13 in a developing boundary layer | V7 | body-fitted OpenFOAM on the TMR grid; TMR CFL3D/FUN3D | 1 (wall-resolved, f = 0.5), 2 |
| V16 | **new:** TMR bump-in-channel 2DBC, SA, Re 3·10⁶, wall function (§11) | K13 with curvature and pressure gradient | V15 | body-fitted OpenFOAM, same wall function; TMR CFL3D/FUN3D | 2 |

Tier 1: hours in total on hydra, the baseline before the AMR merge. Tier 2: days.

## 5. Options

Every case runs with the defaults (secondOrder, ibmConvection, mergeAlpha 0.5, wall law SA,
wallSampleDistance 1.5, ibmConvectionWallLaw off). Off-default, one at a time, only where the
option acts: mergeAlpha 0 (G2, G9, V1, V2, V5, V7, V15), ibmConvection off (G9, V5, V7, V8,
V15, V16), wall law none (V7, V8, V15), ibmConvectionWallLaw on (V7, V8, V15, V16). Full product
only for interacting options: ibmConvection × mergeAlpha in G9, V5 and V7 (the NACA study showed
they interact). secondOrder false
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
- Run times (`timing.py`, 2026-11-16): per run from the solver logs, host, ranks, fluid cells,
  steps, the master's CPU time, pressure iterations per step and the cost per cell and step
  (core-µs, start-up excluded); the pairs of immersed and body-fitted runs of the same problem
  (G9, V2, V4, V5, V7, V12, V15). `latexreport.py <commit> V2=<commit>` writes
  `tables/timings*.tex` and `figures/cost_accuracy` (CPU time against the error). Stored times
  carry the load of their host: compare only runs on one host at one time; the transient pairs
  were re-timed on one host (`~/runs/timing_bench/bench.sh`, finding 16).

## 7. Order of work

1. Driver and gates G1–G9, tier 1 cases (V1–V5, V7 wall-resolved, V12 α = 0, 1, V13, V15
   wall-resolved f = 0.5) on the current branch: the baseline.
2. Merge `feature_amr`; accepted when the gates pass and tier 1 agrees with the baseline.
3. C1, C2 and the free-slip zero-shear fix, each accepted by its gates; then tier 2 (V10 SPH and
   V11 BUB before V6 MIG, V14 before V10, V7 sweeps, V15 complete and V16 before V8, V9, V12
   complete). C1, C2 done; zero shear done (`zeroShear`, default on since 2190505, 2026-11-15,
   §8 prerequisites).

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
  of each body, in `forces.dat` and the log. **Done (2026-11-10):** in `moments.dat` next to
  `forces.dat` (same columns; the readers of `forces.dat` assume six columns per body) and the
  log. It needed the transpose stress of a rotating wall (finding 10). V1 checks the exact
  Taylor-Couette torque on both cylinders, V12 the mean C_M against the body-fitted solution.
- **C2 tangential wall velocity on a static geometry:** a body whose motion leaves it invariant
  (sphere rotating about a diameter, body of revolution about its axis, cylinder sliding along its
  axis) keeps its geometry fixed, has no wall flux, and only supplies a wall velocity to the
  no-slip traction. Today such a body counts as moving and ibmSimpleFoam refuses it. The
  invariance is checked at start-up (level set at the cut points after a virtual small motion),
  fatal otherwise. Used here for the rotating sphere and the sliding pipe.
  **Done (8108c14, 2026-11-14):** opt-in key `invariant yes` in the motion dictionary (user
  decision), any motion type; checked at start-up (virtual motion moving the wall by one cell)
  and at every step (actual state), |Δφ|/h at the edge crossings ≤ `invarianceTolerance` (1e-6);
  fatal for alpha-field bodies. G1–G9 bit-identical to 8d416fd. Acceptance (scripts and runs in
  `~/runs/c2`): Taylor-Couette invariant = moving to 4–5 digits (ibmPimpleFoam), ibmSimpleFoam
  torque equal to 1e-6; sliding pipe (MIG-G1 prototype, Re 2, V = 0.7) = static pipe shifted by
  −V to round-off, u L2 and G order 2.0 (G −1.18, −0.28, −0.070 % at h = 1/16, 1/32, 1/64: G
  misses 0.2 % at the base h = 1/32, refine at the pipe wall); spherical Couette (user decision
  2026-11-14, in place of the periodic-box MIG-G2 for C2: R = 0.5 in Ro = 1, Stokes, exact
  torque 8πμΩR³Ro³/(Ro³ − R³)) torque −1.27, −0.30, −0.13 % at d/h = 16, 32, 48, order 2.07, u L2
  order 2.05; V12 α = 2 steady (ibmSimpleFoam, U relaxed 0.7; 0.9 cycles at D/h = 40, where the
  steady state is only weakly damped) C_L −5.124, −5.315, −5.360, C_D 0.3040, 0.2924, 0.2897,
  C_M 0.4922, 0.4709, 0.4646 at D/h = 40, 80, 160 against the transient means at t = 100–150
  −5.146, −5.321, −5.359; 0.3106, 0.2924, 0.2915; 0.4871, 0.4688, 0.4645 (body-fitted −5.375,
  0.2887, 0.4633). Non-invariant motions (sideways, tilted by 1e-3, rotation about a non-axis, a
  prescribed motion turning sideways later) are fatal. The periodic-box MIG-G2 (translating,
  stress-free and rotating sphere) stays a V6 gate, run with C3.
- **Zero-shear free slip** (`doc/freeSlipPlan.md`): **done; `zeroShear` default true since
  2190505 (f4f2366, 361c23c, 8d60565 on 2026-11-14 behind the switch; default and baselines
  2026-11-15)**. Wall flux of the Laplacian
  form ν (A_w/V)(H_w·w_t,w + ω × n + (∂u_n/∂n) n) on free-slip walls of smooth level-set bodies:
  curvature tensor from the level-set Hessian, Robin factor (I + d H_w)⁻¹ for the wall value,
  ∂u_n/∂n = −tr(P∇uP) explicit. G1–G9 bit-identical with the switch off. Exact rigid rotations
  (Taylor–Couette convex/concave, transient and steady; 3-D spherical shell) and a manufactured
  solution with a varying slip and a normal velocity gradient converge at second order; V2 free
  slip linked 3.9641, 3.9895, 3.9974, 3.9990 (N = 41 … 321) against the body-fitted stress-free
  limit ≈ 3.999 (−24.8 % before). Default on: G7 and V2 re-accepted from 2190505 (§12
  amendments). Open: unlinked with alphaMin > 0 (finding 13), triSurface and alpha-field
  bodies, V11.
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
| MIG-R1 | rigid | 0.4 | 2 (1) | 0, 0.2, 0.4, 0.6, 0.75, 0.9 | fit (table 3); β: fig. 4c | f/Re, V, Ω, β |
| MIG-R2 | rigid | 0.4 | 8, 32, 128 | 0.2, 0.4, 0.6, 0.75, 0.9 | fig. 9 | f/Re |
| MIG-R3 | rigid | 0.2, 0.6 | 2, 128 | neutral position (secant on ε) | fit, fig. 10b | ε_eq, V |
| MIG-R4 | rigid | 0.6, 0.8 | 2 | 0.1 | fit (slope at the centre) | sign and value of ∂_ε f |
| MIG-S1 | stress-free | 0.4 | 2 | 0, 0.2, 0.4, 0.6, 0.75, 0.9 | fit (table 4); β reported | f/Re, V, β |
| MIG-S2 | stress-free | 0.2 | 2, 32 | neutral position | fit, fig. 10d (→ √2/4 for small d) | ε_eq, V |

Convergence: MIG-R1 at ε/eps_* = 0.4 and 0.9 and MIG-R2 at Re = 128, ε/eps_* = 0.4 on all three
meshes (Richardson estimate gives the error bar of the other points). Options: the default
(secondOrder, ibmConvection, linking) everywhere; ibmConvection off and mergeAlpha 0 at MIG-R1
ε/eps_* = 0.4, 0.9 and MIG-R2 Re = 128.

Tolerances on the finest mesh (proposed): f/Re within 3 % of max |f/Re| of the curve (0.01 at
d = 0.4), V within 0.5 %, Ω within 1 %, ε_eq within 0.01. β within 5 % where |β α_G| > 10⁻⁴,
reported elsewhere; G and G_0 converged to 10⁻⁷ relative.
β reference (2026-11-17, from the arXiv text of the paper, eq. 2.7 and figures 4c, 7c): the
column "100 β0" of the fit tables holds β0 itself (β = O(1); the earlier reading β0 = value/100
gave β α_G = 3·10⁻⁵…4·10⁻⁴ here, in fact it is 0.003–0.05), and the rigid fit lies below the
paper's own figure 4c by a factor 1.22–1.24 up to 0.75 eps_* (1.31 at 0.9 eps_*), while V0 and
Ω0 agree with figures 4b, 4d to 0.3 %. MIG-R1 judges β on figure 4c
(`data/riveroScheid2018/fig4_rigid_d0.4.dat`, digitised from the arXiv PDF); for MIG-S1 the
stress-free fit is the only reference at d = 0.4 (figure 8c is a colour map), so β is reported
there. First point (rigid, d = 0.4, 0.4 eps_*, Re 2, d/h 16/32/48): β 0.601, 0.611, 0.614
(extrapolated 0.617) against 0.613 (figure) and 0.496 (fit); f/Re 0.318, 0.300, 0.296
(extrapolated 0.292, fit 0.293); V 1.6660 (fit 1.6647), Ω 0.9165 (fit 0.9129, figure 0.9145).

Report: tables of the points against the reference with the error and the Richardson bar; plots
in the form of figures 4a–d (f/Re, V, β, Ω against ε/eps_*, d = 0.4, both interfaces), figure 9
(f/Re at Re = 8–128) and figures 10b/d (neutral positions) with the reference curves.

Cost (estimate): finest mesh ~0.6 M cells half domain; about 40 converged points × 2–3 outer steps
plus the convergence and option runs: a few days on hydra, the linear and Re = 128 sets first.

## 9. Candidate cases for the gaps (proposed 2026-11-05, all adopted 2026-11-06)

Literature search for independent references, each case chosen to stress one gap. "Exact" and
"spectral" references allow an order-of-accuracy statement; correlations only plausibility.
The user supplied the papers of the sphere, bubble, rotating and oscillating cylinder on
2026-11-06: these four are specified in §10 (V10–V13). The user adopted the others on the same
day: Hiemenz/Homann as gate G9, DFG 2-D as the new V5, DFG 3-D as V14, TMR 2DZP and 2DBC as V15
and V16, specified in §11.

| gap | candidate | stresses | reference | remark |
|---|---|---|---|---|
| K7 at a flat wall | **Hiemenz stagnation flow** onto a plane IBM wall, tilted tan α = 0, 1/3, 1; 3-D: Homann flow | K7 K3 (no curvature, no separation) | exact NS solution: f''' + f f'' − f'² + 1 = 0, f''(0) = 1.232588 (Homann 1.311938); Dirichlet data from it on the outer boundary | 2-D, minutes: a gate (order of u, p, wall shear) |
| K7 at a curved wall, laminar | **DFG 2D-1** (Re 20, steady) and **2D-2** (Re 100, periodic), Schäfer & Turek 1996 | K7 K3 K5 K8, separation, small lift | 2D-1: C_D 5.57953523384, C_L 0.010618948146, Δp 0.11752016697 (Nabh 1998, spectral); 2D-2: official reference file of the FeatFlow page (C_D, C_L max, St ≈ 0.3) | replaces the body-fitted reference of V5 |
| 3-D, single curvature | **DFG 3D-1Z** (Re 20, cylinder in a square channel), 3D-3Z (Re(t) ≤ 100) | 3-D K1 K3 K8 | 3D-1Z: c_D 6.18533, c_L 0.009401 (Braack & Richter 2006); 3D-3Z: c_D max 3.2978 (FeatFlow level 4; the periodic 3D-2Z has no reference) | channel walls body-fitted, only the cylinder immersed |
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
- As implemented (`validation/suite/v10.py`, 2026-11-16): background h0 = 16/(D/h) (cubic cells on
  the spec domain), four boxes fixed in space, [−5.2, 11] × [−5.4, 5.4] × [0, 5.4], [−2.8, 6.5] ×
  [−3, 3] × [0, 3], [−1.6, 4.2] × [−1.8, 1.8] × [0, 1.8], [−1, 3] × [−1.2, 1.2] × [0, 1.2] (~7
  cells of the coarser level between them at D/h 24): 0.27, 0.64, 2.2, 5.1 M cells at D/h 24, 32,
  48, 64. Re 230 on D/h 64 for the onset (reported against the same two-point estimate from J&P's
  |C_L|: linear |C_L|² gives ~180 from their values too). Once at D/h 32: full domain (symmetry
  plane through the cut cells: +8e-5 in C_D at D/h 16), lateral 15 D, ibmConvection off,
  mergeAlpha 0. θ_s from the wall shear the solver applies (ν c V u_P/|S_w| per wall cell,
  area-averaged in 2° bands of θ, spread over 1–3° bands reported): the wall cells cluster on the
  grid planes x = const, ~h/(R sin θ) apart; least-squares u_θ at h and 2h was 1–1.5° low and did
  not converge. Re 50, 100 at D/h 48: θ_s 139.20°, 126.90° (J&P 139.39°, 127.03°).
- Cost: 5.1 M cells at D/h 64; the D/h 64 chain on HoreKa (4 nodes: SIMPLE 0.46 s per iteration,
  Re 250 converged in 4036 iterations; Re 300 4 s per step), D/h ≤ 48 on hydra.

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
- As implemented (`validation/suite/v11.py`, 2026-11-17, code 2190505): `v10.setup` with
  `wall="freeSlip"`; D/h 24, 32, 48 for all Re (no separation), 8, 16, 32 ranks. Surface stress
  and vorticity from an alpha-weighted quadratic fit of u within 4 h of points on the sphere
  (θ 15–165°, four meridians). User decision (2026-11-17): C_D judged on the Richardson
  extrapolation within 2 % with an observed order ≥ 1.5 (D/h 48 does not resolve the surface
  layer: +8.5 % at Re 300; 2 % on one mesh would need D/h ≈ 128, ~40 M cells); the finest mesh
  reported; the tangential stress must decrease with the mesh. Once at D/h 32: Re 100
  ibmConvection off, mergeAlpha 0, zeroShear off; Re 50 lateral 15 D.

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
  with C2, `invariant yes`; SIMPLEC with U relaxed 0.7: 0.9 cycles at D/h = 40).
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

## 11. Gate G9 and cases DFG, TMR (V5, V14–V16, specified 2026-11-06)

Reference data in the repository, all from public sources:

- `validation/reference_hiemenz.py`: Hiemenz and Homann solutions by `solve_bvp` (tolerance 10⁻¹⁰),
  `exact(points, a, nu, origin, normal)` for any wall orientation; `check()` gives f''(0) =
  1.232588 and 1.311938 and displacement constants 0.6479 and 0.5689, as in the literature, and a
  finite-difference residual of the momentum equations of 10⁻⁸ on a tilted wall.
- `validation/reference_dfg.py` with `validation/data/dfg/`: 2D-1 (Nabh 1998), 2D-2 (FeatFlow
  "new reference results" 2014, levels 4–6, one cycle of level 6), 3D-1Z (Braack & Richter 2006;
  FeatFlow levels 1–4), 3D-3Z (FeatFlow level 4 time series), `cycle_statistics()` for the 2D-2
  quantities of our force history.
- `validation/reference_tmr.py` with `validation/data/tmr/` (TMR files unchanged, README with the
  sources): convergence of C_f and forces on the five TMR grids, C_f(x), C_p(x), profiles, ν_t.

Cross-checks: 2D-2 levels 5 and 6 agree to 0.03 % in C_D,max and 0.01 % in St; their Δp(t0 + T/2)
= 2.483 is inside the 1996 interval 2.46–2.50, while C_L,max = 0.987 is just below the 1996 interval
0.99–1.01, so the new files, not the intervals, are the reference. FeatFlow level 4 of 3D-1Z
agrees with Braack & Richter to 0.01 % (C_D) and 0.04 % (C_L). The observed orders of the TMR
flat-plate C_f(0.97) from the files (1.98 CFL3D, 1.34 FUN3D) reproduce the TMR uncertainty table.

### G9: stagnation flow onto a plane wall (Hiemenz 2-D, Homann 3-D)

Exact steady NS solutions with a wall-normal inflow, no curvature and no separation: the
convection at the wall (K7) and the wall traction (K3) at a flat wall in isolation, where V2 has
no convection and the NACA case cannot separate K7 from K12/K13.

- a = 1, ν = 0.01 (99 % boundary-layer thickness ≈ 0.24 in 2-D), `plane` body (noSlip) through the
  origin, rho = 1. 2-D domain [−1, 1] × [−0.5, 1], N = 20, 40, 80, 160 cells per unit length
  (5–40 cells across the boundary layer); tan α = 0 with the wall at cut fractions f = 0.1, 0.5, 0.9
  of the first cell, and tan α = 1/3, 1 (normal (−sin α, cos α)).
- 3-D (Homann): [−0.5, 0.5]² × [−0.4, 0.6], N = 20, 40, 80, tan α = 0 (f = 0.5) and 1/3 (tilted
  about x).
- Outer boundary: U fixedValue from `exact()` on every face, p fixedFluxPressure with the exact
  value at one reference cell; steady ibmSimpleFoam (SIMPLEC), converged to 10⁻¹⁰.
- Body-fitted control (2-D, tan α = 0, wall as a mesh patch, same solver without bodies): the
  error of the outer boundary treatment and of the interior scheme alone, against which the IBM
  error is read.
- Metrics: L2 (α-weighted) and L∞ errors of u and p in the fluid cells; L2 error of the wall shear
  in the cut cells against a s √(aν) f''(0); the wall force against the exact integral of
  p0 − a²s²/2 and of the shear over the wall segment inside the domain.
- Tolerances: order ≥ 1.8 in L2 of u and p, ≥ 0.9 in the wall shear (as GQ in G2), errors within
  1 % of the stored baseline; L∞ and the ratio IBM/body-fitted error reported.
- Options: ibmConvection on/off × mergeAlpha 0/0.5 at every N (the place where K7 is seen alone).
- Cost: 2-D 77 k cells at N = 160 (minutes); 3-D 0.5 M cells at N = 80 (≈ 15 min on 8 cores).

### V5 DFG 2D-1 and 2D-2: cylinder in a channel

Reference: Schäfer & Turek 1996; 2D-1 C_D = 5.57953523384, C_L = 0.010618948146, Δp =
0.11752016697 (Nabh 1998, spectral), recirculation length L_a in [0.0842, 0.0852]; 2D-2 (FeatFlow
level 6): St 0.30184, C_D 3.16426–3.22739, C_L −1.02129–0.98657, Δp(t0 + T/2) 2.48322, and the
cycle C_D(t), C_L(t), Δp(t). Definitions in `reference_dfg.py` (U_mean = 2/3 U_m, D = 0.1).

- Channel walls and inflow/outflow are mesh patches (body-fitted), the cylinder is a `cylinder`
  body (noSlip) at (0.2, 0.2): only the cylinder is immersed. Inflow parabolic profile; outflow U
  zeroGradient, p = 0 (the do-nothing condition of the benchmark agrees for a developed profile;
  Δp is a difference and does not see the level).
- Mesh: background h = 0.01 (220 × 41), refinement by 2 in [0.05, 1.2] × [0, 0.41] (the channel
  walls and the near wake) and in nested boxes around the cylinder; D/h = 20, 40, 80 (tier 1), 160
  (tier 2) at the cylinder. The cylinder centre is not on a mesh line: the asymmetry that produces
  C_L is the benchmark's, and errors of the cut geometry show up in the small lift.
- 2D-1: ibmSimpleFoam to 10⁻¹⁰. 2D-2: ibmPimpleFoam from rest, Courant ≤ 0.5, to t = 30, the last
  cycle between minima of C_L evaluated by `cycle_statistics()`; Δt halved once at D/h = 80.
- Quantities: 2D-1 C_D, C_L, Δp, L_a (zero of u on the axis y = 0.2 behind the cylinder, minus
  0.25); 2D-2 St, C_D and C_L extremes, Δp(t0 + T/2), and the cycle against the FeatFlow cycle
  (plot).
- Tolerances at D/h = 80 (tier 1): 2D-1 C_D and Δp within 0.5 %, C_L within 5 %, L_a inside the
  interval, observed order of C_D ≥ 1.5 from 20, 40, 80; 2D-2 St within 0.5 %, C_D extremes within
  0.5 %, C_L extremes within 1 %, Δp within 1 %. At D/h = 160 (tier 2): C_D, Δp, St within 0.2 %.
- Diagnostic, not judged: pimpleFoam/simpleFoam on a body-fitted O-grid of the same resolution
  (the share of the deviation that belongs to the finite-volume scheme rather than to the IBM).
- Options: ibmConvection × mergeAlpha full product (§5).
- Stresses convection, traction and pressure at a curved wall with separation, the force
  integration, and the cut geometry through the lift. Relies on V2 and G9.
- Cost: 2-D, ≈ 0.1 M cells at D/h = 80; 2D-2 to t = 30 a few hours on 8 cores.

### V14 DFG 3D-1Z and 3D-3Z: cylinder across a square channel

Reference: 3D-1Z c_D = 6.18533, c_L = 0.009401 (Braack & Richter 2006); FeatFlow levels 1–4
(OpenFOAM in the same study, body-fitted, level 4: 0.06 % in c_D, 3.5 % in c_L). 3D-3Z: FeatFlow
level 4 c_D(t), c_L(t), c_D,max 3.2978, c_L,min −0.010999. c = 2F/(U_mean² D H), U_mean =
4/9 U_m.

- Channel [0, 2.5] × [0, 0.41]², all four walls body-fitted, cylinder axis z through (0.5, 0.2)
  spanning the full height: the immersed wall meets the body-fitted walls z = 0 and z = H. Cells
  with both a wall patch face and a cut face exist only here: before the run, G1 is extended by a
  cylinder crossing the domain boundary (Σ αV, closure Σ S_k = S_w in those cells).
- Mesh (as implemented, user 2026-11-15, finding 14): background H/n in y, z (n = 16, 24, 32, 48;
  2.5/round(2.5 n/H) in x), one level in [0.25, 1.4] × [0, H] × [0, H] and one more in the finest
  box [0.35, 0.85] × [0.05, 0.36] × [0, H], both symmetric about the channel centre: D/h = 15.6,
  23.4, 31.2, 46.8 (0.32, 1.06, 2.5, 8.5 M cells). Spec before: finest box [0.4, 0.75] × [0.1, 0.3]
  × [0, H], D/h = 16, 24, 32, 48 on HoreKa if the order is unclear.
- 3D-1Z steady ibmSimpleFoam; 3D-3Z ibmPimpleFoam t ∈ [0, 8], inflow 2.25 sin(πt/8), fixed Δt from
  Courant ≤ 0.5 at |u| = 3.75 (checked from the log; |u| = 3 gave Co 0.60), on D/h = 15.6, 23.4,
  46.8, and Δt/2 on 15.6 (user 2026-11-15: the Δt = 0.005 written before is Courant 2.3–3.5;
  D/h 46.8 user 2026-11-16, finding 15).
- Quantities: c_D, c_L; Δp between (0.45, 0.2, 0.205) and (0.55, 0.2, 0.205) reported (no sourced
  reference); 3D-3Z c_D(t), c_L(t), their extremes and the L2 difference to the FeatFlow series.
- Tolerances: 3D-1Z c_D within 1 %, c_L within 10 % at D/h = 46.8 (user 2026-11-15, finding 14;
  before: c_L within 5 % at D/h = 32); 3D-3Z c_D,max within 1 %, c_L,min within 10 % at D/h = 46.8
  (user 2026-11-16; FeatFlow level 2: 0.9 %, OpenFOAM level 3: 7.3 %), c_D,max order REPORT.
- Stresses 3-D cut cells on a singly curved wall, the junction of an immersed and a body-fitted
  wall, 3-D forces. Relies on V5; runs before V10 SPH (curvature in one direction before two).
- Cost: tier 2; 3D-3Z 6000 (D/h 15.6) to 28000 (D/h 46.8) steps, the finest on 8.5 M cells,
  128 ranks, about a day.

### V15 TMR 2DZP: turbulent flat plate (SA)

Reference: TMR SA verification (compressible M = 0.2: C_f(0.97) = 0.002706, C_D = 0.00286 on the
finest grid; incompressible FUN3D and SC/Tetra: C_f(0.97) = 0.002729, 0.9 % higher). TMR warns that
the SA of OpenFOAM (2014) did not match its equations: the primary reference is therefore
body-fitted OpenFOAM (simpleFoam + SpalartAllmaras, the model of ibmSpalartAllmaras) on the TMR
grids 137 × 97, 273 × 193, 545 × 385, and the TMR values check that reference.

- U = 1, ν = 2·10⁻⁷ (Re 5·10⁶ per unit length), plate x ∈ [0, 2], domain [−0.33333, 2] × [0, 1];
  inflow U = 1 and ν̃ = 3ν, top slip (TMR: moving it to 0.48 changes results < 0.2 %), outflow
  p = 0.
- IBM geometry on the TMR grid itself (tensor-product x_i, y_j): one copy of the first cell layer
  added below y = 0 and the grid shifted up by f Δy₁, so that the plate, a `box` body (noSlip) for
  x ≥ 0 below y = 0, cuts that cell with fluid fraction f = 0.25, 0.5, 0.75. The mesh bottom is a
  symmetry patch for x < 0, (1 − f)Δy₁ (≤ 2·10⁻⁶) below the plane of the plate. Body-fitted and IBM
  runs then share every mesh line but the wall row: the comparison isolates the wall treatment.
  Variant once: a freeSlip `box` for x < 0 at the height of the plate (two touching bodies; needs
  the G1 check that their common face is not a wall).
- Wall-resolved (y⁺ of the first fluid cell centre ≤ 1, wall law none and SA law) and wall
  function (the same grid coarsened in y to y⁺ ≈ 30–100 at the sample, SA law, compared with
  body-fitted Spalding/SA law at the same first-cell height). Tilted plate (STL body, freestream
  rotated by α, tan α = 1/10) with the wall function only, tier 2.
- Quantities: C_f(x), C_f(0.97), C_D (one side, reference length 2), u⁺(y⁺) at x = 0.97008 and
  1.90334, ν_t/ν at x = 0.97, max ν_t/ν(x).
- Tolerances, IBM against body-fitted on the same grid: wall-resolved C_f(0.97) and C_D within
  1 %, C_f(x) within 2 % on x ∈ [0.1, 1.9], max ν_t/ν within 2 %, spread over f within 1 %; wall
  function C_f within 3 %. Body-fitted against TMR incompressible C_f(0.97): REPORT (a model
  implementation difference, not an IBM error).
- Stresses the SA terms with the IBM wall distance and the wall law in a boundary layer whose
  thickness and sample y⁺ change along the wall (the channel V7 has neither), and the
  leading-edge singularity in a cut cell (excluded from the tolerance by x ≥ 0.1). Relies on V7.
- Cost: 545 × 385 = 0.2 M cells, steady SA, ≈ 1 h; tier 1 runs 273 × 193, f = 0.5, wall-resolved,
  IBM and body-fitted.

### V16 TMR 2DBC: turbulent bump in a channel (SA)

Reference: TMR SA verification on grids up to 1409 × 641: C_L = 0.02490, C_D = 0.003572 (C_Dp
0.000382, C_Dv 0.003191), C_f(0.6321975, 0.75, 0.8678025) = 0.005185, 0.006149, 0.002678 (CFL3D;
FUN3D within 0.3 %), C_p(x), C_f(x), u profiles at x = 0.75 and 1.20148, ν_t/ν at 0.75 (all
compressible M = 0.2).

- U = 1, ν = 1/(3·10⁶), wall x ∈ [0, 1.5] with the bump y = 0.05 sin⁴(πx/0.9 − π/3) on
  0.3 ≤ x ≤ 1.2 as an STL body (2-D extruded; leading edge as in V15), symmetry on the bottom for
  x ∈ [−25, 0) and (1.5, 26.5], top symmetry at y = 5, inflow ν̃ = 3ν, outflow p = 0.
- Wall-resolved IBM would need cells of 10⁻⁶ height following the bump (y⁺ ≤ 1 over a bump 0.05
  high): not feasible on a Cartesian mesh, so wall function only. Band around the wall
  x ∈ [−0.1, 1.6], 0 ≤ y − y_bump ≤ 0.08 with Δy = 4·10⁻⁴ (y⁺ ≈ 30–100 at the sample), Δx = 2·10⁻³,
  and Δy halved and doubled once (three meshes).
- Body-fitted reference: OpenFOAM with the same wall function on the TMR grid coarsened at the wall
  to the same first-cell height (as V8); body-fitted wall-resolved on the TMR 705 × 321 grid
  against the TMR values (REPORT, model implementation).
- Quantities: C_p(x), C_f(x), C_f at the three stations, C_L, C_Dp, C_Dv, u profiles, ν_t/ν at 0.75.
- Tolerances, IBM wall function against body-fitted wall function: C_p within 0.005 on
  x ∈ [0.1, 1.4], C_f at the three stations within 3 %, C_L within 2 %; C_Dp (≈ 10 % of the drag)
  reported. Against TMR wall-resolved: REPORT (wall-function model error).
- Stresses the wall law under curvature and favourable/adverse pressure gradient (the sample
  direction and distance on a curved wall, K13), between the flat plate and the NACA. Relies on
  V15.
- Cost: ≈ 0.2 M cells, steady SA, hours; tier 2.

## 12. Implementation (step 7.1, started 2026-11-06)

Driver `validation/suite/suite.py` (`list`, `run <ids>|gates|tier1 [--jobs n] [--quick]`,
`run <id> --reuse <case root>` to re-evaluate finished runs, `accept <ids>`, `report`); one
module per gate or case (`g1.py` … `v15.py`) with its setup, metrics, reference, tolerances,
the components stressed and relied on, and a regression tolerance against the accepted baseline
(`validation/suite/baseline/`). Results are keyed by the commit of the code under test (the last
commit touching `src`, `applications`, `tutorials`; suite-only commits keep the key), JSON in
`validation/suite/results/<commit>/`, case directories in `~/runs/suite/<commit>/<id>/` (never
overwritten: a rerun goes to `<id>.1`, …). Case copies take only the git-tracked files of a
tutorial (its directory may hold the output of a local run).

Errors against exact or reference fields are evaluated at the **value points** of the cells
(cell centre + `ibmValueOffset`: the group centroid of a linked cell; `ibmCentroidOffset` when
unlinked). Evaluating a linked cell at its own centroid shows an O(h) error that is not there
(Taylor–Couette linked u: order 1.39 at the centroids, 2.05 at the value points).

### Amendments to the tolerances of §3–§4 (after the first runs)

- **G1**: the fluid volume comes from sub-sampling (quadrature), not a discretisation:
  bounded (≤ 1e-4 of the body volume; observed 8e-5 2-D, 1.3e-5 3-D) instead of an order; the
  order ≥ 1.8 is required of the wall area (observed 1.9–4.0).
- **G2**: GL and GQ unlinked only (with linking a group has one gradient, the per-cell test of
  `ibmOrderTest` is undefined: a linking-aware GL check is missing; the linked pressure gradient
  is tested in the flow by G3 and V1). GQ (first order with geometric noise, pair orders 0.93,
  1.01, 0.73) fitted over all four meshes, ≥ 0.8.
- **G6**: the budget is judged per step; a step whose scalar solve hits maxIter is reported
  apart (`galilean_zeroFlux`, t = 11.925, the known case of `doc/secondOrderStudy.md`).
  Regression (2026-11-11, user decision): the budget errors are round-off (10⁻¹²–10⁻¹⁰), so a
  relative tolerance flags noise (2b1adc3: steady 9.1e-11 -> 4.1e-11 = 9 units at 5 % + 1e-12);
  tolerance 5 % + a tenth of the metric's gate (10⁻¹⁰ transient, 10⁻⁹ steady).
- **G7** (2026-11-11, user decision): force and moment per body, |F − F_base| ≤ 10⁻⁶
  max(|F_base|, |M_base|/L), |M − M_base| ≤ 10⁻⁶ max(|M_base|, L |F_base|), L the largest distance
  of the body's wall cells from their area centroid. The force-only relative check failed on
  Taylor–Couette, whose net force is a cancellation (|F| ≈ 4e-8 against a torque of 0.28: 2b1adc3
  changed it by 5e-13 = 1.3e-5 relative, 3e-12 of the load scale). Baselines from before store no
  moments: the moment is then reported, and this run's moment sets the force scale.
  G7 re-accepted from 8d416fd (2026-11-14, user decision) to store the moments; its forces equal
  the e7d742a baseline within 1e-8 of the load scale.
- **V12** (2026-11-12, user decision): St is regressed only where the wake sheds (α < 1.9, the
  rule of the St reference check); at α = 2 the lift decays (C_Lrms 0.005) and its dominant
  frequency jumped 0.364 -> 0.410 between 31d1827 and 2b1adc3 with C_Lrms equal to 4e-5.
- **V5, V12, V13 re-accepted from 8d416fd** (2026-11-14, user decision). At 5d9caa5, V12 α = 2
  D/h = 40 depended on the decomposition (C_Lrms 0.148 with the old one, 0.020 with the link
  constraint). The flux of a processor face between two linked-group masters is stored on both
  ranks; the projection cannot remove a difference between the two copies, and
  `cutCellGeometry::ddtCorr` fed it back with gain interpolate(rAU)/Δt > 1 (rAU ≈ Δt/(1.5 ᾱ) in a
  merged small-α group), so it grew from round-off by ~10 % per step. 8d416fd puts the old fluid
  fraction (group mean) in ddtCorr, the coefficient of U⁰ in ddt(α, U) as ρ in the compressible
  ddtCorr(ρ, U, φ), and makes the flux of processor and cyclic faces single-valued after the
  pressure correction (`syncCoupledFlux`, ibmPimpleFoam and ibmSimpleFoam). α = 2 D/h = 40 to
  t = 150: serial, old and constrained decomposition agree (C_D 0.31028, C_L −5.14907, C_Lrms
  0.04917). All transient results change at the discretisation level, every coarse-mesh value
  towards the reference: V5 2D-2 D/h = 40 C_Lmax 0.9946 -> 0.9904 (0.9866); V12 α = 1 D/h = 40
  C_D 1.140 -> 1.123 (1.084), α = 2 D/h = 40 C_D 0.355 -> 0.311, D/h = 80 0.2934 -> 0.2924 (0.29);
  V13 D/h = 40 c_d 2.141 -> 2.121 (2.09). Exception: V12 α = 1 D/h = 40 St 0.1923 -> 0.1928
  (0.1921). Changes are 1–2 % at D/h = 40 and below 1 % on the finer meshes; every check against
  the references passes. Gates, V1–V4, V7 and V15 pass against their baselines (steady and
  G7 tutorial results unchanged within 10⁻⁸).
- **G7 and V2 re-accepted from 2190505** (2026-11-15, user decision; 2194276): zero shear on by
  default. V2 free slip 2.8403, 2.9155, 2.9877, 3.0108 -> 3.9641, 3.9895, 3.9974, 3.9990
  (N = 41 … 321; body-fitted stress-free 4.00318 at M = 256, −0.105 % at N = 321, observed order
  2.4), no slip unchanged within 1e-12; the free-slip V2 check is now against the stress-free
  drag with the no-slip tolerance 0.3 %. G7: `cylinderRe1_pimple_freeSlip` Fx 0.3593 -> 0.4924
  (+37 %), `twoCylinders/smallFreeSlip` 0.1972 -> 0.2736 (+39 %) and its no-slip neighbour
  `bigNoSlip` 0.6425 -> 0.6243 (−2.8 %: Ubar fixed, driving gradient 0.07136 -> 0.07630, the
  summed Fx equals gradient × fluid volume before and after); the spanwise moment of the
  symmetric free-slip cylinder −3.4e-16 -> −3.9e-10 (1e-8 of its scale, solver tolerance); 11
  out-of-plane round-off values (≤ 1e-18) from 361c23c. G1–G6, G8, G9 bit-identical.
- **G8**: cylinder drag (the lift is zero by symmetry, solver noise) and the channel driving
  gradient u_τ², each run required stationary (the gradient changes ≤ 1e-10 over the last 100
  iterations): `residualControl` cannot trigger in this 1-D problem (the Uy and p residuals are
  normalised by ≈ 0 and stay O(0.1) at the converged state). Same criterion in V7.
- **G9**: u L2 over all fluid cells and p L2 in the interior (farther than 0.1 from the outer
  faces) ≥ 1.8; the pressure where the immersed wall meets the outer boundary, max norms and the
  options are reported. Outer boundary: exact U on every face (outflow faces scaled by 1 + O(h²)
  so that the θ-weighted boundary flux closes), p zeroGradient (fixedFluxPressure is evaluated
  before the solver sets it).
- **V1**: u L2 ≥ 1.8 linked and unlinked; p reported (the cut-cell pressure converges below
  second order, open: `doc/secondOrderStudy.md`).
- **V2**: no slip within 0.3 % of the body-fitted M = 256 drag (5.83958); free slip (zero
  shear, the default since 2190505) within 0.3 % of the body-fitted stress-free drag (4.00318,
  M = 256); before, reported against slip (3.04843) and zero shear. Richardson lines are nan:
  `common.richardson` needs equal refinement ratios (N = 41 … 321: 1.976, 1.988, 1.994).
- **V3**: the notebook implements the method of 2026-09 (first order, unlinked, no convective
  correction). With those options the force history now differs from it by 2.5e-3 (2e-5 on
  2026-09-08; finding 5): reported, as are the defaults and the Galilean equivalence, until the
  reference is regenerated; the independent moving-body reference is V13.
- **V5**: mesh family k (background 0.01/2^k, two nested levels: D/h = 40, 80, 160); 2D-2 from
  rest to t = 16 on 8 ranks, PCG (GAMG stalls on the linked cut-cell pressure: near 1e-6 with
  the default agglomeration, near 2e-9 with pair agglomeration), statistics over the last 3 s.
  2D-1: CG preconditioned by GAMG with the DICGaussSeidel smoother (DIC-PCG hits maxIter 1000 short of relTol
  0.01 at D/h = 160: 24 s per iteration against 3.5 s; 13 iterations per solve); D/h = 160 on 8
  ranks (5661 SIMPLE iterations at D/h = 80; 20335 at D/h = 160). Pressure difference from a
  quadratic least-squares fit over the fluid cells within 4 h of each point (within 2.5 h: ~10
  cells for 6 coefficients, 11 % low). 2D-2 with 3 PIMPLE outer correctors (finding 6).
- **Linear solvers of the immersed cases** (2026-11-11, `doc/crossProcessorLinksHandout.md`):
  the pressure of linked cells needs the incomplete-Cholesky smoother (DICGaussSeidel), not pair
  agglomeration; GAMG with Gauss-Seidel stalls. Steady immersed cases use CG preconditioned by
  GAMG(DICGaussSeidel) for p and smoothSolver for U and nuTilda, linked or not. A linked group
  split over processors slows every solver (V5 D/h 80, 12 split groups: 46 instead of 6.5
  iterations per solve), so decomposed static cases (V5, V12, NACA) keep each group on one
  processor: `ibmSetGeometry -linkZone` and `preserveFaceZones` in decomposeParDict. All
  immersed cases with a Krylov solver (V5, V12, V13, V15, the cylinderRe1, AMR and NACA
  tutorials) wrap their preconditioners in `linkedGroup`, which is its base where no group is
  split. G8 and the moving-body V13 are decomposed without the constraint and rely on it
  (2026-11-11, README §4), with which split groups no longer cost iterations: V13 D/h = 40,
  DIC-PCG 206 to 207 iterations per solve with 1 to 10 split links (204 with none; 276 to 535
  and the limit of 1000 without the wrapper), PBiCGStab-DILU for U 4.3 to 4.9 (8.5 to 57; 3.1
  with none). V13 uses GAMG (coarsest level to relTol 0.01) as the base since
  2026-11-12: DIC-PCG ran into maxIter 1000 at D/h = 160, split groups or not. Its c_d and c_i
  move by up to 6e-5 with the base: at these tolerances the pressure solver sets the fifth
  digit (D/h = 40: DIC 2.14122, GAMG 2.14117, with tolerance 1e-12 and relTol 1e-4 2.14118 and
  2.14117).
- **V12, V13**: nested static refinement (h0 = 3.2, 7–9 levels), 8 ranks; V13 at Courant ≤ 0.8.
- **V15**: the IBM plate gives no C_f(x) yet (no wall-shear output): C_D against body-fitted;
  C_f(0.97) of the body-fitted run against TMR reported. Leading edge half a cell (5e-4)
  downstream of the mesh line x = 0.

### Findings

1. ~~G8: the linked wall-resolved SA channel depends on the decomposition~~ **resolved
   (5008932): not the method, the iteration.** With ν̃ relaxed by 0.9 SIMPLEC ends in an exact
   two-iteration cycle of the driving gradient, body-fitted (6.5e-5) and immersed (2e-4 to
   6.9e-4) alike; where it stops, and in parallel the cycle itself (the Gauss–Seidel smoothers
   act per processor), set the "result". With ν̃ 0.7 (U 0.9) both converge to round-off (Ux, ν̃
   residuals 1e-15), serial = 4 = 6 ranks to 12 digits (G = 0.00296311949746), and the
   immersed walls (f = 0.3) match the body-fitted channel to 1.8e-5 (0.16 % closer than the cycle
   values). `tutorials/turbulence/channelSA/makeCase.py` now writes ν̃ 0.7; its body-fitted
   cases had also been given the `cutCellCorrected` schemes since de164be (af25a0b).
2. ~~G9: pressure at the junction of an immersed wall with the outer boundary~~ **resolved
   (6619272 level wall, e7d742a 45° wall).** Did not
   converge (max |p − p_exact| ≈ 0.07 at N = 20…160 where the level wall meets the outflow
   boundary; body-fitted corner 2e-3 and converging). Relevant to V14 (cylinder spanning the
   channel walls) and V15 (plate meeting the outflow).
   **Level wall: the gate's boundary data.** A cut outer face carries theta U_b.S; G9 took U_b
   at the geometric face centre, which at f = 0.5 lies on the wall (U_b = 0: the wet half of the
   face carried no flux) and at f = 0.1 behind it: an O(1) mass defect of the corner cell at
   every N. With U_b at the centroid of the wet part the junction max converges (f = 0.5: 2.1e-2,
   8.0e-3, 4.0e-3, 2.3e-3; f = 0.1: 3.8e-2 ... 3.8e-3; was 0.07). A fixedValue on a cut face is
   the value over its wet part (README). **Open: the 45° wall** where it meets the top (inflow)
   boundary in an acute (45°) fluid wedge: one small cut cell (alpha 0.19 at N = 160) at
   -1.5e-2, -1.8e-2, -2.6e-2 (N = 40, 80, 160), the same without linking and without
   ibmConvection; with the exact p prescribed on the outer faces too (instead of zeroGradient)
   -1.1e-2 at every N. Its neighbours and the obtuse (135°) junction at the bottom converge
   (1.6e-3); the L2 norms and force orders are unaffected. Not V14/V15 (right-angle junctions).
   **45° wall: the boundary gradient of cut cells (e7d742a).** A fixedValue patch differences
   v_b − v_P over the mesh delta (cell centre to face centre), but v_b is the value at the
   wet-face centroid and v_P at the value point (fluid or group centroid): where the offset has a
   component normal to the boundary the viscous flux of the face is wrong by an O(1) factor
   (right-angle junctions: the offset is parallel to the boundary, no error). `cutCellCorrected`
   now corrects fixed-value faces of non-coupled patches like the internal ones, with
   d = x'_f − x_v. Junction cell: u error 6.6e-3, 1.9e-3, 3.9e-4 (N = 40, 80, 160; was 3.2e-2,
   1.7e-2, 7.8e-3), p error 1.8e-2, 1.2e-2, 5.4e-3 (first order like every zeroGradient-p
   boundary cell, e.g. the top-left corner 8.0e-3, 4.0e-3, 1.9e-3); with the exact p on the
   outer faces the junction is no longer the maximum (1.4e-2 ... 1.7e-3, top row). G1–G8 and
   every tutorial force bit-identical.
3. ~~G9: wall at 45°: interior p max norm stalls near 6e-3~~ **resolved (e2a9780): it is finding
   2.** The 45° wall meets the bottom of the Hiemenz box at x = -0.465, 0.135 from the left face;
   the "interior" excluded the faces other than the bottom, so the junction counted. The max
   error at N = 80, 160 sits within 0.007 of the bottom. Interior = 0.1 from every outer face but
   a body-fitted wall: p Linf 2.8e-2, 6.9e-3, 2.1e-3, 5.1e-4 (order 1.88, was 0.06).
4. ~~G1: the wall area depends on nSubSamples~~ **resolved (f562ca6): cut cells closed by the
   sampling.** The wall area comes from the exact face fractions (closure), not from the samples;
   but a cut cell whose fluid is smaller than one sample got alpha = 0 and was closed like a solid
   cell, which moved its wall segment onto the neighbours' faces and dropped its volume. Not only
   nSubSamples 8: the default 32 left 2-11 such cells per G1 mesh, 64 still 2-3. A candidate cell
   with a wet face now keeps half a sample's fraction (centroid: wet-area-weighted wet face
   centres). 2-D area independent of nSubSamples 8...64, second order (8.93e-4, 2.20e-4, 5.45e-5,
   1.41e-5). Static tutorials bit-identical; moving bodies change (c8c771f, G7: 2.0e-5
   oscillating, 5.8e-4 translating; G6 translating budget 6.5e-12 -> 3.7e-10 per step, within
   its 1e-9; not traced, presumably the cells with half a sample's fraction). No nSubSamples minimum needed
   for the area; the volume error still scales with h/nSubSamples. The inverted (concave)
   cylinder still varied by ~1e-5 of the area: f562ca6 looped over the wrong faces (a `point c`
   shadowed the cell; c8c771f), so for it no missed cell was caught. Since c8c771f the 2-D wall
   area is identical to 12 digits for nSubSamples 8...64, outside and inside, N = 20...160.
5. ~~V3: the first-order moving path no longer reproduces the notebook~~ **resolved (41e9947).**
   Bisect over 11 builds with fixed case files and every option pinned: the force difference
   (2.5e-3) is the same back to 38e21af, the commit that matched: V3 took the maximum over sampled
   steps including the impulsive first two (the first step differs by 0.19 % since 2026-09-08,
   README). The final fields changed at 7f50469: fillUncovered (intended; p of uncovered cells
   from the fluid neighbours). reference_moving.py --fillUncovered ports it (without it the stored
   references are reproduced bit for bit): force 1.8e-5 / 3.4e-5 after the first 10 steps, u
   2.5e-7 / 5.7e-7 (oscillating / translating); V3 checks 1e-4 and 1e-5. Galilean equivalence
   (translating vs fixed cylinder, N = 81): +0.55 % at first order, +0.08 % at the defaults.
6. **One PISO pass per step is first order in time** (V5 2D-2, V12) **resolved**. DFG 2D-2 at D/h = 40
   (nOuterCorrectors 1, Courant 0.8, backward): C_L amplitude +12.8 %, St −1.2 %; at dt/2
   +6.7 %, −0.7 %; with 3 outer correctors (same dt) +0.8 %, −0.3 % (the remaining spatial
   error). D/h = 80 with one pass (+5.9 %) carries the error of its time step, not of its mesh.
   The convecting flux of a single pass is that of the previous step. **Generic PISO, not the
   immersed wall** (D2, `bodyfitted.py dfg`, O-grid M = 64, 44k cells, same Δt): one pass
   C_Lmax +12.9 %, St −1.28 %; dt/2 +7.0 %, −0.82 %; 3 outer correctors +1.1 %, −0.34 % (C_Dmax
   +1.04 / +0.53 / +0.02 %). The one-pass error halves with Δt and matches the immersed one.
   V5 2D-2 and V12 use 3 outer correctors (V13: 2), so does the one time-resolved tutorial
   (`amr/cylinderRe100`); the others are Re = 1 marched to a steady state with Δt ≈ 1. With them V5
   passes at D/h = 80 (St −0.04 %, C_Dmax +0.14 %, C_Lmax −0.07 %, Δp −0.35 %).
7. ~~V12 α = 1: mean C_L −2.482 (D/h = 80) against −2.265 (figure 11 digitised)~~ **resolved
   (668471e): the reference.**, moving away
   on refinement (−2.452 at D/h = 40); α = 2 (steady) converges to its reference (−5.600,
   −5.469 → −5.351), α = 0 passes. **Not finding 6**: with 3 outer correctors (D/h = 80, to
   t = 150) the mean over the same 9 whole cycles is −2.464 (one pass −2.482); C_D, C_L r.m.s.
   and St agree at the same α. Remaining candidates: the reference (read off a bitmap curve; no
   independent tabulated value at Re = 200 found), the rotating-wall treatment at moderate α,
   domain and boundary conditions. Proposed: a body-fitted rotating cylinder at α = 1.
   Body-fitted pimpleFoam on O-grids, same domain and boundary conditions (`bodyfitted.py rot`,
   M = 64 / 128 cells per quarter circle, 61k / 243k cells): α = 1 mean C_L −2.4984 / −2.4994,
   C_D 1.083, C_Lrms 0.4946, St 0.1931; the IBM (3 outer correctors) converges to it: −2.375,
   −2.465 at D/h = 40, 80 (1.4 % short). Figure 11 is 10 % off both. α = 0: body-fitted C_D
   1.3206, St 0.1949, C_Lrms 0.4732 (IBM D/h = 80: 1.3229, 0.1949, 0.4734); α = 2: C_L −5.375,
   C_D 0.2887 (IBM D/h = 40, 80, 160: −5.073, −5.303, −5.358; C_D 0.354, 0.294, 0.292). V12
   checks α = 1 mean C_L against the body-fitted M 128 (`data/rotatingCylinder_bodyFitted`)
   and reports figure 11 and the body-fitted comparison for every α. Not checked: whether the
   shared domain/boundary setup or the digitisation separates both from Mittal & Kumar.
8. ~~V12 statistics over a window that is not a whole number of cycles~~ **resolved**: the last
   50 time units are 9.6 periods at α = 1; the partial period biased mean C_L by −0.021 (α = 1)
   and +0.015 (α = 0, reference 0, tolerance 0.02). `meshes.periodic_statistics` now averages
   over the whole C_L cycles inside the window; V12 re-evaluated, verdicts unchanged.
9. ~~G9: wall shear force order 0.91 for the tilted Homann wall (0.023 → 0.043 at N = 80 with
   the boundary data of 6619272)~~ **resolved: the metrics.** (a) The exact Homann force came from
   a staircase quadrature over the box, 0.85 % off (viscous) and 0.8 % (pressure); now the wall
   polygon, exact. (b) Most of the error was the wall-normal part ν u_n/d of the discrete wall
   flux: exact normal viscous stress on a no-slip wall is zero, the discrete one is first order
   (Homann −0.146, −0.076, −0.039 of |F_v| at N = 20, 40, 80); body-fitted the same (OpenFOAM's
   wall patch, level Hiemenz −0.037 … −0.0042; tangential error IBM / body-fitted 0.0323 /
   0.0305, 0.0088 / 0.0083, 0.0023 / 0.0023, 0.0006 / 0.0006). The old boundary data gave a
   positive normal error that partly cancelled it. G9 checks the tangential part (orders 1.79 to
   1.93; threshold raised from 0.9 to 1.5, 2026-11-10) and reports the normal part. (c) The pressure force depends on the pressure constant
   (fixed in one reference cell): now with the mean error removed as in the norms (finest
   errors 1e-5 … 2e-3; Homann was 4 %).
10. **The wall traction is ν ∂u/∂n, not the stress** (C1, 2026-11-10). The no-slip traction
   `ν (u_P − u_b)/d` is the wall flux of the Laplacian form of the momentum equation, which is
   right for the momentum (div (∇u)ᵀ = 0) and for the force on a closed rigid body, but the moment
   needs the stress `ν (∇u + (∇u)ᵀ)·n`. On the no-slip wall of a rigid body `(∇u)ᵀ·n = n × ω`;
   its moment on a rotating circular cylinder is `−2πνωR²` per unit length. Taylor-Couette
   (V1, inner cylinder rotating): torque −1.963 against the exact −3.534 (−45 %, at every N)
   without it, outer cylinder (at rest) second order already; with it both second order:
   inner 7.4e-3, 1.8e-3, 4.7e-4, 1.1e-4, outer 5.4e-3, 1.4e-3, 3.7e-4, 7.9e-5 (N = 40…320,
   linked); secondOrder false 1.8e-2, 7.7e-3 (inner) at N = 80, 160. `forces()` adds
   `ν ω × S_k` (zero sum on a closed body) and its moment `ω tr T − ω·T` exactly from the wall
   pieces. V12 body-fitted C_M (OpenFOAM forces, full stress): 0.2296 / 0.2295 (α = 1),
   0.4635 / 0.4633 (α = 2) at M = 64 / 128; figure 11 ((C_P − C_D)/α) 0.238 and 0.468.
11. ~~V12 α = 2, D/h = 40 depends on the decomposition (5d9caa5)~~ **resolved (8d416fd)**: the
   two copies of a processor-face flux drift apart, amplified by ddtCorr without the fluid
   fraction; ddtCorr with α⁰ and a single-valued coupled-face flux (§12 amendments, "V5, V12,
   V13 re-accepted from 8d416fd").
12. ~~2-D transient runs with a linked free-slip wall: U_z grows without bound~~ **resolved
   (361c23c)**. `fvMatrix::H()` sets the components in empty directions to zero; `mergeH`
   added the link-face terms back for every component, so the merged H_z of a group became its
   internal coupling. In a group of a small cut cell without a wall coefficient (free slip;
   α = 0.004 slave under a full master) H/A exceeded U_z by 1.28 per corrector sweep: round-off
   of grad(p)_z (1e-32) grew to U_z = 1e31 in 6 time units (Taylor-Couette, N = 160,
   ibmPimpleFoam, old and zero-shear free slip alike); it crashed a concave free-slip case at
   t = 22. No-slip walls hid it (A_w/(d V) in the merged A); the in-plane solution never
   depended on U_z. Found by elimination (convection, ddtCorr, ibmConvection, consistent,
   friction lag ruled out by fresh runs) and by printing H, A, the face coefficients of the
   growing group (`~/runs/zeroshear/root`). G7: 11 out-of-plane round-off metrics (Fz, Mx, My
   ≤ 1e-18) changed, everything else bit-identical.
13. **Open: unlinked free slip with alphaMin > 0** (2026-11-14). Absorbing slivers (α < 0.02)
   into the body closes faces across the slip flow: steps that drag on it. Concave
   Taylor-Couette with a zero-shear outer wall, unlinked: u L2 0.082 (N = 40), 0.059 (N = 80)
   with alphaMin 0.02, 4.3e-4 with alphaMin 1e-9; V2 zero shear unlinked 3.9747, 3.9936, 4.0032,
   4.0119 (N = 41 … 321), rising past the reference; with alphaMin 1e-9 3.9748, 3.9936, 3.9980,
   3.9991, as linked. The old free slip has the same defect.
   Linked runs (alphaMin 0) are not affected.
14. **V14 3D-1Z: the duct upstream sets C_D, the cut pattern sets C_L** (2026-11-15, 2190505,
   runs `~/runs/v14`, `~/runs/suite/2190505/V14`). With the boxes of the spec C_D = 6.1214,
   6.1625, 6.1705 at D/h 15.6, 23.4, 31.2 (order 3.6, extrapolated −0.17 %): 16–32 cells across
   the duct between the inlet and the boxes. Larger boxes symmetric about the channel centre:
   6.2204, 6.2000, 6.1928, order 1.96, extrapolated 6.18323 (−0.034 %; FeatFlow level 4
   6.18465). The junction of the cylinder with the walls z = 0, H needs nothing special: G1
   closure 0 there, second-order C_D. C_L is the same with both boxes (+142, −40, −19 %;
   −6.8 % at D/h 46.8): the lift, 0.15 % of the drag, follows the cut pattern; DFG 2D-1 on the
   same y-mesh (no junction) gives +114, −57, +29, +2.3, −4.8 % at D/h 15.6 … 62.4. Mesh,
   D/h 46.8 and the C_L tolerance of 10 % there decided by the user.
15. **V14 3D-3Z: the time-step dependence is the Rhie–Chow rAU, a resolution effect** (2026-11-16,
   2190505, runs `~/runs/v14/pimple3d`, `~/runs/v14/slicepimple`). V14.1 (Courant 0.60): c_D,max
   3.4512 (+4.65 %), 3.3760 (+2.37 %) at D/h 15.6, 23.4 (order 1.66), Δt/2 at 15.6 −0.47 %; c_L,min
   −0.0069, −0.0164 (−0.0110): the lift unresolved, as in 3D-1Z. ibmPimpleFoam run to the steady
   state differs from ibmSimpleFoam (SIMPLEC) by the rAU in the Rhie–Chow flux (Δt against the
   relaxation floor); the ddtCorr, the junction and the hanging nodes are excluded (ddtCorr off:
   −0.81 %; uniform mesh: −0.53 %). Re 44, PIMPLE − SIMPLE c_D at Co 0.6 / 0.3: 3-D −0.51 / −0.73 %
   (D/h 15.6), −0.17 % (23.4); 2-D slice −0.56 / −0.81, −0.055 / −0.091, −0.002 / −0.006 % at
   D/h 15.6, 31.2, 62.4: the same in 2-D and 3-D, order ≥ 3 in h, both solvers share the limit.
   Δt from |u| = 3.75; c_D,max and c_L,min judged at D/h 46.8 (user).
16. **Cost of the immersed wall against body-fitted runs** (2026-11-16, `timing.py`, report
   §Cost). Host and load dominate the stored times: the same work took 4.5–6.5× longer in the
   tier-1 run on eridanus than on io2/hydra (V4 vs V2). Same mesh and solver (G9): 0.8–1.5× per
   iteration, 11–25 % fewer SIMPLE iterations, 0.7–1.3× in total; V15 same grid 1.04× per
   iteration (2026-11-11: 3.9×, then 1.34× with the cached gradient); V7 (100 cells) 1.45×.
   Same host (cetus, v2412, continued from the final states): V12 α = 1 immersed D/h 80 21–22 µs
   per cell and step against body-fitted M 64 38 µs (829 against 1838 pressure iterations: the
   stretched O-grid), V5 2D-2 1.0×. Equal accuracy is what differs: steady cylinders (V2, V4)
   2–2.5× the CPU time of the O-grid; V12 mean C_L −1.20 % (α = 1, D/h 80, 27 core-h) and
   −0.30 % (α = 2, D/h 160, 286 core-h) against −0.04 % and −0.007 % of body-fitted M 64
   (35, 38 core-h). The steady immersed cases G9, V2, V4 use PCG/DIC, whose iterations double
   with N (V2: 21 … 208 per step), where V15 with a GAMG base needs 2: a case setting. Decision (user, 2026-11-16): the GAMG
   base stays as it is (V5, V13, V15), G9/V2/V4 and their baselines are not changed; the
   immersed cost is to come down with AMR (degrees of freedom only where needed).

17. **V11 BUB: the zero-shear free slip gives the clean-bubble drag in 3-D** (2026-11-17, 2190505,
   `~/runs/suite/2190505/V11.1`, report §V11). C_D at D/h 24, 32, 48: Re 50 0.6825, 0.6774,
   0.6748; Re 100 0.3978, 0.3858, 0.3782; Re 200 0.2403, 0.2229, 0.2102; Re 300 0.1841, 0.1645,
   0.1498; orders 3.1, 2.4, 1.9, 1.9; extrapolated +1.5, +1.2, −0.2, −1.0 % from Blanco &
   Magnaudet (PASS, user decision above). Max tangential stress (μU/a) at D/h 48: 0.095, 0.21,
   0.43, 0.64, order 0.9–1.5; surface vorticity/(2u_θ/a) at the equator 0.95 … 0.80; no reversed
   flow on the axis. Lateral 15 D −0.02 % (Re 50), ibmConvection off −0.07 %, zeroShear off
   C_D 0.216 (−44 %, vorticity ratio 0.41: the former free slip). mergeAlpha 0 (unlinked,
   alphaMin 0.02, 264 slivers absorbed) cycles with period 2 between 0.3787 and 0.3868 (mean
   −0.8 %) through 20000 iterations: finding 13 in 3-D, reported, not judged. Baseline accepted from
   2190505 (user, 2026-11-17).

### Decisions (2026-11-07)

- **D1** Tier-1 baseline accepted from 5008932 (`validation/suite/baseline/`), V12 with its FAIL;
  V12 is re-accepted after D3/D4.
  Done (2026-11-10): baseline re-accepted from e7d742a (G1–G9, V1–V5, V7, V13, V15; V1–V5, V7,
  V13, V15 bit-identical to c8c771f; G9 with the metrics of finding 9) and V12 from c8c771f
  (no cut cell touches an outer boundary, e7d742a cannot change it): every V12 check against
  the literature and the body-fitted references passes, α = 2 at D/h = 160 C_L −5.357 (the same
  case on 128 instead of 48 ranks agrees to 6 digits). G7 takes the moving-body change of
  finding 4.
- **D2** (finding 6) Body-fitted DFG 2D-2 with one pass at D/h = 40 first, to separate generic
  PISO from lagged wall terms; then 2–3 outer correctors in the transient tutorials and docs.
  Done: generic PISO (finding 6); `amr/cylinderRe100` to 3 outer correctors, README.
- **D3** (finding 7) Body-fitted rotating cylinder at α = 1 (O-grid, same domain, BCs and
  rotating-wall velocity) as independent reference; the check stays FAIL until then.
  Done (668471e): the IBM agrees with it, figure 11 does not (finding 7).
- **D4** V12 α = 2: add D/h = 160 (steady, locally).
  Done: −5.073, −5.303, −5.358 at D/h = 40, 80, 160 (body-fitted −5.375): converging.
- **D5** (finding 2) Junction pressure stays REPORT in G9; fixed before V14–V16 (tier 2).
- **D6** (finding 3) 45° wall: locate the cells of the max error and isolate the cut topology now.
- **D7** (finding 4) Minimum nSubSamples (16 in 2-D) in tutorials and docs with a warning below
  it now; wall area independent of the sub-sampling later.
- **D8** (finding 5) Bisect V3 over bc08f21, b01303a, 7f50469; regenerate the notebook reference
  if the change is intended.
- **D9** Merge feature_amr right after D1; rerun tier 1 on the merge against the baseline; the
  D2–D8 work continues on the merged branch.
