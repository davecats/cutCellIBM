# Validation suite: design (in progress)

Started 2026-11-05. Principles: `doc/developerNotes.md`. Each case states which components it
**stresses** and which it **relies on**; a gate is a cheap check with a tolerance that must pass
before the cases that rely on it are run or interpreted. The inventory of the existing cases, the
component list and the option matrix are still to be written; this file so far holds the cases
already specified.

## Case MIG: inertial migration of a sphere in a pipe (Rivero-Rodriguez & Scheid 2018)

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

Not reproducible with the code: everything with Ca > 0 (§3.2–3.3, figures 11–19, table 5,
eq. C2) needs a deformable interface with surface tension.

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
