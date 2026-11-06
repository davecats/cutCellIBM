# NASA TMR, Spalart-Allmaras verification data (cases V15 2DZP, V16 2DBC)

Copied unchanged on 2026-11-06 from the TMR mirror https://tmbwg.github.io/turbmodels/ (the
former turbmodels.larc.nasa.gov now redirects to nasa.gov), file name = path with `/` → `_`:

| file | source path | content |
|---|---|---|
| `FlatPlate_SA_cf_convergence.dat` | FlatPlate/SA/cf_convergence.dat | C_f(x = 0.97008) on the 5 grids, CFL3D, FUN3D (M = 0.2) |
| `FlatPlate_SA_drag_convergence.dat` | FlatPlate/SA/drag_convergence.dat | C_D (reference length 2) |
| `FlatPlate_validation_cf_incomp_results_sa.dat` | FlatPlate_validation/cf_incomp_results_sa.dat | C_f(0.97), incompressible FUN3D and SC/Tetra |
| `FlatPlate_SA_cf_plate.dat` | FlatPlate/SA/cf_plate.dat | C_f(x), finest grid 545 x 385 |
| `FlatPlate_SA_flatplate_u+y+.dat`, `FlatPlate_SA_flatplate_u.dat` | FlatPlate/SA/ | profiles at x = 0.97008, 1.90334 (CFL3D) |
| `FlatPlate_SA_mut_0.97.dat` | FlatPlate/SA/mut_0.97.dat | nu_t/nu at x = 0.97 |
| `Bump_SA_force_convergence.dat` | Bump/SA/force_convergence.dat | C_L, C_D, C_Dp, C_Dv (reference length 1.5) |
| `Bump_SA_cf_convergence_sa.dat` | Bump/SA/cf_convergence_sa.dat | C_f at x = 0.6321975, 0.75, 0.8678025 |
| `Bump_SA_cf_bump.dat`, `Bump_SA_cp_bump.dat` | Bump/SA/ | C_f(x), C_p(x), finest grid 1409 x 641 |
| `Bump_SA_bump_u.dat` | Bump/SA/bump_u.dat | u(y - y0) at x = 0.75, 1.20148 (CFL3D) |
| `Bump_SA_mut_0.75.dat` | Bump/SA/mut_0.75.dat | nu_t/nu at x = 0.75 |

Reader and constants: `validation/reference_tmr.py`. Grids (for the body-fitted reference runs):
https://tmbwg.github.io/turbmodels/flatplate_grids.html, bump_grids.html.
