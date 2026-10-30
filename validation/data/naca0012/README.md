# NACA 0012, Re = 1e6, Spalart-Allmaras: body-fitted vs immersed

`results.txt`: output of `validation/naca0012_compare.py` (lift, drag, convergence) and
`validation/naca0012_surface.py` (C_p, C_f; cell values and `--image 3e-3`) for cases of
`tutorials/turbulence/naca0012/makeCase.py` (the case directories were not kept):

    makeCase.py tmr_wr_a0 --kind bf --mesh wr --grid n0012_449-129.p2dfmt --endTime 10000
    makeCase.py tmr_wf_a0 --kind bf --mesh wf --grid n0012_449-129.p2dfmt --endTime 10000
    makeCase.py ibm_a0    --kind ibm --R 500 --mergeAlpha 0.5 --endTime 10000 [--secondOrder]

and the same with `--alpha 10`. The TMR grid n0012_449-129.p2dfmt is in
https://www.nasa.gov/wp-content/uploads/2026/02/naca0012-grids.zip (see
https://tmbwg.github.io/turbmodels/naca0012_grids.html). Body-fitted runs: `simpleFoam` +
`SpalartAllmaras`, `nutUSpaldingWallFunction` on the wall; IBM runs: `ibmSimpleFoam` +
`ibmSpalartAllmaras`, `wallDist ibmLevelSet`, Spalding wall function, cell linking.
