# Spalart-Allmaras channel: body-fitted vs immersed walls

`results.txt`: output of `validation/channelSA_compare.py` for the cases of
`tutorials/turbulence/channelSA/makeCase.py` (the case directories were not kept). Body-fitted
runs with `simpleFoam` + `SpalartAllmaras` + `nutUSpaldingWallFunction`; IBM runs with
`ibmSimpleFoam` + `ibmSpalartAllmaras` + `wallDist ibmLevelSet`, walls as two `plane` bodies cut at
the fluid fraction f (bottom) and 1 - f (top) of their cell.
