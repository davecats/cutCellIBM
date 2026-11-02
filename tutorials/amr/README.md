# Adaptive mesh refinement tutorials

`cylinderRe100`: laminar flow past a cylinder at Re = 100 (D = 1, U = 1). The background
mesh is uniform, 40 × 24 cells of h = D/2. At the first time step `ibmPimpleFoam` refines it
to h = D/32 at the cylinder (`maxLevel 4`, two finest cell layers on either side of the
wall) and to D/16 within half a diameter of it. Every 10 steps it then refines the wake
where |ω| > 0.5 (the `vorticity` function object) to D/16 and coarsens where the wake has
moved on.

```
./Allrun
```

The cylinder sits 0.01 D above the centre line, so the shedding starts within about
40 convective times. `refinementLevel` in every time directory shows the mesh levels, and
`postProcessing/ibmForces` the forces (per unit depth 0.5).

Settings, criteria and limitations: `doc/amr.md`.
