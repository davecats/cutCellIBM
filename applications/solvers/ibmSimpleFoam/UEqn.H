    // Momentum predictor on the cut mesh:
    //   theta-weighted convection (phi carries the wet face area),
    //   theta-weighted laplacian (no wall term: zero tangential traction),
    //   the no-slip wall traction nu Awall/dWall on the diagonal,
    //   and a blanking diagonal in the solid cells.

    tmp<volScalarField> tnuEff(turbulence->nuEff());
    const volScalarField& nuEff = tnuEff();

    // Viscosity of the no-slip wall traction: nuEff (= nu) for laminar flow, the
    // immersed-wall counterpart of a nut wall function for turbulent flow
    tmp<volScalarField> tnuWall;
    if (ibmTurbulent)
    {
        tnuWall = ibm.wallViscosity(U, turbulence->nu());
    }
    const volScalarField& nuWall = ibmTurbulent ? tnuWall() : nuEff;

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      - fvm::laplacian(ibm.thetaInterpolate(nuEff), U)
      + fvm::Sp(nuWall*ibm.noSlipCoeff(), U)
      + fvm::Sp(nuEff*ibm.blankCoeff(), U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    if (ibmTurbulent)
    {
        // Transpose part of the viscous stress, div(nuEff dev2(T(grad U))), on the
        // wet faces, with the velocity gradient over the fluid and the wall
        // velocity on the immersed wall (zero for a constant viscosity)
        const volTensorField gradU
        (
            ibm.fluidGrad
            (
                U,
                ibm.moving()
              ? vectorField(ibm.Ub().primitiveField())
              : vectorField(mesh.nCells(), Zero)
            )
        );
        UEqn -= fvc::div
        (
            ibm.thetaInterpolate(nuEff)
           *(mesh.Sf() & fvc::interpolate(dev2(T(gradU))))
        );
    }

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve(UEqn == -ibm.grad(p));

        fvOptions.correct(U);
    }
