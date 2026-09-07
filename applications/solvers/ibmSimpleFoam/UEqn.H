    // Momentum predictor on the cut mesh:
    //   theta-weighted convection (phi carries the wet face area),
    //   theta-weighted laplacian (no wall term: zero tangential traction),
    //   the no-slip wall traction nu Awall/dWall on the diagonal,
    //   and a blanking diagonal in the solid cells.

    tmp<volScalarField> tnuEff(turbulence->nuEff());
    const volScalarField& nuEff = tnuEff();

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      - fvm::laplacian(ibm.thetaInterpolate(nuEff), U)
      + fvm::Sp(nuEff*ibm.noSlipCoeff(), U)
      + fvm::Sp(nuEff*ibm.blankCoeff(), U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve(UEqn == -ibm.grad(p));

        fvOptions.correct(U);
    }
