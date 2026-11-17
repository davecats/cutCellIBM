    // Momentum predictor on the cut mesh:
    //   theta-weighted convection (phi carries the wet face area),
    //   theta-weighted laplacian (no wall term of its own; the stress-free
    //   wall flux of free-slip walls, zeroShear, added below),
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

    // Convection (cutCellGeometry::convection); turbulent: the wall law in the
    // wall cells with a wall function on the faces with ibmConvection
    tmp<volScalarField> tnuLam;
    if (ibmTurbulent) tnuLam = turbulence->nu();
    const volScalarField* nuLam = ibmTurbulent ? &tnuLam() : nullptr;

    tmp<fvVectorMatrix> tUEqn
    (
        ibm.convection(phi, U, nuLam)
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
                ibm.wallVelocity()
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

    if (ibm.wallVelocity())
    {
        // Invariant bodies (static geometry, wall velocity u_b): the u_b
        // half of the no-slip wall traction nu (u_b - u_P) Awall/d, and the
        // blanking source driving the solid interior to u_b
        UEqn -=
            ibmTurbulent
          ? ibm.movingWallSource(nuWall, nuEff)
          : ibm.movingWallSource(nuEff);
    }

    if (ibm.zeroShear())
    {
        // Stress-free free-slip walls: the wall flux of the Laplacian form,
        // nu (H.(u - u_b)_t + omega x n + (du_n/dn) n) (zeroShear false:
        // du/dn = 0, doc/freeSlipPlan.md)
        UEqn += ibm.zeroShearSource(nuEff, U);
    }

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        if (ibm.linking())
        {
            // Small cut cells linked to their masters for the solution only:
            // pEqn takes the coefficients of the merged control volumes from
            // the unlinked matrix (ibm.mergeA, ibm.mergeH)
            // (the residual redistributed in the groups, which the
            // pressure corrector makes uniform: see cutCellGeometry::linkCells)
            fvVectorMatrix UEqnLinked(UEqn == -ibm.grad(p));
            ibm.linkCells(UEqnLinked, true);
            solve(UEqnLinked);
        }
        else
        {
            solve(UEqn == -ibm.grad(p));
        }

        fvOptions.correct(U);
    }
