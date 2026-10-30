{
    // A and H of the merged control volumes for linked small cut cells
    tmp<volScalarField> tA(UEqn.A());
    tmp<volVectorField> tH(UEqn.H());
    ibm.mergeA(UEqn, tA.ref());
    ibm.mergeH(UEqn, tH.ref());
    volScalarField rAU(1.0/tA);
    volVectorField HbyA(constrainHbyA(rAU*tH, U, p));
    surfaceScalarField phiHbyA("phiHbyA", ibm.flux(HbyA));
    adjustPhi(phiHbyA, U, p);

    tmp<volScalarField> rAtU(rAU);

    if (simple.consistent())
    {
        // SIMPLEC mobility 1/(A - H1); the floor only guards against a
        // vanishing denominator (simpleFoam has none)
        tmp<volScalarField> H1(UEqn.H1());
        ibm.mergeA(UEqn, H1.ref());
        rAtU = 1.0/max(1.0/rAU - H1, 1e-3/rAU);
        phiHbyA +=
            ibm.thetaInterpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
        volVectorField gradp(ibm.grad(p));
        ibm.groupMean(gradp);
        HbyA -= (rAU - rAtU())*gradp;
    }

    tUEqn.clear();

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAtU());

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(ibm.thetaInterpolate(rAtU()), p)
         ==
            fvc::div(phiHbyA)
        );

        ibm.constrainPressure(pEqn);
        ibm.linkCells(pEqn);
        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    #include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    // Pressure force over the cell volume, merged for linked small cut cells
    volVectorField gradp(ibm.grad(p));
    ibm.groupMean(gradp);
    U = HbyA - rAtU()*gradp;
    U.correctBoundaryConditions();
    fvOptions.correct(U);
}
