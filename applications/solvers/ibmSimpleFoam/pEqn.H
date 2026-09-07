{
    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phiHbyA("phiHbyA", ibm.flux(HbyA));
    adjustPhi(phiHbyA, U, p);

    tmp<volScalarField> rAtU(rAU);

    if (simple.consistent())
    {
        // SIMPLEC mobility 1/(A - H1); the floor only guards against a
        // vanishing denominator (simpleFoam has none)
        rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 1e-3/rAU);
        phiHbyA +=
            ibm.thetaInterpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
        HbyA -= (rAU - rAtU())*ibm.grad(p);
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
    U = HbyA - rAtU()*ibm.grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
}
