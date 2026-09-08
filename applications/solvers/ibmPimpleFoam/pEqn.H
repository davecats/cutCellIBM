volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA("phiHbyA", ibm.flux(HbyA));

if (pimple.ddtCorr())
{
    phiHbyA += fvc::interpolate(rAU)*ibm.ddtCorr(U, phi);
}

if (p.needReference())
{
    adjustPhi(phiHbyA, U, p);
}

tmp<volScalarField> rAtU(rAU);

if (pimple.consistent())
{
    // SIMPLEC mobility 1/(A - H1), floored at the transient diagonal
    // alpha/deltaT as in the reference implementation
    rAtU = 1.0/max(1.0/rAU - UEqn.H1(), ibm.alpha()/runTime.deltaT());
    phiHbyA +=
        ibm.thetaInterpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*ibm.grad(p);
}

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAtU());

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(ibm.thetaInterpolate(rAtU()), p)
     ==
        fvc::div(phiHbyA)
    );

    if (ibm.moving())
    {
        // The moving wall enters the algorithm here and nowhere else: the
        // corrected flux is made to satisfy div(phi) = -S, S = u_b . Sw
        pEqn -= ibm.massSource();
    }

    ibm.constrainPressure(pEqn);
    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

#include "ibmContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAtU*ibm.grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
