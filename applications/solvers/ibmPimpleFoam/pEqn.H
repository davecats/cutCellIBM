// A and H of the merged control volumes for linked small cut cells
tmp<volScalarField> tA(UEqn.A());
tmp<volVectorField> tH(UEqn.H());
// Gradient increments of ibmMeanVelocityForce in this corrector
fv::ibmMeanVelocityForce::addCorrectorSources(fvOptions, tH.ref());
ibm.mergeA(UEqn, tA.ref());
ibm.mergeH(UEqn, tH.ref());
volScalarField rAU(1.0/tA);
volVectorField HbyA(constrainHbyA(rAU*tH, U, p));
surfaceScalarField phiHbyA("phiHbyA", ibm.flux(HbyA));

if (pimple.ddtCorr())
{
    phiHbyA += ibm.ddtCorr(rAU, U, phi);
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
    tmp<volScalarField> H1(UEqn.H1());
    ibm.mergeA(UEqn, H1.ref());
    rAtU = 1.0/max(1.0/rAU - H1, ibm.alpha()/runTime.deltaT());
    phiHbyA +=
        ibm.thetaInterpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    volVectorField gradp(ibm.grad(p));
    ibm.groupMean(gradp);
    HbyA -= (rAU - rAtU())*gradp;
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
    ibm.linkCells(pEqn);
    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
        ibm.syncCoupledFlux(phi);
    }
}

#include "ibmContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

// Pressure force over the cell volume, merged for linked small cut cells
volVectorField gradp(ibm.grad(p));
ibm.groupMean(gradp);
U = HbyA - rAtU*gradp;
U.correctBoundaryConditions();
fvOptions.correct(U);
