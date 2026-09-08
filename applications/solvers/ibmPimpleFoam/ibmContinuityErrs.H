// Continuity error of the cut-cell flux: div(phi) + S must vanish, S being
// the wall flux of a moving body (zero for a body at rest)
{
    volScalarField contErr(fvc::div(phi));
    if (ibm.moving())
    {
        contErr += ibm.massSource();
    }

    scalar sumLocalContErr = runTime.deltaTValue()*
        mag(contErr)().weightedAverage(mesh.V()).value();

    scalar globalContErr = runTime.deltaTValue()*
        contErr.weightedAverage(mesh.V()).value();
    cumulativeContErr += globalContErr;

    Info<< "time step continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr
        << endl;
}
