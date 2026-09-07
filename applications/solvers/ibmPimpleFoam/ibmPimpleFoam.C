/*---------------------------------------------------------------------------*\
Application
    ibmPimpleFoam

Description
    Transient solver for incompressible, laminar flow around immersed bodies
    represented by the cut-cell immersed boundary method of cutCellGeometry.
    PIMPLE (merged PISO-SIMPLE) pressure-velocity coupling as in pimpleFoam,
    on a static mesh.

    The transient Rhie-Chow correction (ddtCorr) is available in its Euler
    form, theta-consistent. Turbulence models are constructed for their
    effective viscosity only; simulationType laminar is supported.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "cutCellGeometry.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible laminar flow with cut-cell"
        " immersed boundaries (PIMPLE)."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        ibm.reportForces(p, U, turbulence->nuEff());

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}
