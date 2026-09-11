/*---------------------------------------------------------------------------*\
Application
    ibmPimpleFoam

Description
    Transient solver for incompressible, laminar flow around immersed bodies
    represented by the cut-cell immersed boundary method of cutCellGeometry.
    PIMPLE (merged PISO-SIMPLE) pressure-velocity coupling as in pimpleFoam,
    on a static mesh.

    Bodies may move (doc/moving_ibm.pdf): the mesh is fixed and the
    fractions are rebuilt every time step; continuity acquires the wall
    flux as a source, the unsteady term uses the swept old fluid volume, and
    the no-slip wall and the blanking are driven to the body velocity. The
    static path is untouched: every moving-body operation sits behind
    ibm.moving().

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
#include "ibmScalarTransportList.H"

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

        if (ibm.moving())
        {
            // Move the bodies and rebuild theta, alpha, Sw, the mass source
            // and the swept old fluid volume (alpha.oldTime())
            ibm.update();

            // Put the carried flux back on the geometry the step is about
            // to use: a face that has just opened has no history, one that
            // has just closed carries a history that no longer exists. For
            // a fluid moving rigidly with the body this lands exactly on
            // the fixed point (free stream preserved from the first outer
            // iteration). phi.oldTime() keeps the previous Rhie-Chow flux.
            phi = ibm.flux(U);

            pRefCell = ibm.pressureRefCell(pRefCellUser);
        }

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

        // Passive scalars, if constant/scalarTransportProperties exists
        scalars.solve();

        ibm.reportForces(p, U, turbulence->nuEff());

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}
