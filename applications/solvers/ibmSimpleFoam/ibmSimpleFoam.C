/*---------------------------------------------------------------------------*\
Application
    ibmSimpleFoam

Description
    Steady-state solver for incompressible, laminar flow around immersed
    bodies represented by the cut-cell immersed boundary method of
    cutCellGeometry (fluid area and volume fractions, closure-defined wall
    segment, no-slip or free-slip wall). SIMPLE/SIMPLEC pressure-velocity
    coupling as in simpleFoam.

    Turbulence models are constructed for their effective viscosity but
    their own transport equations are not immersed-boundary aware, so only
    simulationType laminar is supported.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "cutCellGeometry.H"
#include "ibmScalarTransportList.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state solver for incompressible laminar flow with cut-cell"
        " immersed boundaries (SIMPLE/SIMPLEC)."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Adaptive refinement (constant/dynamicMeshDict, ibmRefineFvMesh):
        // the fields and the flux have been mapped, the flux rebuilt on the
        // changed faces and the cut-cell geometry rebuilt. The cells that
        // have become solid are blanked and the reference cell found again.
        mesh.update();
        if (mesh.topoChanging())
        {
            if (!ibm.moving())
            {
                ibm.zeroSolid(U);
                ibm.zeroSolid(p);
            }
            setRefCell(p, simple.dict(), pRefCell, pRefValue);
            ibm.checkRefCell(pRefCell);
        }

        // --- Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "pEqn.H"
        }

        laminarTransport.correct();
        turbulence->correct();

        // Passive scalars, if constant/scalarTransportProperties exists
        scalars.solve();

        ibm.reportForces
        (
            p,
            U,
            ibmTurbulent
          ? ibm.wallViscosity(U, turbulence->nu())()
          : turbulence->nuEff()()
        );

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}
