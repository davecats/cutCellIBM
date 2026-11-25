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
    ibm.moving(). A body whose motion leaves it invariant (rotation about
    its axis, sliding along it) keeps a static geometry and only supplies
    the wall velocity (ibm.wallVelocity()).

    Free bodies (motion type free, see freeMotion) move under the force and
    moment of the flow: their state is predicted at the start of every time
    step from the loads at the old time level and, with implicit coupling,
    corrected at the start of every further PIMPLE outer corrector from the
    loads of the flow just solved, the geometry rebuilt each time
    (update(false)); nOuterCorrectors sets the number of coupling
    iterations. With couplingTolerance (PIMPLE dictionary, default 0: never)
    the body state is kept for the rest of the step once its change between
    correctors (relative to its velocity scale) is below it: the remaining
    outer correctors converge the flow without rebuilding the geometry. The
    outer loop of a viscous step contracts slowly (~0.9 per corrector at
    nu dt/h^2 ~ 4), the coupling iterate much faster, and the rebuild is the
    largest part of a corrector.

    The transient Rhie-Chow correction (ddtCorr) is available in its Euler
    form, theta- and alpha-consistent; the flux of a coupled face is made
    single-valued after the pressure correction. Turbulence models are constructed for their
    effective viscosity only; simulationType laminar is supported.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "cutCellGeometry.H"
#include "ibmMeanVelocityForce.H"
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
    #include "createDynamicFvMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // Loads on the free bodies, with the viscosities of reportForces
    auto freeBodies = [&](const word& stage) -> scalar
    {
        const tmp<volScalarField> tnuMolecular(turbulence->nu());
        const tmp<volScalarField> tnu
        (
            ibmTurbulent
          ? ibm.wallViscosity(U, turbulence->nu())
          : turbulence->nuEff()
        );
        if (stage == "predict")
        {
            ibm.predictFreeBodies(p, U, tnu(), &tnuMolecular());
            return 0;
        }
        else if (stage == "correct")
        {
            return ibm.correctFreeBodies(p, U, tnu(), &tnuMolecular());
        }
        ibm.reportFreeBodyLoads(p, U, tnu(), &tnuMolecular());
        return 0;
    };

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

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
            setRefCell(p, pimple.dict(), pRefCellUser, pRefValue);
            pRefCell =
                ibm.moving()
              ? ibm.pressureRefCell(pRefCellUser)
              : pRefCellUser;
        }

        // Free bodies: the state at the new time from the loads at the old
        if (ibm.freeBodies())
        {
            freeBodies("predict");
        }

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

            // Cells just uncovered carry the solid pressure (zero), a jump
            // of the size of the pressure level for the momentum predictor
            ibm.fillUncovered(p);
        }
        else
        {
            // Invariant bodies (motion with "invariant yes"): the wall
            // velocity at the new time, static geometry. No-op without them.
            ibm.update();
        }

        // Coupling iterations of this step: until the change of the body
        // state is below couplingTolerance
        const scalar couplingTol =
            pimple.dict().getOrDefault<scalar>("couplingTolerance", 0);
        bool coupling = true;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            // Implicit coupling: correct the free bodies from the flow of
            // the previous outer corrector and rebuild the geometry within
            // the step (the old geometry of the step is kept)
            if (ibm.implicitFreeBodies() && !pimple.firstIter() && coupling)
            {
                const scalar change = freeBodies("correct");
                Info<< "Free bodies: change of the coupling iterate "
                    << change << endl;
                if (change < couplingTol)
                {
                    // The geometry follows this last correction; the state
                    // is kept for the rest of the step
                    coupling = false;
                    Info<< "Free bodies: coupling converged (tolerance "
                        << couplingTol << "), state kept for the rest of "
                        << "the step" << endl;
                }
                ibm.update(false);
                phi = ibm.flux(U);
                pRefCell = ibm.pressureRefCell(pRefCellUser);
                ibm.fillUncovered(p);
            }

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

        {
            // molecular viscosity: transpose stress on rotating walls
            const tmp<volScalarField> tnuMolecular(turbulence->nu());
            ibm.reportForces
            (
                p,
                U,
                ibmTurbulent
              ? ibm.wallViscosity(U, turbulence->nu())()
              : turbulence->nuEff()(),
                &tnuMolecular()
            );
        }

        if (ibm.freeBodies())
        {
            freeBodies("final");
        }

        runTime.write();
        ibm.writeFreeBodies();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}
