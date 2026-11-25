/*---------------------------------------------------------------------------*\
    cutCellIBM: bodies moving freely under the flow (freeMotion), part of
    cutCellGeometry

    The solver calls, in every time step,
        predictFreeBodies()   before update(): the force of the fields at the
                              old time level gives the predicted state;
        correctFreeBodies()   implicit coupling only, at the start of every
                              PIMPLE outer corrector after the first: the
                              force of the flow just solved corrects the
                              state, followed by update(false);
        writeFreeBodies()     after runTime.write().
\*---------------------------------------------------------------------------*/

#include "cutCellGeometry.H"
#include "freeMotion.H"
#include "fvMesh.H"
#include "Time.H"
#include "IFstream.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::freeMotion* Foam::cutCellGeometry::freeMotionOf(const label bodyi)
{
    return dynamic_cast<freeMotion*>(&bodies_[bodyi].motion());
}


const Foam::freeMotion* Foam::cutCellGeometry::freeMotionOf
(
    const label bodyi
) const
{
    return dynamic_cast<const freeMotion*>(&bodies_[bodyi].motion());
}


bool Foam::cutCellGeometry::freeBodies() const
{
    forAll(bodies_, bodyi)
    {
        if (freeMotionOf(bodyi)) return true;
    }
    return false;
}


bool Foam::cutCellGeometry::implicitFreeBodies() const
{
    forAll(bodies_, bodyi)
    {
        const freeMotion* fm = freeMotionOf(bodyi);
        if (fm && fm->implicit()) return true;
    }
    return false;
}


void Foam::cutCellGeometry::initFreeBodies()
{
    if (!freeBodies())
    {
        return;
    }

    // State of a restart
    const fileName stateFile
    (
        mesh_.time().timePath()/"uniform"/"ibmFreeBodies"
    );
    autoPtr<dictionary> statePtr;
    {
        IFstream is(stateFile);
        if (is.good())
        {
            statePtr.reset(new dictionary(is));
        }
    }

    forAll(bodies_, bodyi)
    {
        freeMotion* fm = freeMotionOf(bodyi);
        if (!fm) continue;

        const ibmBody& b = bodies_[bodyi];
        scalar V = -1;
        symmTensor I(Zero);
        if (!b.massProperties(V, I))
        {
            V = -1;
            I = Zero;
        }
        fm->setMassProperties(V, I);

        if (statePtr && statePtr->found(b.name()))
        {
            fm->readState(statePtr->subDict(b.name()));
            Info<< "    body " << b.name() << ": free motion state read from "
                << stateFile << endl;
        }

        Info<< "    body " << b.name() << ": free motion, "
            << freeMotion::couplingNames[fm->coupling()] << " coupling, "
            << "density " << fm->density() << ", volume " << fm->volume()
            << ", inertia " << fm->inertia()
            << (fm->implicit()
                ? ", added mass " + Foam::name(fm->addedMass())
                : std::string())
            << endl;
    }
}


Foam::vector Foam::cutCellGeometry::meanPressureGradient() const
{
    // Set by ibmMeanVelocityForce: the force per unit mass that drives the
    // flow, i.e. minus the mean pressure gradient the periodic pressure of
    // the solver leaves out
    const word name("ibmMeanPressureGradient");
    if (mesh_.foundObject<uniformDimensionedVectorField>(name))
    {
        return mesh_.lookupObject<uniformDimensionedVectorField>(name).value();
    }
    return Zero;
}


void Foam::cutCellGeometry::freeBodyLoads
(
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& nu,
    const volScalarField* nuMolecularPtr,
    List<vector>& F,
    List<vector>& M
) const
{
    List<vector> Fp, Fv, Mp, Mv;
    forces(p, U, nu, Fp, Fv, Mp, Mv, nuMolecularPtr);
    const vector G(meanPressureGradient());

    F.setSize(bodies_.size());
    M.setSize(bodies_.size());
    forAll(bodies_, bodyi)
    {
        F[bodyi] = Fp[bodyi] + Fv[bodyi];
        M[bodyi] = Mp[bodyi] + Mv[bodyi];
        const freeMotion* fm = freeMotionOf(bodyi);
        if (fm)
        {
            F[bodyi] += fm->bodyForce();
            if (fm->meanPressureGradient())
            {
                F[bodyi] += fm->volume()*G;
            }
        }
    }
}


void Foam::cutCellGeometry::reportFreeBodies
(
    const word& stage,
    const List<vector>& F,
    const List<vector>& M
) const
{
    forAll(bodies_, bodyi)
    {
        const freeMotion* fm = freeMotionOf(bodyi);
        if (!fm) continue;
        const ibmBody& b = bodies_[bodyi];
        const scalar t = mesh_.time().value();
        Info<< "IBM free body " << b.name() << " " << stage
            << ": t " << t
            << " displacement " << fm->displacement(t)
            << " velocity " << fm->velocity(t)
            << " angularVelocity " << fm->omega()
            << " force " << F[bodyi]
            << " moment " << M[bodyi]
            << " change " << fm->residual();
        if (fm->implicit())
        {
            Info<< " addedMass " << fm->effectiveAddedMass()
                << " addedInertia " << fm->effectiveAddedInertia();
        }
        if (fm->rotates())
        {
            Info<< " orientation " << fm->orientation(t);
        }
        Info<< endl;
    }
}


void Foam::cutCellGeometry::predictFreeBodies
(
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& nu,
    const volScalarField* nuMolecularPtr
)
{
    List<vector> F, M;
    freeBodyLoads(p, U, nu, nuMolecularPtr, F, M);
    const scalar dt = mesh_.time().deltaTValue();
    forAll(bodies_, bodyi)
    {
        freeMotion* fm = freeMotionOf(bodyi);
        if (fm)
        {
            fm->predict(dt, F[bodyi], M[bodyi]);
        }
    }
    reportFreeBodies("predictor", F, M);
}


Foam::scalar Foam::cutCellGeometry::correctFreeBodies
(
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& nu,
    const volScalarField* nuMolecularPtr
)
{
    List<vector> F, M;
    freeBodyLoads(p, U, nu, nuMolecularPtr, F, M);
    const scalar dt = mesh_.time().deltaTValue();
    scalar change = 0;
    forAll(bodies_, bodyi)
    {
        freeMotion* fm = freeMotionOf(bodyi);
        if (fm && fm->implicit())
        {
            fm->correct(dt, F[bodyi], M[bodyi]);
            change = max(change, fm->residual());
        }
    }
    reportFreeBodies("corrector", F, M);
    return change;
}


void Foam::cutCellGeometry::reportFreeBodyLoads
(
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& nu,
    const volScalarField* nuMolecularPtr
) const
{
    if (!freeBodies())
    {
        return;
    }
    List<vector> F, M;
    freeBodyLoads(p, U, nu, nuMolecularPtr, F, M);
    reportFreeBodies("final", F, M);
}


void Foam::cutCellGeometry::writeFreeBodies() const
{
    if (!freeBodies() || !mesh_.time().writeTime())
    {
        return;
    }
    IOdictionary dict
    (
        IOobject
        (
            "ibmFreeBodies",
            mesh_.time().timeName(),
            "uniform",
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    forAll(bodies_, bodyi)
    {
        const freeMotion* fm = freeMotionOf(bodyi);
        if (fm)
        {
            dict.add(bodies_[bodyi].name(), fm->state());
        }
    }
    dict.regIOobject::write();
}


void Foam::cutCellGeometry::bodyVolumes(List<scalar>& Vb) const
{
    Vb.setSize(bodies_.size());
    Vb = 0;
    const scalarField& V = mesh_.V();
    forAll(cellBody_, celli)
    {
        const label bodyi = cellBody_[celli];
        if (bodyi >= 0)
        {
            Vb[bodyi] += (1 - alpha_[celli])*V[celli];
        }
    }
    forAll(Vb, bodyi)
    {
        reduce(Vb[bodyi], sumOp<scalar>());
    }
}


// ************************************************************************* //
