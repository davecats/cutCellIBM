#include "ibmMeanVelocityForce.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(ibmMeanVelocityForce, 0);
    addToRunTimeSelectionTable(option, ibmMeanVelocityForce, dictionary);
}
}


Foam::fv::ibmMeanVelocityForce::ibmMeanVelocityForce
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    meanVelocityForce(sourceName, modelType, dict, mesh),
    Vfluid_(0)
{
    const scalarField& cv = mesh_.V();
    const volScalarField& a = alpha();
    forAll(cells_, i)
    {
        Vfluid_ += a[cells_[i]]*cv[cells_[i]];
    }
    reduce(Vfluid_, sumOp<scalar>());
    Info<< "    fluid volume of the selection = " << Vfluid_ << nl << endl;
}


const Foam::volScalarField& Foam::fv::ibmMeanVelocityForce::alpha() const
{
    return mesh_.lookupObject<volScalarField>("ibmAlpha");
}


Foam::scalar Foam::fv::ibmMeanVelocityForce::magUbarAve
(
    const volVectorField& U
) const
{
    const scalarField& cv = mesh_.V();
    const volScalarField& a = alpha();
    scalar magUbarAve = 0.0;
    forAll(cells_, i)
    {
        const label celli = cells_[i];
        magUbarAve += (flowDir_ & U[celli])*a[celli]*cv[celli];
    }
    reduce(magUbarAve, sumOp<scalar>());
    return magUbarAve/Vfluid_;
}


void Foam::fv::ibmMeanVelocityForce::correct(volVectorField& U)
{
    const scalarField& rAU = rAPtr_();
    const scalarField& cv = mesh_.V();
    const volScalarField& a = alpha();

    // Response of the bulk velocity to a unit change of the force: the
    // source of cell P is gradP alpha V, so du_P = rAU alpha dGradP
    scalar rAUave = 0.0;
    forAll(cells_, i)
    {
        const label celli = cells_[i];
        rAUave += rAU[celli]*a[celli]*a[celli]*cv[celli];
    }
    reduce(rAUave, sumOp<scalar>());
    rAUave /= Vfluid_;

    const scalar magUbarAve = this->magUbarAve(U);

    // Accumulate: correct() is called after the momentum predictor and
    // after every pressure corrector, and each increment must be kept
    const scalar dGradPnew = relaxation_*(mag(Ubar_) - magUbarAve)/rAUave;
    dGradP_ += dGradPnew;
    forAll(cells_, i)
    {
        const label celli = cells_[i];
        U[celli] += flowDir_*rAU[celli]*a[celli]*dGradPnew;
    }
    U.correctBoundaryConditions();

    const scalar gradP = gradP0_ + dGradP_;

    Info<< "Pressure gradient source: uncorrected Ubar = " << magUbarAve
        << ", pressure gradient = " << gradP << endl;

    writeProps(gradP);
}


void Foam::fv::ibmMeanVelocityForce::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    volVectorField::Internal Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldi] + "Sup",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(eqn.dimensions()/dimVolume, Zero)
    );

    const scalar gradP = gradP0_ + dGradP_;
    const volScalarField& a = alpha();

    forAll(cells_, i)
    {
        const label celli = cells_[i];
        Su[celli] = flowDir_*gradP*a[celli];
    }

    eqn += Su;
}


void Foam::fv::ibmMeanVelocityForce::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    this->addSup(eqn, fieldi);
}


void Foam::fv::ibmMeanVelocityForce::constrain
(
    fvMatrix<vector>& eqn,
    const label
)
{
    // Response of a cell to a uniform change of the source: (A - H1)^-1,
    // guarded against a vanishing denominator
    const volScalarField A(eqn.A());
    tmp<volScalarField> trA(1.0/max(A - eqn.H1(), 1e-3*A));

    if (!rAPtr_)
    {
        rAPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    name_ + ":rA",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                trA
            )
        );
    }
    else
    {
        rAPtr_() = trA;
    }

    gradP0_ += dGradP_;
    dGradP_ = 0.0;
}
