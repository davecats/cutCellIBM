#include "ibmMeanVelocityForce.H"
#include "fvMatrices.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Time.H"
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
    cellSetOption(sourceName, modelType, dict, mesh),
    Ubar_(coeffs_.get<vector>("Ubar")),
    flowDir_(Zero),
    relaxation_(coeffs_.getOrDefault<scalar>("relaxation", 1)),
    gradP0_(0),
    dGradP_(0),
    rAPtr_(nullptr),
    Vfluid_(0),
    steady_(false)
{

    coeffs_.readEntry("fields", fieldNames_);
    if (fieldNames_.size() != 1)
    {
        FatalErrorInFunction
            << "settings are:" << fieldNames_ << exit(FatalError);
    }
    fv::option::resetApplied();

    {
        const word scheme
        (
            mesh_.ddtScheme("ddt(ibmAlpha," + fieldNames_[0] + ')')
        );
        steady_ = (scheme == "steadyState");
    }

    if (mag(Ubar_) > 0)
    {
        flowDir_ = Ubar_/mag(Ubar_);
        if (coeffs_.found("flowDir"))
        {
            flowDir_ = normalised(coeffs_.get<vector>("flowDir"));
        }
    }
    else
    {
        if (!coeffs_.found("flowDir"))
        {
            FatalIOErrorInFunction(coeffs_)
                << "Ubar is zero: the direction of the force must be given"
                << " with flowDir" << exit(FatalIOError);
        }
        flowDir_ = normalised(coeffs_.get<vector>("flowDir"));
    }

    // Read the initial force from file if it exists (restart)
    IFstream propsFile
    (
        mesh_.time().timePath()/"uniform"/(name_ + "Properties")
    );
    if (propsFile.good())
    {
        Info<< "    Reading pressure gradient from file" << endl;
        dictionary propsDict(propsFile);
        propsDict.readEntry("gradient", gradP0_);
    }

    const scalarField& cv = mesh_.V();
    const volScalarField& a = alpha();
    forAll(cells_, i)
    {
        Vfluid_ += a[cells_[i]]*cv[cells_[i]];
    }
    reduce(Vfluid_, sumOp<scalar>());

    Info<< "    Ubar = " << Ubar_ << ", flow direction " << flowDir_ << nl
        << "    initial pressure gradient = " << gradP0_ << nl
        << "    fluid volume of the selection = " << Vfluid_ << nl << endl;
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


void Foam::fv::ibmMeanVelocityForce::writeProps(const scalar gradP) const
{
    if (mesh_.time().writeTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                name_ + "Properties",
                mesh_.time().timeName(),
                "uniform",
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        propsDict.add("gradient", gradP);
        propsDict.regIOobject::write();
    }
}


void Foam::fv::ibmMeanVelocityForce::correct(volVectorField& U)
{
    const scalarField& rAU = rAPtr_();
    const scalarField& cv = mesh_.V();
    const volScalarField& a = alpha();

    // Response of the bulk velocity to a unit change of the force: the
    // source of cell P is gradP alpha V, so du_P = rAU alpha dGradP. The
    // fluid volume of a moving body changes in time, so it is re-summed.
    scalar rAUave = 0.0;
    scalar Vfluid = 0.0;
    forAll(cells_, i)
    {
        const label celli = cells_[i];
        rAUave += rAU[celli]*a[celli]*a[celli]*cv[celli];
        Vfluid += a[celli]*cv[celli];
    }
    reduce(rAUave, sumOp<scalar>());
    reduce(Vfluid, sumOp<scalar>());
    Vfluid_ = Vfluid;
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

    Info<< "Pressure gradient source: uncorrected Ubar = " << magUbarAve
        << ", pressure gradient = " << gradP() << endl;

    writeProps(gradP());
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

    const volScalarField& a = alpha();
    forAll(cells_, i)
    {
        const label celli = cells_[i];
        Su[celli] = flowDir_*gradP()*a[celli];
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
    // Response of a cell to a uniform change of the source: (A - H1)^-1.
    // The row sum A - H1 can vanish or turn negative in cells a moving
    // body sweeps through by more than they hold (alpha^n < 0), so it is
    // floored at the transient diagonal alpha/deltaT, as in the pressure
    // equation; in steady state (relaxed equations) at a fraction of A.
    const volScalarField A(eqn.A());
    tmp<volScalarField> tfloor(1e-3*A);
    if (!steady_)
    {
        tfloor = max(tfloor(), alpha()/mesh_.time().deltaT());
    }
    const volScalarField AmH1(A - eqn.H1());
    tmp<volScalarField> trA(1.0/max(AmH1, tfloor()));

    if (debug)
    {
        label nFloored = 0;
        forAll(AmH1, celli)
        {
            if (AmH1[celli] < tfloor()[celli]) ++nFloored;
        }
        reduce(nFloored, sumOp<label>());
        Info<< "    " << type() << ": mobility floored in " << nFloored
            << " cells, min(A - H1)/max(A) = "
            << gMin(AmH1.primitiveField())/gMax(A.primitiveField()) << endl;
    }

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


bool Foam::fv::ibmMeanVelocityForce::read(const dictionary& dict)
{
    NotImplemented;
    return false;
}
