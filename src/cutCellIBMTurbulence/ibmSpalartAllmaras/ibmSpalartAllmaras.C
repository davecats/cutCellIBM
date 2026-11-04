#include "ibmSpalartAllmaras.H"
#include "cutCellGeometry.H"
#include "fvOptions.H"
#include "bound.H"

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
ibmSpalartAllmaras<BasicTurbulenceModel>::ibmSpalartAllmaras
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    SpalartAllmaras<BasicTurbulenceModel>
    (
        alpha, rho, U, alphaRhoPhi, phi, transport, propertiesName, type
    ),
    ibm_
    (
        this->mesh_.template lookupObject<cutCellGeometry>("ibmProperties")
    ),
    negative_
    (
        Switch::getOrAddToDict("negativeNuTilda", this->coeffDict_, true)
    ),
    Ct3_
    (
        dimensioned<scalar>::getOrAddToDict("Ct3", this->coeffDict_, 1.2)
    ),
    Cn1_
    (
        dimensioned<scalar>::getOrAddToDict("Cn1", this->coeffDict_, 16)
    )
{
    const word method
    (
        this->mesh_.schemesDict().subOrEmptyDict("wallDist")
            .template getOrDefault<word>("method", "meshWave")
    );
    if (method != "ibmLevelSet")
    {
        FatalErrorInFunction
            << type << " needs the distance to the immersed walls: set "
            << "wallDist { method ibmLevelSet; } in system/fvSchemes "
            << "(found " << method << ")" << exit(FatalError);
    }

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


template<class BasicTurbulenceModel>
tmp<vectorField>
ibmSpalartAllmaras<BasicTurbulenceModel>::wallVelocity() const
{
    if (ibm_.moving())
    {
        return tmp<vectorField>(new vectorField(ibm_.Ub().primitiveField()));
    }
    return tmp<vectorField>(new vectorField(this->mesh_.nCells(), Zero));
}


template<class BasicTurbulenceModel>
bool ibmSpalartAllmaras<BasicTurbulenceModel>::read()
{
    if (SpalartAllmaras<BasicTurbulenceModel>::read())
    {
        negative_.readIfPresent("negativeNuTilda", this->coeffDict());
        Ct3_.readIfPresent(this->coeffDict());
        Cn1_.readIfPresent(this->coeffDict());
        return true;
    }
    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField>
ibmSpalartAllmaras<BasicTurbulenceModel>::DnuTildaEffNeg() const
{
    // fn = 1 for chi >= 0 (exactly: the base diffusivity), and
    // (Cn1 + chi^3)/(Cn1 - chi^3) for chi < 0 (denominator >= Cn1)
    const volScalarField chi(this->chi());
    const volScalarField chiNeg(min(chi, dimensionedScalar(dimless, Zero)));
    const volScalarField fn
    (
        pos0(chi) + neg(chi)*(Cn1_ + pow3(chiNeg))/(Cn1_ - pow3(chiNeg))
    );
    return tmp<volScalarField>::New
    (
        "DnuTildaEff",
        (this->nu() + this->nuTilda_*fn)/this->sigmaNut_
    );
}


template<class BasicTurbulenceModel>
void ibmSpalartAllmaras<BasicTurbulenceModel>::correctNut()
{
    // nut = nuTilda fv1 where nuTilda >= 0 (fv1 of chi >= 0, which is
    // regular), zero where it is negative
    const volScalarField chi(this->chi());
    this->nut_ =
        max(this->nuTilda_, dimensionedScalar(this->nuTilda_.dimensions(), Zero))
       *this->fv1(max(chi, dimensionedScalar(dimless, Zero)));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void ibmSpalartAllmaras<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    {
        fv::options& fvOptions(fv::options::New(this->mesh_));

        eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

        volScalarField& nuTilda = this->nuTilda_;
        const volScalarField& alpha = ibm_.alpha();
        const volScalarField nu(this->nu());

        // chi >= 0 in the closure functions of the standard model (fv1 is
        // singular at chi^3 = -Cv1^3); unchanged where nuTilda >= 0
        const volScalarField chi
        (
            max(this->chi(), dimensionedScalar(dimless, Zero))
        );
        const volScalarField fv1(this->fv1(chi));

        // Vorticity from the velocity gradient over the fluid, with the
        // wall velocity on the immersed wall. With linked small cut cells
        // the gradients are averaged over each merged group; with its
        // common wall distance and its linked values, every cell of a group
        // then has the source density of the merged cell
        // (over the merged control volume for linked small cut cells)
        volTensorField gradU(ibm_.fluidGrad(this->U_, wallVelocity()()));
        ibm_.groupAverage(gradU.primitiveFieldRef());
        const volScalarField::Internal Omega
        (
            ::sqrt(scalar(2))*mag(skew(gradU()))
        );
        const volScalarField::Internal Stilda
        (
            max
            (
                Omega
              + this->fv2(chi(), fv1())*nuTilda()/sqr(this->kappa_*this->y_),
                this->Cs_*Omega
            )
        );

        // nuTilda vanishes on the immersed no-slip walls
        volVectorField gradNuTilda
        (
            ibm_.fluidGrad(nuTilda, scalarField(this->mesh_.nCells(), Zero))
        );
        ibm_.groupAverage(gradNuTilda.primitiveFieldRef());

        // SA-neg where nuTilda < 0: production with Omega and (1 - Ct3),
        // destruction -Cw1 (nuTilda/y)^2 (fw replaced by -1). Where
        // nuTilda >= 0 both factors are the standard ones to the bit.
        const volScalarField::Internal production
        (
            negative_
          ? volScalarField::Internal
            (
                pos0(nuTilda())*Stilda
              + neg(nuTilda())*(1 - Ct3_)*Omega
            )
          : Stilda
        );
        const volScalarField::Internal fw(this->fw(Stilda));
        const volScalarField::Internal fwNeg
        (
            negative_
          ? volScalarField::Internal
            (
                pos0(nuTilda())*fw
              - neg(nuTilda())*dimensionedScalar(dimless, 1)
            )
          : fw
        );

        tmp<fvScalarMatrix> nuTildaEqn
        (
            fvm::ddt(alpha, nuTilda)
          + fvm::div(this->phi_, nuTilda)
          - fvm::laplacian
            (
                ibm_.thetaInterpolate
                (
                    negative_ ? DnuTildaEffNeg() : this->DnuTildaEff()
                ),
                nuTilda
            )
          + fvm::Sp(nu/this->sigmaNut_*ibm_.noSlipCoeff(), nuTilda)
          + fvm::Sp(nu/this->sigmaNut_*ibm_.blankCoeff(), nuTilda)
          - this->Cb2_/this->sigmaNut_*alpha*magSqr(gradNuTilda)
         ==
            this->Cb1_*alpha()*production*nuTilda()
          - fvm::Sp
            (
                this->Cw1_*alpha()*fwNeg*nuTilda()/sqr(this->y_),
                nuTilda
            )
          + fvOptions(nuTilda)
        );

        if (ibm_.ibmConvection())
        {
            nuTildaEqn.ref() += ibm_.convectionSource(this->phi_, nuTilda);
        }

        nuTildaEqn.ref().relax();

        // Wall cells of walls the mesh does not resolve: nuTilda towards
        // the wall solution of the model, kappa u_tau y, which holds down to
        // the wall (fv1 is built for it), with the weight w of the wall law
        // (cutCellGeometry::wallLaw): a conductance w/(1 - w) times the
        // diagonal, so w = 0 (resolved walls) leaves the equation unchanged
        // and w = 1 fixes the value
        {
            scalarField uTau, w;
            ibm_.wallLaw(this->U_, nu, uTau, w);
            const scalarField& y = this->y_.field();
            // With the coefficients of the coupled faces: the same weight
            // in every decomposition
            const scalarField diag0(ibm_.coupledDiag(nuTildaEqn()));
            scalarField& diag = nuTildaEqn.ref().diag();
            scalarField& source = nuTildaEqn.ref().source();
            for (const label celli : ibm_.wallCells())
            {
                if (w[celli] <= 0) continue;
                const scalar c = diag0[celli]*w[celli]/max(1 - w[celli], 1e-3);
                diag[celli] += c;
                source[celli] += c*this->kappa_.value()*uTau[celli]*y[celli];
            }
        }

        ibm_.linkCells(nuTildaEqn.ref());
        fvOptions.constrain(nuTildaEqn.ref());
        solve(nuTildaEqn);
        fvOptions.correct(nuTilda);
        if (!negative_)
        {
            bound(nuTilda, dimensionedScalar(nuTilda.dimensions(), Zero));
        }
        nuTilda.correctBoundaryConditions();
    }

    this->correctNut();
}

} // End namespace RASModels
} // End namespace Foam
