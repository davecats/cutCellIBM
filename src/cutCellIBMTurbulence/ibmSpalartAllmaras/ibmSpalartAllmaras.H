/*---------------------------------------------------------------------------*\
    cutCellIBM: cut-cell immersed boundary method for OpenFOAM

Class
    Foam::RASModels::ibmSpalartAllmaras

Description
    Spalart-Allmaras RAS model (OpenFOAM's SpalartAllmaras, unchanged
    coefficients and closure functions) whose nuTilda equation sees the
    immersed walls of the cut-cell geometry, assembled as the passive scalars
    of the library:

        ddt(alpha nuTilda) + div(phi, nuTilda)
      - div(theta (nu + nuTilda)/sigma grad nuTilda)       (wet faces only)
      + nu/sigma Awall/dWall nuTilda                        (nuTilda_w = 0)
      - Cb2/sigma alpha |grad nuTilda|^2
      = alpha Cb1 Stilda nuTilda - alpha Cw1 fw (nuTilda/y)^2

    with alpha, theta the fluid fractions, phi the solver's wet-face flux,
    the Dirichlet wall term on the no-slip walls (zero gradient on free-slip
    walls), blanking of the solid cells, and the gradients of U and nuTilda
    averaged over the fluid with their wall values (the wall velocity, 0).
    The wall distance y must come from the immersed walls:
    \\verbatim
    system/fvSchemes:   wallDist { method ibmLevelSet; }
    \\endverbatim
    The wall shear of the momentum equation (wall functions) is a matter of
    the solver (cutCellGeometry::wallViscosity); nut is nuTilda fv1 as in
    the body-fitted model.

    On walls the mesh does not resolve (wallFunction spalartAllmaras),
    nuTilda of the wall cells is drawn towards the model's own wall
    solution, nuTilda = kappa u_tau y, which holds through the viscous
    sublayer and the log layer alike, with the weight w of the wall law from
    the cell size in wall units (cutCellGeometry::wallLaw; hPlusResolved,
    hPlusModelled in ibmProperties): w = 0 leaves resolved walls alone,
    w = 1 fixes the value where the wall-cell gradients cannot carry the
    production of the log layer.

    Negative nuTilda (negativeNuTilda true, the default): the SA-neg model
    of Allmaras, Johnson & Spalart (ICCFD7-1902, 2012). Where nuTilda < 0
    the production is Cb1 (1 - Ct3) Omega nuTilda, the destruction
    -Cw1 (nuTilda/y)^2 (it drives nuTilda back to zero), the diffusivity
    (nu + fn nuTilda)/sigma with fn = (Cn1 + chi^3)/(Cn1 - chi^3), and
    nut = 0; nuTilda is not clipped. Where nuTilda >= 0 the model is the
    standard one, unchanged to the bit. Clipping (negativeNuTilda false,
    the base model) leaves a residual that does not converge where the
    solution wants to undershoot, e.g. at the edge of a wake on a stretched
    mesh.
    \verbatim
    ibmSpalartAllmarasCoeffs { negativeNuTilda true; Ct3 1.2; Cn1 16; }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef ibmSpalartAllmaras_H
#define ibmSpalartAllmaras_H

#include "SpalartAllmaras.H"

namespace Foam
{

class cutCellGeometry;

namespace RASModels
{

template<class BasicTurbulenceModel>
class ibmSpalartAllmaras
:
    public SpalartAllmaras<BasicTurbulenceModel>
{
protected:

        //- The cut-cell geometry of the immersed bodies
        const cutCellGeometry& ibm_;

        //- SA-neg: negative nuTilda allowed (no clipping)
        Switch negative_;

        //- SA-neg coefficients
        dimensionedScalar Ct3_;
        dimensionedScalar Cn1_;

        //- Wall velocity of every cell (zero for bodies at rest)
        tmp<vectorField> wallVelocity() const;

        //- Diffusivity of nuTilda, (nu + fn nuTilda)/sigma (SA-neg)
        tmp<volScalarField> DnuTildaEffNeg() const;

        //- nut = nuTilda fv1, zero where nuTilda < 0
        virtual void correctNut();

private:

    ibmSpalartAllmaras(const ibmSpalartAllmaras&) = delete;
    void operator=(const ibmSpalartAllmaras&) = delete;

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    TypeName("ibmSpalartAllmaras");

    ibmSpalartAllmaras
    (
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const transportModel& transport,
        const word& propertiesName = turbulenceModel::propertiesName,
        const word& type = typeName
    );

    virtual ~ibmSpalartAllmaras() = default;

    //- Re-read the coefficients
    virtual bool read();

    //- Solve the nuTilda equation on the cut mesh and update nut
    virtual void correct();
};

} // End namespace RASModels
} // End namespace Foam

#ifdef NoRepository
    #include "ibmSpalartAllmaras.C"
#endif

#endif
