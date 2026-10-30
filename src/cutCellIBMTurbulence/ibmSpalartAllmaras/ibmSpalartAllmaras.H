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

        //- Wall velocity of every cell (zero for bodies at rest)
        tmp<vectorField> wallVelocity() const;

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

    //- Solve the nuTilda equation on the cut mesh and update nut
    virtual void correct();
};

} // End namespace RASModels
} // End namespace Foam

#ifdef NoRepository
    #include "ibmSpalartAllmaras.C"
#endif

#endif
