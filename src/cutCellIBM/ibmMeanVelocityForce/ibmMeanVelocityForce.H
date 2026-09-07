/*---------------------------------------------------------------------------*\
    cutCellIBM: cut-cell immersed boundary method for OpenFOAM

Class
    Foam::fv::ibmMeanVelocityForce

Description
    meanVelocityForce (constant flow rate forcing) for the cut-cell immersed
    boundary: the bulk velocity is averaged over the fluid volume alpha V
    and the body force acts on the fluid volume of each cell only, so that
    solid and cut cells receive no spurious momentum. The fluid fraction is
    the field ibmAlpha registered by cutCellGeometry.

    The force increment is computed with the mobility of a *uniform* velocity
    change, 1/(A - H1) (the SIMPLEC coefficient), instead of 1/A. For a
    uniform force the neighbour couplings of the momentum matrix cancel, so
    this is the exact response of the bulk velocity to the force; 1/A
    underestimates it by the diffusion number nu dt/h^2 when the viscous
    term is implicit, which makes the standard controller overshoot.

    Usage, in constant/fvOptions:
    \verbatim
    momentumSource
    {
        type            ibmMeanVelocityForce;
        selectionMode   all;
        fields          (U);
        Ubar            (0.5 0 0);
        relaxation      1.0;
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef ibmMeanVelocityForce_H
#define ibmMeanVelocityForce_H

#include "meanVelocityForce.H"

namespace Foam
{
namespace fv
{

class ibmMeanVelocityForce
:
    public meanVelocityForce
{
    // Private data

        //- Fluid volume of the selected cells
        scalar Vfluid_;


    // Private member functions

        const volScalarField& alpha() const;


protected:

    //- Fluid-volume average of the velocity along the flow direction
    virtual scalar magUbarAve(const volVectorField& U) const;


public:

    TypeName("ibmMeanVelocityForce");

    ibmMeanVelocityForce
    (
        const word& sourceName,
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh
    );

    virtual ~ibmMeanVelocityForce() = default;


    // Member Functions

        virtual void correct(volVectorField& U);

        virtual void addSup(fvMatrix<vector>& eqn, const label fieldi);

        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Store the mobility of a uniform velocity change, 1/(A - H1)
        virtual void constrain(fvMatrix<vector>& eqn, const label fieldi);
};

} // End namespace fv
} // End namespace Foam

#endif
