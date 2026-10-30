/*---------------------------------------------------------------------------*\
    cutCellIBM: cut-cell immersed boundary method for OpenFOAM

Class
    Foam::patchDistMethods::ibmLevelSet

Description
    Wall distance for turbulence models that see immersed walls: the
    smaller of the distance to the wall patches of the mesh (meshWave, if
    there are any) and the distance to the immersed bodies, taken from their
    level sets (distance functions for every level-set body type). The
    distance is evaluated at the fluid centroids (the centre of a cut cell
    may lie in the solid) and floored at a small positive value in the
    solid. Requires the cut-cell geometry
    (constant/ibmProperties) to be constructed before the turbulence model,
    as the ibm solvers do.

    Selected in system/fvSchemes:
    \verbatim
    wallDist
    {
        method ibmLevelSet;
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef ibmLevelSetPatchDistMethod_H
#define ibmLevelSetPatchDistMethod_H

#include "patchDistMethod.H"

namespace Foam
{
namespace patchDistMethods
{

class ibmLevelSet
:
    public patchDistMethod
{
    //- No copy construct
    ibmLevelSet(const ibmLevelSet&) = delete;

    //- No copy assignment
    void operator=(const ibmLevelSet&) = delete;

public:

    TypeName("ibmLevelSet");

    ibmLevelSet
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const labelHashSet& patchIDs
    );

    //- Correct the given distance-to-wall field
    virtual bool correct(volScalarField& y);

    //- Correct the distance-to-wall and the normal-to-wall fields
    virtual bool correct(volScalarField& y, volVectorField& n);
};

} // End namespace patchDistMethods
} // End namespace Foam

#endif
