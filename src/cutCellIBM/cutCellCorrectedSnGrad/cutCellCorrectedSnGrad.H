/*---------------------------------------------------------------------------*\
    cutCellIBM: cut-cell immersed boundary method for OpenFOAM

Class
    Foam::fv::cutCellCorrectedSnGrad

Description
    Surface-normal gradient consistent with the fluid centroids of cut
    cells: the implicit part uses the cut-cell delta coefficients 1/(n.d) of
    cutCellGeometry (d the vector between the fluid centroids of the two
    cells) and the explicit correction is g_f.n - (g_f.d)/(n.d) with g_f the
    interpolated least-squares gradient on the centroids (over-relaxed
    non-orthogonal form). On uncut cells and when the geometry is not run
    with "secondOrder true" it reduces to the plain corrected scheme on an
    orthogonal mesh.

    Select with
        snGradSchemes    { default cutCellCorrected; }
        laplacianSchemes { default Gauss linear cutCellCorrected; }
    so that the momentum and pressure Laplacians, their fluxes and
    fvc::snGrad are all consistent.

\*---------------------------------------------------------------------------*/

#ifndef cutCellCorrectedSnGrad_H
#define cutCellCorrectedSnGrad_H

#include "snGradScheme.H"

namespace Foam
{
namespace fv
{

template<class Type>
class cutCellCorrectedSnGrad
:
    public snGradScheme<Type>
{
    //- No copy assignment
    void operator=(const cutCellCorrectedSnGrad&) = delete;

public:

    TypeName("cutCellCorrected");

    cutCellCorrectedSnGrad(const fvMesh& mesh)
    :
        snGradScheme<Type>(mesh)
    {}

    cutCellCorrectedSnGrad(const fvMesh& mesh, Istream&)
    :
        snGradScheme<Type>(mesh)
    {}

    virtual ~cutCellCorrectedSnGrad() = default;

    //- Cut-cell delta coefficients
    virtual tmp<surfaceScalarField> deltaCoeffs
    (
        const GeometricField<Type, fvPatchField, volMesh>&
    ) const;

    //- Explicit correction is used when the geometry is second order
    virtual bool corrected() const;

    //- Explicit correction
    virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
    correction(const GeometricField<Type, fvPatchField, volMesh>&) const;
};

} // End namespace fv
} // End namespace Foam

#ifdef NoRepository
    #include "cutCellCorrectedSnGrad.C"
#endif

#endif
