#include "cutCellCorrectedSnGrad.H"
#include "cutCellGeometry.H"
#include "fvMesh.H"

namespace
{
    const Foam::cutCellGeometry& geometry(const Foam::fvMesh& mesh)
    {
        return mesh.lookupObject<Foam::cutCellGeometry>("ibmProperties");
    }
}


template<class Type>
Foam::tmp<Foam::surfaceScalarField>
Foam::fv::cutCellCorrectedSnGrad<Type>::deltaCoeffs
(
    const GeometricField<Type, fvPatchField, volMesh>&
) const
{
    const cutCellGeometry& ibm = geometry(this->mesh());
    if (ibm.secondOrder())
    {
        return ibm.cutDeltaCoeffs();
    }
    return this->mesh().nonOrthDeltaCoeffs();
}


template<class Type>
bool Foam::fv::cutCellCorrectedSnGrad<Type>::corrected() const
{
    return geometry(this->mesh()).secondOrder();
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fv::cutCellCorrectedSnGrad<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return geometry(this->mesh()).snGradCorrection(vf);
}
