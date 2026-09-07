#include "cutCellGeometry.H"

template<class Type>
void Foam::cutCellGeometry::zeroSolid
(
    GeometricField<Type, fvPatchField, volMesh>& f
) const
{
    Field<Type>& fI = f.primitiveFieldRef();
    for (const label celli : solidCells_)
    {
        fI[celli] = Zero;
    }
    f.correctBoundaryConditions();
}
