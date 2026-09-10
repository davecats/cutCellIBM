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


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::cutCellGeometry::snGradCorrection
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceField;

    tmp<SurfaceField> tcorr
    (
        new SurfaceField
        (
            IOobject
            (
                "cutCellSnGradCorr(" + vf.name() + ')',
                vf.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<Type>(vf.dimensions()/dimLength, Zero)
        )
    );
    SurfaceField& corr = tcorr.ref();
    corr.setOriented();

    if (!secondOrder_)
    {
        return tcorr;
    }

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        const tmp<volVectorField> tg(lsqGrad(vf.component(cmpt)));
        const volVectorField& g = tg();

        surfaceScalarField c
        (
            IOobject("c", vf.instance(), mesh_),
            mesh_,
            dimensionedScalar(vf.dimensions()/dimLength, Zero)
        );

        forAll(own, facei)
        {
            if (theta_[facei] <= 0 || solidCells_.test(own[facei]) || solidCells_.test(nei[facei]))
            {
                continue;
            }
            const vector n(Sf[facei]/magSf[facei]);
            const vector& d = dCentroid_[facei];
            const vector gf(0.5*(g[own[facei]] + g[nei[facei]]));
            c[facei] = (gf & n) - (gf & d)*cutDeltaCoeffs_[facei];
        }
        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& fvp = mesh_.boundary()[patchi];
            if (!fvp.coupled()) continue;
            const labelUList& fc = fvp.faceCells();
            const vectorField gN(g.boundaryField()[patchi].patchNeighbourField());
            const fvsPatchScalarField& tp = theta_.boundaryField()[patchi];
            const fvsPatchVectorField& dp = dCentroid_.boundaryField()[patchi];
            const fvsPatchScalarField& dcp = cutDeltaCoeffs_.boundaryField()[patchi];
            const fvsPatchVectorField& Sfp = Sf.boundaryField()[patchi];
            const fvsPatchScalarField& magSfp = magSf.boundaryField()[patchi];
            fvsPatchScalarField& cp = c.boundaryFieldRef()[patchi];
            forAll(fvp, i)
            {
                if (tp[i] <= 0 || solidCells_.test(fc[i])) continue;
                const vector n(Sfp[i]/magSfp[i]);
                const vector gf(0.5*(g[fc[i]] + gN[i]));
                cp[i] = (gf & n) - (gf & dp[i])*dcp[i];
            }
        }

        corr.replace(cmpt, c);
    }

    return tcorr;
}
