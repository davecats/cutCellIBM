#include "cutCellGeometry.H"
#include "surfaceInterpolate.H"

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


template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::cutCellGeometry::fluidGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const Field<Type>& wallValue
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradField;

    tmp<GradField> tg
    (
        new GradField
        (
            IOobject
            (
                "fluidGrad(" + vf.name() + ')',
                vf.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<GradType>(vf.dimensions()/dimLength, Zero),
            "zeroGradient"
        )
    );
    Field<GradType>& gI = tg.ref().primitiveFieldRef();

    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tvff
    (
        fvc::interpolate(vf)
    );
    const GeometricField<Type, fvsPatchField, surfaceMesh>& vff = tvff();
    const surfaceVectorField& Sf = mesh_.Sf();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    forAll(own, facei)
    {
        const GradType v(theta_[facei]*(Sf[facei]*vff[facei]));
        gI[own[facei]] += v;
        gI[nei[facei]] -= v;
    }
    forAll(mesh_.boundary(), patchi)
    {
        const labelUList& fc = mesh_.boundary()[patchi].faceCells();
        const fvsPatchScalarField& tp = theta_.boundaryField()[patchi];
        const fvsPatchVectorField& Sfp = Sf.boundaryField()[patchi];
        const fvsPatchField<Type>& vfp = vff.boundaryField()[patchi];
        forAll(tp, i)
        {
            gI[fc[i]] += tp[i]*(Sfp[i]*vfp[i]);
        }
    }

    // Wall value through closure: sum_f theta_f Sf + Sw = 0, so a uniform
    // field equal to the wall value has zero gradient
    const vectorField& Sw = Sw_.primitiveField();
    const scalarField& V = mesh_.V();
    forAll(gI, celli)
    {
        if (alpha_[celli] > 0)
        {
            gI[celli] = (gI[celli] + Sw[celli]*wallValue[celli])
                       /(alpha_[celli]*V[celli]);
        }
        else
        {
            gI[celli] = Zero;
        }
    }

    tg.ref().correctBoundaryConditions();
    return tg;
}


template<class Type>
void Foam::cutCellGeometry::linkCells(fvMatrix<Type>& eqn) const
{
    if (mergeAlpha_ <= 0)
    {
        return;
    }

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    scalarField& diag = eqn.diag();
    const bool symmetric = eqn.symmetric();
    scalarField& upper = eqn.upper();
    scalarField* lowerPtr = symmetric ? nullptr : &eqn.lower();

    forAll(linkFace_, celli)
    {
        const label facei = linkFace_[celli];
        if (facei < 0) continue;
        const label other = (own[facei] == celli) ? nei[facei] : own[facei];

        // A conductance with the sign of the diagonal (positive for
        // -laplacian, negative for the pressure laplacian)
        const scalar sgn = (diag[celli] + diag[other] >= 0) ? 1 : -1;
        const scalar c =
            mergePenalty_*max(mag(diag[celli]), mag(diag[other]));

        diag[celli] += sgn*c;
        diag[other] += sgn*c;
        upper[facei] -= sgn*c;
        if (lowerPtr)
        {
            (*lowerPtr)[facei] -= sgn*c;
        }
    }
}


template<class Type>
void Foam::cutCellGeometry::groupAverage(Field<Type>& f) const
{
    if (mergeAlpha_ <= 0)
    {
        return;
    }

    const scalarField& alpha = alpha_.primitiveField();
    const scalarField& V = mesh_.V();
    Field<Type> sum(f.size(), Zero);
    scalarField vol(f.size(), Zero);
    forAll(linkMaster_, celli)
    {
        const label m = linkMaster_[celli];
        if (m < 0) continue;
        sum[m] += alpha[celli]*V[celli]*f[celli];
        vol[m] += alpha[celli]*V[celli];
    }
    forAll(vol, m)
    {
        if (vol[m] > 0)
        {
            sum[m] += alpha[m]*V[m]*f[m];
            vol[m] += alpha[m]*V[m];
            f[m] = sum[m]/vol[m];
        }
    }
    forAll(linkMaster_, celli)
    {
        const label m = linkMaster_[celli];
        if (m >= 0)
        {
            f[celli] = f[m];
        }
    }
}


template<class Type>
void Foam::cutCellGeometry::groupMean
(
    GeometricField<Type, fvPatchField, volMesh>& f
) const
{
    if (mergeAlpha_ <= 0)
    {
        return;
    }

    const scalarField& V = mesh_.V();
    Field<Type>& fI = f.primitiveFieldRef();
    Field<Type> sum(fI.size(), Zero);
    scalarField vol(fI.size(), Zero);
    forAll(linkMaster_, celli)
    {
        const label m = (linkMaster_[celli] >= 0) ? linkMaster_[celli] : celli;
        sum[m] += V[celli]*fI[celli];
        vol[m] += V[celli];
    }
    forAll(linkMaster_, celli)
    {
        const label m = (linkMaster_[celli] >= 0) ? linkMaster_[celli] : celli;
        fI[celli] = sum[m]/vol[m];
    }
    f.correctBoundaryConditions();
}


template<class Type>
void Foam::cutCellGeometry::mergeA
(
    const fvMatrix<Type>& eqn,
    volScalarField& A
) const
{
    if (mergeAlpha_ <= 0)
    {
        return;
    }

    const scalarField& V = mesh_.V();
    const scalarField& upper = eqn.upper();
    const scalarField& lower = eqn.lower();
    scalarField& AI = A.primitiveFieldRef();
    scalarField sum(AI.size(), Zero);
    scalarField vol(AI.size(), Zero);
    forAll(linkMaster_, celli)
    {
        const label m = (linkMaster_[celli] >= 0) ? linkMaster_[celli] : celli;
        sum[m] += V[celli]*AI[celli];
        vol[m] += V[celli];
        const label facei = linkFace_[celli];
        if (linkMaster_[celli] >= 0 && facei >= 0)
        {
            // The coefficients of the internal face drop out of the merged
            // diagonal (A_i holds -sum of the off-diagonals of its row)
            sum[m] += upper[facei] + lower[facei];
        }
    }
    forAll(linkMaster_, celli)
    {
        const label m = (linkMaster_[celli] >= 0) ? linkMaster_[celli] : celli;
        AI[celli] = sum[m]/vol[m];
    }
    A.correctBoundaryConditions();
}


template<class Type>
void Foam::cutCellGeometry::mergeH
(
    const fvMatrix<Type>& eqn,
    GeometricField<Type, fvPatchField, volMesh>& H
) const
{
    if (mergeAlpha_ <= 0)
    {
        return;
    }

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const scalarField& V = mesh_.V();
    const scalarField& upper = eqn.upper();
    const scalarField& lower = eqn.lower();
    const Field<Type>& psi = eqn.psi().primitiveField();
    Field<Type>& HI = H.primitiveFieldRef();
    Field<Type> sum(HI.size(), Zero);
    scalarField vol(HI.size(), Zero);
    forAll(linkMaster_, celli)
    {
        const label m = (linkMaster_[celli] >= 0) ? linkMaster_[celli] : celli;
        sum[m] += V[celli]*HI[celli];
        vol[m] += V[celli];
        const label facei = linkFace_[celli];
        if (linkMaster_[celli] >= 0 && facei >= 0)
        {
            // H_own holds -upper psi_nei, H_nei -lower psi_own: inside the
            // merged cell, not a coupling to a neighbour
            sum[m] += upper[facei]*psi[nei[facei]] + lower[facei]*psi[own[facei]];
        }
    }
    forAll(linkMaster_, celli)
    {
        const label m = (linkMaster_[celli] >= 0) ? linkMaster_[celli] : celli;
        HI[celli] = sum[m]/vol[m];
    }
    H.correctBoundaryConditions();
}
