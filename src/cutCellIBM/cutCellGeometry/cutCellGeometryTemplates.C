#include "cutCellGeometry.H"
#include "surfaceInterpolate.H"
#include "fvcSurfaceIntegrate.H"
#include "convectionScheme.H"
#include "syncTools.H"
#include "PstreamBuffers.H"

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
    const volVectorField& C = mesh_.C();
    const surfaceVectorField& Cf = mesh_.Cf();
    const vectorField& off = valueOffset_.primitiveField();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        const tmp<volScalarField> tvc(vf.component(cmpt));
        const volScalarField& vc = tvc();
        const tmp<volVectorField> tg(lsqGrad(vc, true));
        const volVectorField& g = tg();

        surfaceScalarField c
        (
            IOobject("c", vf.instance(), mesh_),
            mesh_,
            dimensionedScalar(vf.dimensions()/dimLength, Zero)
        );

        // Only the faces where the bracket does not vanish (snGradFaces_;
        // on the other open faces d is parallel to n, n.d = 1/cutDeltaCoeffs)
        for (const label facei : snGradFaces_)
        {
            const vector n(Sf[facei]/magSf[facei]);
            const vector& d = dCentroid_[facei];
            const vector gf(0.5*(g[own[facei]] + g[nei[facei]]));
            c[facei] = (gf & n) - (gf & d)*cutDeltaCoeffs_[facei];
        }
        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& fvp = mesh_.boundary()[patchi];
            if (!fvp.coupled())
            {
                // A fixed value on a cut face is the value at its wet-face
                // centroid x'_f, the cell value is at its value point x_v:
                // the boundary condition differences them over the mesh
                // delta (cell centre to face centre), n.d_mesh. With
                // d = x'_f - x_v the gradient is
                //   (v_b - v_P - g_P.(d - (n.d) n))/(n.d)
                if (!vf.boundaryField()[patchi].fixesValue()) continue;
                const labelUList& fc = fvp.faceCells();
                const fvsPatchScalarField& tp = theta_.boundaryField()[patchi];
                const fvsPatchVectorField& Sfp = Sf.boundaryField()[patchi];
                const fvsPatchScalarField& magSfp = magSf.boundaryField()[patchi];
                const fvsPatchVectorField& Cfp = Cf.boundaryField()[patchi];
                const fvsPatchVectorField& wop =
                    wetFaceOffset_.boundaryField()[patchi];
                const scalarField& dcMesh = fvp.deltaCoeffs();
                const scalarField& vb = vc.boundaryField()[patchi];
                fvsPatchScalarField& cp = c.boundaryFieldRef()[patchi];
                forAll(fvp, i)
                {
                    const label o = fc[i];
                    if (tp[i] <= 0 || solidCells_.test(o)) continue;
                    if (tp[i] >= 1 && magSqr(off[o]) == 0) continue;
                    const vector n(Sfp[i]/magSfp[i]);
                    const vector d(Cfp[i] + wop[i] - C[o] - off[o]);
                    const scalar nd = max(n & d, 0.1/dcMesh[i]);
                    cp[i] =
                        (vb[i] - vc[o])*(1.0/nd - dcMesh[i])
                      - (g[o] & (d - (n & d)*n))/nd;
                }
                continue;
            }
            const labelUList& fc = fvp.faceCells();
            const vectorField gN(g.boundaryField()[patchi].patchNeighbourField());
            const fvsPatchVectorField& dp = dCentroid_.boundaryField()[patchi];
            const fvsPatchScalarField& dcp = cutDeltaCoeffs_.boundaryField()[patchi];
            const fvsPatchVectorField& Sfp = Sf.boundaryField()[patchi];
            const fvsPatchScalarField& magSfp = magSf.boundaryField()[patchi];
            fvsPatchScalarField& cp = c.boundaryFieldRef()[patchi];
            for (const label i : snGradPatchFaces_[patchi])
            {
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
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::cutCellGeometry::convectionSource
(
    const surfaceScalarField& phi,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const volScalarField* nu
) const
{
    tmp<fv::convectionScheme<Type>> tcs
    (
        fv::convectionScheme<Type>::New
        (
            mesh_, phi, mesh_.divScheme("div(" + phi.name() + ',' + vf.name() + ')')
        )
    );
    return fvc::surfaceIntegrate
    (
        convectionCorrection(phi, vf, tcs().interpolate(phi, vf)(), nu)
    );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::cutCellGeometry::convectionCorrection
(
    const surfaceScalarField& phi,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& vSchemeFace,
    const volScalarField* nu
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceField;

    tmp<SurfaceField> tcorr
    (
        new SurfaceField
        (
            IOobject
            (
                "ibmConvectionCorr(" + vf.name() + ')',
                vf.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<Type>(phi.dimensions()*vf.dimensions(), Zero)
        )
    );
    if (!ibmConvection_)
    {
        return tcorr;
    }
    SurfaceField& corr = tcorr.ref();

    const direction nCmpt = pTraits<Type>::nComponents;
    PtrList<volVectorField> g(nCmpt);
    for (direction cmpt = 0; cmpt < nCmpt; ++cmpt)
    {
        g.set(cmpt, lsqGrad(vf.component(cmpt), true).ptr());
    }

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const volVectorField& C = mesh_.C();
    const surfaceVectorField& Cf = mesh_.Cf();
    const vectorField& off = valueOffset_.primitiveField();

    // Wall cells with a wall function: the wall law across the wall (velocity)
    wallProfileData prof;
    if (nu && ibmConvectionWallLaw_)
    {
        wallProfile(vf, *nu, prof);
    }
    const bool law = prof.valid();

    // Value of cell c (value point x_c + off_c) extrapolated to the point x
    auto extrapolate = [&]
    (
        const Type& v, const label c, const vector& r
    )
    {
        Type e(v);
        for (direction cmpt = 0; cmpt < nCmpt; ++cmpt)
        {
            setComponent(e, cmpt) += g[cmpt][c] & r;
        }
        return e;
    };

    forAll(own, facei)
    {
        const label o = own[facei], n = nei[facei];
        if (theta_[facei] <= 0 || solidCells_.test(o) || solidCells_.test(n))
        {
            continue;
        }
        const label u = (phi[facei] >= 0) ? o : n;
        if (theta_[facei] >= 1 && magSqr(off[u]) == 0)
        {
            continue;
        }
        const vector r(Cf[facei] + wetFaceOffset_[facei] - C[u] - off[u]);
        Type e(extrapolate(vf[u], u, r));
        if (law && prof.w()[u] > 0)
        {
            e = wallProfileValue
            (
                e, vf[u], r, prof.n()[u], prof.Ub()[u], dWall_[u],
                prof.uTau()[u], prof.w()[u], (*nu)[u]
            );
        }
        corr[facei] = phi[facei]*(e - vSchemeFace[facei]);
    }

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];
        if (!fvp.coupled()) continue;
        const labelUList& fc = fvp.faceCells();
        const Field<Type> vN(vf.boundaryField()[patchi].patchNeighbourField());
        const vectorField offN
        (
            valueOffset_.boundaryField()[patchi].patchNeighbourField()
        );
        PtrList<vectorField> gN(nCmpt);
        for (direction cmpt = 0; cmpt < nCmpt; ++cmpt)
        {
            gN.set
            (
                cmpt,
                new vectorField(g[cmpt].boundaryField()[patchi].patchNeighbourField())
            );
        }
        const vectorField d(fvp.delta());
        const fvsPatchScalarField& tp = theta_.boundaryField()[patchi];
        const fvsPatchScalarField& php = phi.boundaryField()[patchi];
        const fvsPatchVectorField& Cfp = Cf.boundaryField()[patchi];
        const fvsPatchVectorField& wop = wetFaceOffset_.boundaryField()[patchi];
        const fvsPatchField<Type>& vsp = vSchemeFace.boundaryField()[patchi];
        fvsPatchField<Type>& cp = corr.boundaryFieldRef()[patchi];
        scalarField uTauN, wN, yN, nuN;
        vectorField nN, UbN;
        if (law)
        {
            uTauN = prof.uTau().boundaryField()[patchi].patchNeighbourField();
            wN = prof.w().boundaryField()[patchi].patchNeighbourField();
            yN = dWall_.boundaryField()[patchi].patchNeighbourField();
            nuN = nu->boundaryField()[patchi].patchNeighbourField();
            nN = prof.n().boundaryField()[patchi].patchNeighbourField();
            UbN = prof.Ub().boundaryField()[patchi].patchNeighbourField();
        }
        forAll(fvp, i)
        {
            const label c = fc[i];
            if (tp[i] <= 0 || solidCells_.test(c)) continue;
            const point xf(Cfp[i] + wop[i]);
            Type v;
            if (php[i] >= 0)
            {
                if (tp[i] >= 1 && magSqr(off[c]) == 0) continue;
                const vector r(xf - C[c] - off[c]);
                v = extrapolate(vf[c], c, r);
                if (law && prof.w()[c] > 0)
                {
                    v = wallProfileValue
                    (
                        v, vf[c], r, prof.n()[c], prof.Ub()[c], dWall_[c],
                        prof.uTau()[c], prof.w()[c], (*nu)[c]
                    );
                }
            }
            else
            {
                if (tp[i] >= 1 && magSqr(offN[i]) == 0) continue;
                const vector r(xf - (C[c] + d[i]) - offN[i]);
                v = vN[i];
                for (direction cmpt = 0; cmpt < nCmpt; ++cmpt)
                {
                    setComponent(v, cmpt) += gN[cmpt][i] & r;
                }
                if (law && wN[i] > 0)
                {
                    v = wallProfileValue
                    (
                        v, vN[i], r, nN[i], UbN[i], yN[i], uTauN[i], wN[i],
                        nuN[i]
                    );
                }
            }
            cp[i] = php[i]*(v - vsp[i]);
        }
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
void Foam::cutCellGeometry::groupSum(UList<Type>& f) const
{
    List<Type> s(slotGroup_.size(), Zero);
    forAll(groupSlot_, celli)
    {
        if (groupSlot_[celli] >= 0) s[groupSlot_[celli]] += f[celli];
    }

    groupSlotSum(s);

    forAll(groupSlot_, celli)
    {
        if (groupSlot_[celli] >= 0) f[celli] = s[groupSlot_[celli]];
    }
}


template<class Type>
void Foam::cutCellGeometry::groupSlotSum(List<Type>& s) const
{
    if (Pstream::parRun())
    {
        // Partial sums to the processor of the master, which adds them in
        // processor order, then the totals back
        PstreamBuffers toMaster(Pstream::commsTypes::nonBlocking);
        forAll(sendSlots_, proci)
        {
            if (sendSlots_[proci].size())
            {
                UOPstream os(proci, toMaster);
                os << List<Type>(UIndirectList<Type>(s, sendSlots_[proci]));
            }
        }
        toMaster.finishedSends();
        forAll(recvSlots_, proci)
        {
            if (recvSlots_[proci].size())
            {
                UIPstream is(proci, toMaster);
                const List<Type> v(is);
                forAll(v, k)
                {
                    s[recvSlots_[proci][k]] += v[k];
                }
            }
        }

        PstreamBuffers fromMaster(Pstream::commsTypes::nonBlocking);
        forAll(recvSlots_, proci)
        {
            if (recvSlots_[proci].size())
            {
                UOPstream os(proci, fromMaster);
                os << List<Type>(UIndirectList<Type>(s, recvSlots_[proci]));
            }
        }
        fromMaster.finishedSends();
        forAll(sendSlots_, proci)
        {
            if (sendSlots_[proci].size())
            {
                UIPstream is(proci, fromMaster);
                const List<Type> v(is);
                forAll(v, k)
                {
                    s[sendSlots_[proci][k]] = v[k];
                }
            }
        }
    }
}


template<class Type>
void Foam::cutCellGeometry::fillUncovered
(
    GeometricField<Type, fvPatchField, volMesh>& f
) const
{
    if (!moving_ || solidOld_.size() != mesh_.nCells())
    {
        return;
    }

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    Field<Type>& fI = f.primitiveFieldRef();

    // Known: fluid before and now; to fill: fluid now, solid before
    boolList known(mesh_.nCells(), false);
    label nFill = 0;
    forAll(known, celli)
    {
        if (solidCells_.test(celli)) continue;
        known[celli] = !solidOld_.test(celli);
        if (!known[celli]) ++nFill;
    }
    if (returnReduce(nFill, sumOp<label>()) == 0)
    {
        return;
    }

    for (label layer = 0; layer < 10; ++layer)
    {
        f.correctBoundaryConditions();
        boolList nbrKnown;
        syncTools::swapBoundaryCellList(mesh_, known, nbrKnown);

        Field<Type> sum(fI.size(), Zero);
        labelList n(fI.size(), 0);
        forAll(own, facei)
        {
            if (theta_[facei] <= 0) continue;
            const label o = own[facei], nb = nei[facei];
            if (known[o] && !known[nb]) { sum[nb] += fI[o]; ++n[nb]; }
            if (known[nb] && !known[o]) { sum[o] += fI[nb]; ++n[o]; }
        }
        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& fvp = mesh_.boundary()[patchi];
            if (!fvp.coupled()) continue;
            const labelUList& fc = fvp.faceCells();
            const fvsPatchScalarField& tp = theta_.boundaryField()[patchi];
            const Field<Type> fN(f.boundaryField()[patchi].patchNeighbourField());
            forAll(fvp, i)
            {
                const label bFacei = fvp.start() + i - mesh_.nInternalFaces();
                if (tp[i] > 0 && nbrKnown[bFacei] && !known[fc[i]])
                {
                    sum[fc[i]] += fN[i];
                    ++n[fc[i]];
                }
            }
        }

        label nLeft = 0;
        forAll(known, celli)
        {
            if (solidCells_.test(celli) || known[celli]) continue;
            if (n[celli] > 0)
            {
                fI[celli] = sum[celli]/n[celli];
                known[celli] = true;
            }
            else
            {
                ++nLeft;
            }
        }
        if (returnReduce(nLeft, sumOp<label>()) == 0)
        {
            break;
        }
    }
    f.correctBoundaryConditions();
}


template<class Type>
Foam::tmp<Foam::scalarField>
Foam::cutCellGeometry::coupledDiag(const fvMatrix<Type>& eqn) const
{
    tmp<scalarField> td(new scalarField(eqn.diag()));
    scalarField& d = td.ref();
    forAll(mesh_.boundary(), patchi)
    {
        if (!mesh_.boundary()[patchi].coupled()) continue;
        const labelUList& fc = mesh_.boundary()[patchi].faceCells();
        const Field<Type>& ic = eqn.internalCoeffs()[patchi];
        forAll(fc, i)
        {
            d[fc[i]] += cmptAv(ic[i]);
        }
    }
    return td;
}


template<class Type>
void Foam::cutCellGeometry::linkCells
(
    fvMatrix<Type>& eqn,
    const bool groupUniform
) const
{
    if (mergeAlpha_ <= 0)
    {
        return;
    }

    if (groupUniform)
    {
        // The residual A psi - b of every group cell at the current psi,
        // less its volume share of the residual of its group, added to the
        // source: a psi uniform in the groups that satisfies the merged
        // equations then satisfies every linked row (the sum over a group is
        // unchanged and the penalty carries nothing)
        const scalarField& V = mesh_.V();
        const Field<Type> r((eqn & eqn.psi())().primitiveField()*V);
        Field<Type> rGroup(r.size(), Zero);
        scalarField vol(r.size(), Zero);
        forAll(linkGroup_, celli)
        {
            if (linkGroup_[celli] < 0) continue;
            rGroup[celli] = r[celli];
            vol[celli] = V[celli];
        }
        groupSum(rGroup);
        groupSum(vol);
        Field<Type>& source = eqn.source();
        forAll(linkGroup_, celli)
        {
            if (linkGroup_[celli] < 0) continue;
            source[celli] += r[celli] - V[celli]/vol[celli]*rGroup[celli];
        }
    }

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    scalarField& diag = eqn.diag();
    const bool symmetric = eqn.symmetric();
    scalarField& upper = eqn.upper();
    scalarField* lowerPtr = symmetric ? nullptr : &eqn.lower();

    // The conductances from the diagonal before linking (independent of the
    // order of the links, and of the decomposition), and that of the cell
    // across the coupled faces
    const scalarField diag0(coupledDiag(eqn));
    scalarField nbrDiag0;
    syncTools::swapBoundaryCellList(mesh_, diag0, nbrDiag0);

    // A conductance with the sign of the diagonal (positive for -laplacian,
    // negative for the pressure laplacian)
    auto conductance = [&](const scalar dO, const scalar dN, scalar& sgn)
    {
        sgn = (dO + dN >= 0) ? 1 : -1;
        return mergePenalty_*max(mag(dO), mag(dN));
    };

    forAll(linkFace_, celli)
    {
        const label facei = linkFace_[celli];
        if (facei < 0 || facei >= mesh_.nInternalFaces()) continue;
        const label other = (own[facei] == celli) ? nei[facei] : own[facei];
        scalar sgn;
        const scalar c = conductance(diag0[celli], diag0[other], sgn);
        diag[celli] += sgn*c;
        diag[other] += sgn*c;
        upper[facei] -= sgn*c;
        if (lowerPtr)
        {
            (*lowerPtr)[facei] -= sgn*c;
        }
    }

    // Links over coupled faces: each side adds the conductance to its own
    // row (internalCoeffs to the diagonal, -boundaryCoeffs is the
    // off-diagonal coefficient of the cell across), as the coupled patches
    // of gaussLaplacianScheme
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];
        if (!fvp.coupled()) continue;
        const labelUList& fc = fvp.faceCells();
        Field<Type>& ic = eqn.internalCoeffs()[patchi];
        Field<Type>& bc = eqn.boundaryCoeffs()[patchi];
        forAll(fvp, i)
        {
            const label bFacei = fvp.start() + i - mesh_.nInternalFaces();
            if (!boundaryLink_[bFacei]) continue;
            scalar sgn;
            const scalar c =
                conductance(diag0[fc[i]], nbrDiag0[bFacei], sgn);
            ic[i] += sgn*c*pTraits<Type>::one;
            bc[i] += sgn*c*pTraits<Type>::one;
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
    forAll(linkGroup_, celli)
    {
        if (linkGroup_[celli] < 0) continue;
        sum[celli] = alpha[celli]*V[celli]*f[celli];
        vol[celli] = alpha[celli]*V[celli];
    }
    groupSum(sum);
    groupSum(vol);
    forAll(linkGroup_, celli)
    {
        if (linkGroup_[celli] >= 0 && vol[celli] > 0)
        {
            f[celli] = sum[celli]/vol[celli];
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
    forAll(linkGroup_, celli)
    {
        if (linkGroup_[celli] < 0) continue;
        sum[celli] = V[celli]*fI[celli];
        vol[celli] = V[celli];
    }
    groupSum(sum);
    groupSum(vol);
    forAll(linkGroup_, celli)
    {
        if (linkGroup_[celli] >= 0)
        {
            fI[celli] = sum[celli]/vol[celli];
        }
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
    forAll(linkGroup_, celli)
    {
        if (linkGroup_[celli] < 0) continue;
        sum[celli] = V[celli]*AI[celli];
        vol[celli] = V[celli];
        const label facei = linkFace_[celli];
        if (facei >= 0 && facei < mesh_.nInternalFaces())
        {
            // The coefficients of the internal face drop out of the merged
            // diagonal (A_i holds -sum of the off-diagonals of its row)
            sum[celli] += upper[facei] + lower[facei];
        }
    }
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];
        if (!fvp.coupled()) continue;
        const labelUList& fc = fvp.faceCells();
        const Field<Type>& bc = eqn.boundaryCoeffs()[patchi];
        forAll(fvp, i)
        {
            if (boundaryLink_[fvp.start() + i - mesh_.nInternalFaces()])
            {
                // Each side: the off-diagonal of its own row, -boundaryCoeffs
                sum[fc[i]] -= cmptAv(bc[i]);
            }
        }
    }
    groupSum(sum);
    groupSum(vol);
    forAll(linkGroup_, celli)
    {
        if (linkGroup_[celli] >= 0)
        {
            AI[celli] = sum[celli]/vol[celli];
        }
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
    const GeometricField<Type, fvPatchField, volMesh>& psi = eqn.psi();
    const Field<Type>& psiI = psi.primitiveField();
    Field<Type>& HI = H.primitiveFieldRef();
    Field<Type> sum(HI.size(), Zero);
    scalarField vol(HI.size(), Zero);
    forAll(linkGroup_, celli)
    {
        if (linkGroup_[celli] < 0) continue;
        sum[celli] = V[celli]*HI[celli];
        vol[celli] = V[celli];
        const label facei = linkFace_[celli];
        if (facei >= 0 && facei < mesh_.nInternalFaces())
        {
            // H_own holds -upper psi_nei, H_nei -lower psi_own: inside the
            // merged cell, not a coupling to a neighbour
            sum[celli] +=
                upper[facei]*psiI[nei[facei]] + lower[facei]*psiI[own[facei]];
        }
    }
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];
        if (!fvp.coupled()) continue;
        const labelUList& fc = fvp.faceCells();
        const Field<Type>& bc = eqn.boundaryCoeffs()[patchi];
        const Field<Type> psiN(psi.boundaryField()[patchi].patchNeighbourField());
        forAll(fvp, i)
        {
            if (boundaryLink_[fvp.start() + i - mesh_.nInternalFaces()])
            {
                // H holds +boundaryCoeffs psi_N on a coupled face
                sum[fc[i]] -= cmptMultiply(bc[i], psiN[i]);
            }
        }
    }
    groupSum(sum);
    groupSum(vol);
    forAll(linkGroup_, celli)
    {
        if (linkGroup_[celli] >= 0)
        {
            HI[celli] = sum[celli]/vol[celli];
        }
    }
    H.correctBoundaryConditions();
}


// ************************************************************************* //
