/*---------------------------------------------------------------------------*\
    cutCellIBM: wall shear of the immersed no-slip walls in turbulent flow
\*---------------------------------------------------------------------------*/

#include "cutCellGeometry.H"
#include "DynamicField.H"
#include "cyclicPolyPatch.H"
#include "globalIndex.H"
#include "PstreamBuffers.H"
#include "treeDataCell.H"

namespace
{

// Law of the wall of the Spalart-Allmaras model in an equilibrium boundary
// layer (Allmaras, Johnson & Spalart, ICCFD7-1902, 2012): u+(y+) in closed
// form, u+ = y+ at the wall and ln(y+)/0.41 + 5.0334 in the log layer
Foam::scalar saUPlusOfYPlus(const Foam::scalar yp)
{
    const Foam::scalar
        B = 5.0333908790505579,
        a1 = 8.148221580024245, a2 = -6.9287093849022945,
        b1 = 7.4600876082527945, b2 = 7.468145790401841,
        c1 = 2.5496773539754747, c2 = 1.3301651588535228,
        c3 = 3.599459109332379, c4 = 3.6397531868684494;
    return
        B
      + c1*::log(Foam::sqr(yp + a1) + Foam::sqr(b1))
      - c2*::log(Foam::sqr(yp + a2) + Foam::sqr(b2))
      - c3*::atan2(b1, yp + a1)
      - c4*::atan2(b2, yp + a2);
}

// u+ from Re_y = u+ y+ with the SA law: y+ u+(y+) increases with y+, and
// u+ <= y+ puts the root between sqrt(Re_y) and max(Re_y, 1)
Foam::scalar saUPlus(const Foam::scalar Rey)
{
    using namespace Foam;
    if (Rey <= 0) return 0;

    scalar lo = ::sqrt(Rey), hi = max(Rey, scalar(1));
    for (label it = 0; it < 200 && hi - lo > 1e-12*hi; ++it)
    {
        const scalar mid = 0.5*(lo + hi);
        if (mid*saUPlusOfYPlus(mid) > Rey) hi = mid; else lo = mid;
    }
    return Rey/(0.5*(lo + hi));
}


} // End anonymous namespace


Foam::scalar Foam::cutCellGeometry::wallFunctionUTau
(
    const vector& Us,
    const volScalarField& nu,
    const label celli,
    const scalar y
) const
{
    const label bodyi = cellBody_[celli];
    if
    (
        bodyi < 0
     || bodies_[bodyi].wall() != ibmBody::NOSLIP
     || bodies_[bodyi].wallFunction() == ibmBody::NOWALLFUNCTION
    )
    {
        return 0;
    }

    // Tangential slip of the sample along the wall of the cell
    const vector n(Sw_[celli]/max(Awall_[celli], VSMALL));
    vector slip(Us);
    if (wallVelocity_)
    {
        slip -= UbPtr_()[celli];
    }
    const scalar magUt = mag(slip - (slip & n)*n);
    const scalar up = saUPlus(magUt*y/nu[celli]);
    return (up > SMALL) ? magUt/up : 0;
}


void Foam::cutCellGeometry::calcWallSamples() const
{
    wallSamplesDone_ = true;
    const label nCells = mesh_.nCells();
    wallCellSize_.setSize(nCells);
    wallCellSize_ = 0;
    sampleCell_.setSize(nCells);
    sampleCell_ = -1;
    sampleDist_.setSize(nCells);
    sampleDist_ = 0;
    sampleSlot_.setSize(nCells);
    sampleSlot_ = -1;
    nSampleSlots_ = 0;
    sampleSendCells_.setSize(Pstream::nProcs());
    sampleRecvSlots_.setSize(Pstream::nProcs());
    forAll(sampleSendCells_, proci)
    {
        sampleSendCells_[proci].clear();
        sampleRecvSlots_[proci].clear();
    }

    // Extent of the (whole) wall cell along the wall normal, V over its
    // projected area: the resolution of the wall-normal gradients
    const vectorField& Sf = mesh_.faceAreas();
    const scalarField& V = mesh_.V();
    for (const label celli : wallCells_)
    {
        const vector n(Sw_[celli]/max(Awall_[celli], VSMALL));
        scalar projected = 0;
        for (const label facei : mesh_.cells()[celli])
        {
            projected += mag(Sf[facei] & n);
        }
        wallCellSize_[celli] = 2*V[celli]/max(projected, VSMALL);
    }

    if (wallSampleDistance_ <= 0)
    {
        return;
    }
    forAll(bodies_, bodyi)
    {
        if (!bodies_[bodyi].hasLevelSet())
        {
            FatalErrorInFunction
                << "Body " << bodies_[bodyi].name() << " has no level set: "
                << "the wall-law samples need the distance from the immersed "
                << "walls (level-set bodies), or wallSampleDistance 0"
                << exit(FatalError);
        }
    }

    // The point searches below run on some processors only: build the
    // search structures, whose construction communicates, on all of them
    (void)mesh_.tetBasePtIs();
    (void)mesh_.cellTree();

    // Sample point: wallSampleDistance h from the wall along the normal of
    // the wall cell, through its fluid centroid
    const pointField xc(fluidCentroids());
    DynamicList<label> cells;
    for (const label celli : wallCells_)
    {
        if (wallSampleDistance_*wallCellSize_[celli] > dWall_[celli])
        {
            cells.append(celli);
        }
    }
    scalarField dc;
    wallDistance(pointField(xc, cells), dc);

    // The cell holding it: its value sits at its fluid centroid, or at the
    // merged centroid for a linked cell. A sample in the wall cell's own
    // group (or in the solid) leaves the wall cell sampling itself. The
    // group key is global (the group, or the global cell index), the same
    // on every processor.
    const globalIndex cellIndex(nCells);
    auto group = [&](const label c)
    {
        return (linkGroup_[c] >= 0) ? linkGroup_[c] : -1 - cellIndex.toGlobal(c);
    };
    // Separations of the translational cyclics: a sample beyond one is
    // looked for in the periodic image
    DynamicList<vector> separations;
    for (const polyPatch& pp : mesh_.boundaryMesh())
    {
        const cyclicPolyPatch* cpp = isA<cyclicPolyPatch>(pp);
        if (cpp && cpp->owner() && cpp->separated() && cpp->separation().size())
        {
            separations.append(cpp->separation()[0]);
        }
    }

    // The cell holding a point near a wall cell: searched in its face
    // neighbours, layer by layer, before the whole mesh (the octree of
    // findCell has a few leaves holding all the cells of a mesh refined
    // around a small body, and is then linear in their number)
    const labelListList& cellCells = mesh_.cellCells();
    auto locate = [&](const label celli, const point& p)
    {
        DynamicList<label> front(1, celli);
        DynamicList<label> next;
        labelHashSet visited;
        visited.insert(celli);
        for (label layer = 0; layer < 4; ++layer)
        {
            for (const label c : front)
            {
                if (mesh_.pointInCell(p, c)) return c;
            }
            next.clear();
            for (const label c : front)
            {
                for (const label nb : cellCells[c])
                {
                    if (visited.insert(nb)) next.append(nb);
                }
            }
            front = next;
        }
        return mesh_.findCell(p);
    };

    // The cell of this processor holding a point, also in the periodic
    // images
    auto findLocal = [&](const point& xs)
    {
        label s = mesh_.findCell(xs);
        for (label k = 0; s < 0 && k < separations.size(); ++k)
        {
            s = mesh_.findCell(xs + separations[k]);
            if (s < 0) s = mesh_.findCell(xs - separations[k]);
        }
        return s;
    };

    label nMissed = 0;
    DynamicList<label> found;
    DynamicField<point> pts;
    DynamicList<label> queryCells;
    DynamicField<point> queryPts;
    forAll(cells, i)
    {
        const label celli = cells[i];
        const vector n(Sw_[celli]/max(Awall_[celli], VSMALL));
        const point xs
        (
            xc[celli] - (wallSampleDistance_*wallCellSize_[celli] - dc[i])*n
        );
        label s = locate(celli, xs);
        for (label k = 0; s < 0 && k < separations.size(); ++k)
        {
            s = mesh_.findCell(xs + separations[k]);
            if (s < 0) s = mesh_.findCell(xs - separations[k]);
        }
        if (s < 0 && Pstream::parRun())
        {
            // On another processor (or outside the domain)
            queryCells.append(celli);
            queryPts.append(xs);
        }
        else if (s < 0 || alpha_[s] <= 0)
        {
            ++nMissed;
        }
        else if (group(s) != group(celli))
        {
            sampleCell_[celli] = s;
            found.append(celli);
            pts.append(xc[s]);
        }
    }
    scalarField ds;
    wallDistance(pts, ds);
    forAll(found, i)
    {
        const label s = sampleCell_[found[i]];
        sampleDist_[found[i]] = (linkWallDist_[s] > 0) ? linkWallDist_[s] : ds[i];
    }

    // Samples on other processors: every processor looks for the points the
    // others did not find, the lowest processor holding a point answers
    // with its cell, whether it is fluid, its group key and the distance of
    // its value point from the wall; the values are exchanged in wallLaw
    label nRemote = 0;
    if (Pstream::parRun())
    {
        const label nProcs = Pstream::nProcs();
        const label myProc = Pstream::myProcNo();
        List<pointField> allQueries(nProcs);
        {
            PstreamBuffers toAll(Pstream::commsTypes::nonBlocking);
            if (queryPts.size())
            {
                for (label proci = 0; proci < nProcs; ++proci)
                {
                    if (proci == myProc) continue;
                    UOPstream os(proci, toAll);
                    os << pointField(queryPts);
                }
            }
            labelList recvSizes;
            toAll.finishedSends(recvSizes);
            forAll(recvSizes, proci)
            {
                if (proci == myProc || !recvSizes[proci]) continue;
                UIPstream is(proci, toAll);
                is >> allQueries[proci];
            }
        }

        // (the wall distance of all the hits in one call: the level set
        // of a surface body may communicate)
        List<labelList> cellOf(nProcs);
        List<labelList> keyOf(nProcs);
        List<boolList> fluid(nProcs);
        DynamicList<labelPair> hits;
        DynamicField<point> hitPts;
        for (label proci = 0; proci < nProcs; ++proci)
        {
            if (proci == myProc) continue;
            const pointField& q = allQueries[proci];
            cellOf[proci].setSize(q.size(), -1);
            keyOf[proci].setSize(q.size(), 0);
            fluid[proci].setSize(q.size(), false);
            forAll(q, k)
            {
                const label s = findLocal(q[k]);
                cellOf[proci][k] = s;
                if (s < 0) continue;
                keyOf[proci][k] = group(s);
                fluid[proci][k] = alpha_[s] > 0;
                hits.append(labelPair(proci, k));
                hitPts.append(xc[s]);
            }
        }
        scalarField dHit;
        wallDistance(hitPts, dHit);
        List<scalarField> distOf(nProcs);
        forAll(distOf, proci)
        {
            distOf[proci].setSize(cellOf[proci].size(), 0);
        }
        forAll(hits, j)
        {
            const label proci = hits[j].first(), k = hits[j].second();
            const label s = cellOf[proci][k];
            distOf[proci][k] = (linkWallDist_[s] > 0) ? linkWallDist_[s] : dHit[j];
        }

        PstreamBuffers answers(Pstream::commsTypes::nonBlocking);
        for (label proci = 0; proci < nProcs; ++proci)
        {
            if (proci == myProc || allQueries[proci].empty()) continue;
            UOPstream os(proci, answers);
            os << cellOf[proci] << keyOf[proci] << fluid[proci] << distOf[proci];
        }
        labelList recvSizes;
        answers.finishedSends(recvSizes);

        // The lowest processor that holds the point
        labelList fromProc(queryCells.size(), -1);
        labelList remoteCell(queryCells.size(), -1);
        for (label proci = 0; proci < nProcs; ++proci)
        {
            if (proci == myProc || !recvSizes[proci]) continue;
            UIPstream is(proci, answers);
            const labelList cellOf(is);
            const labelList keyOf(is);
            const boolList fluid(is);
            const scalarField distOf(is);
            forAll(queryCells, k)
            {
                if (fromProc[k] >= 0 || cellOf[k] < 0) continue;
                fromProc[k] = proci;
                const label celli = queryCells[k];
                if (!fluid[k])
                {
                    remoteCell[k] = -2;
                }
                else if (keyOf[k] != group(celli))
                {
                    remoteCell[k] = cellOf[k];
                    sampleDist_[celli] = distOf[k];
                }
            }
        }

        List<DynamicList<label>> request(nProcs);
        List<DynamicList<label>> slots(nProcs);
        forAll(queryCells, k)
        {
            if (fromProc[k] < 0 || remoteCell[k] == -2)
            {
                ++nMissed;
            }
            else if (remoteCell[k] >= 0)
            {
                sampleSlot_[queryCells[k]] = nSampleSlots_;
                request[fromProc[k]].append(remoteCell[k]);
                slots[fromProc[k]].append(nSampleSlots_);
                ++nSampleSlots_;
            }
        }
        nRemote = nSampleSlots_;

        PstreamBuffers requests(Pstream::commsTypes::nonBlocking);
        forAll(request, proci)
        {
            if (request[proci].size())
            {
                UOPstream os(proci, requests);
                os << labelList(request[proci]);
            }
            sampleRecvSlots_[proci].transfer(slots[proci]);
        }
        requests.finishedSends(recvSizes);
        forAll(recvSizes, proci)
        {
            if (proci == myProc || !recvSizes[proci]) continue;
            UIPstream is(proci, requests);
            sampleSendCells_[proci] = labelList(is);
        }
    }

    Info<< "    wall-law samples at " << wallSampleDistance_ << " h: "
        << returnReduce(found.size() + nRemote, sumOp<label>())
        << " wall cells";
    reduce(nRemote, sumOp<label>());
    if (nRemote > 0)
    {
        Info<< " (" << nRemote << " on another processor)";
    }
    reduce(nMissed, sumOp<label>());
    if (nMissed > 0)
    {
        Info<< ", " << nMissed << " sample themselves (sample in the solid "
            << "or outside the domain)";
    }
    Info<< endl;
}


void Foam::cutCellGeometry::wallLaw
(
    const volVectorField& U,
    const volScalarField& nu,
    scalarField& uTau,
    scalarField& weight
) const
{
    uTau.setSize(mesh_.nCells());
    uTau = 0;
    weight.setSize(mesh_.nCells());
    weight = 0;

    const scalar lnLo = ::log(hPlusResolved_);
    const scalar lnHi = ::log(max(hPlusModelled_, hPlusResolved_*(1 + SMALL)));

    if (!wallSamplesDone_)
    {
        calcWallSamples();
    }

    // Velocities of the samples on other processors
    vectorField Uremote(nSampleSlots_, Zero);
    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
        forAll(sampleSendCells_, proci)
        {
            if (sampleSendCells_[proci].size())
            {
                UOPstream os(proci, pBufs);
                os << vectorField
                (
                    UIndirectList<vector>(U.primitiveField(), sampleSendCells_[proci])
                );
            }
        }
        pBufs.finishedSends();
        forAll(sampleRecvSlots_, proci)
        {
            if (sampleRecvSlots_[proci].size())
            {
                UIPstream is(proci, pBufs);
                const vectorField v(is);
                UIndirectList<vector>(Uremote, sampleRecvSlots_[proci]) = v;
            }
        }
    }

    for (const label celli : wallCells_)
    {
        const label s = sampleCell_[celli];
        const label slot = sampleSlot_[celli];
        const bool sampled = (s >= 0 || slot >= 0);
        const scalar utCell =
            wallFunctionUTau(U[celli], nu, celli, dWall_[celli]);
        const scalar utSample =
            !sampled
          ? utCell
          : wallFunctionUTau
            (
                (s >= 0) ? U[s] : Uremote[slot], nu, celli, sampleDist_[celli]
            );
        if (utSample <= 0) continue;

        const scalar x = min
        (
            max
            (
                (::log(utSample*wallCellSize_[celli]/nu[celli]) - lnLo)
               /(lnHi - lnLo),
                0
            ),
            1
        );
        const scalar w = sqr(x)*(3 - 2*x);
        weight[celli] = w;
        uTau[celli] =
            (!sampled || w <= 0)
          ? utCell
          : ::sqrt((1 - w)*sqr(utCell) + w*sqr(utSample));
    }
}


void Foam::cutCellGeometry::wallProfile
(
    const volVectorField& U,
    const volScalarField& nu,
    wallProfileData& d
) const
{
    scalarField uTau, w;
    wallLaw(U, nu, uTau, w);
    if (gMax(w) <= 0)
    {
        return;
    }

    auto field = [&](const word& name, const dimensionSet& dims)
    {
        return new volScalarField
        (
            IOobject(name, mesh_.time().timeName(), mesh_),
            mesh_,
            dimensionedScalar(dims, Zero),
            "zeroGradient"
        );
    };
    d.uTau.reset(field("ibmProfileUTau", dimVelocity));
    d.w.reset(field("ibmProfileWeight", dimless));
    d.uTau->primitiveFieldRef() = uTau;
    d.w->primitiveFieldRef() = w;
    d.uTau->correctBoundaryConditions();
    d.w->correctBoundaryConditions();

    d.n.reset
    (
        new volVectorField
        (
            IOobject("ibmProfileNormal", mesh_.time().timeName(), mesh_),
            mesh_,
            dimensionedVector(dimless, Zero),
            "zeroGradient"
        )
    );
    for (const label celli : wallCells_)
    {
        d.n->primitiveFieldRef()[celli] =
            Sw_[celli]/max(Awall_[celli], VSMALL);
    }
    d.n->correctBoundaryConditions();

    d.Ub.reset
    (
        new volVectorField
        (
            IOobject("ibmProfileUb", mesh_.time().timeName(), mesh_),
            mesh_,
            dimensionedVector(dimVelocity, Zero),
            "zeroGradient"
        )
    );
    if (wallVelocity_)
    {
        d.Ub->primitiveFieldRef() = UbPtr_().primitiveField();
    }
    d.Ub->correctBoundaryConditions();
}


Foam::vector Foam::cutCellGeometry::wallProfileValue
(
    const vector& e,
    const vector& v,
    const vector& r,
    const vector& n,
    const vector& Ub,
    const scalar y,
    const scalar uTau,
    const scalar w,
    const scalar nu
) const
{
    // The wall law along the wall normal from the value point to the
    // wet-face centroid (n out of the fluid: the distance falls along n)
    const scalar upV = saUPlusOfYPlus(y*uTau/nu);
    if (upV <= SMALL)
    {
        return e;
    }
    const scalar yf = max(y - (r & n), scalar(0));
    const scalar ratio = saUPlusOfYPlus(yf*uTau/nu)/upV;

    const vector slip(v - Ub);
    const vector vt(slip - (slip & n)*n);
    const vector eLaw(v + (ratio - 1)*vt + ((e - v) & n)*n);
    return (1 - w)*e + w*eLaw;
}


Foam::tmp<Foam::volScalarField> Foam::cutCellGeometry::wallViscosity
(
    const volVectorField& U,
    const volScalarField& nu
) const
{
    tmp<volScalarField> tnuW(new volScalarField("ibmWallViscosity", nu));
    scalarField& nuW = tnuW.ref().primitiveFieldRef();

    scalarField uTau, w;
    wallLaw(U, nu, uTau, w);

    for (const label celli : wallCells_)
    {
        if (uTau[celli] <= 0) continue;

        // nu_w |u_t|/dWall = u_tau^2 on the cell's own wall. With its own
        // value the law gives nu_w = nu d+/u+(d+) >= nu (u+ <= y+); with a
        // sample further out the wall cell may be faster than the linear
        // law (no floor nu, which would add to the wall shear of the wall
        // law), and nu_w is unbounded as the wall cell's own slip vanishes:
        // the implicit wall drag then holds the cell at rest, a second
        // steady state. On the wall law nu_w stays below nu + kappa u_tau d
        // (u+ >= ln(1 + kappa d+)/kappa), which caps it with d at least the
        // cell size
        vector slip(U[celli]);
        if (wallVelocity_)
        {
            slip -= UbPtr_()[celli];
        }
        const vector n(Sw_[celli]/max(Awall_[celli], VSMALL));
        const scalar magUt = max(mag(slip - (slip & n)*n), VSMALL);
        const scalar nuLaw = sqr(uTau[celli])*dWall_[celli]/magUt;
        nuW[celli] =
            min
            (
                nuLaw,
                nu[celli]
              + 0.41*uTau[celli]*max(wallCellSize_[celli], dWall_[celli])
            );
    }

    // Diagnostics: the kinematic wall shear the wall cells apply,
    // nu_w |u_t|/dWall (per cell the tangential part of the traction that
    // forces() sums), and the weight of the wall law. Registered fields
    // updated on every call and written with the solution, so that they
    // exist at the time where residualControl ends the run (written one
    // step after the converged iteration, which is no write step)
    if (!wallShearPtr_.valid() || wallShearPtr_().size() != mesh_.nCells())
    {
        // (check the old ones out before registering the new ones)
        wallShearPtr_.clear();
        wallLawWeightPtr_.clear();
        wallShearPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "ibmWallShear", mesh_.time().timeName(), mesh_,
                    IOobject::NO_READ, IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(sqr(dimVelocity), Zero),
                "zeroGradient"
            )
        );
        wallLawWeightPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "ibmWallLawWeight", mesh_.time().timeName(), mesh_,
                    IOobject::NO_READ, IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimless, Zero),
                "zeroGradient"
            )
        );
    }
    {
        scalarField& tauW = wallShearPtr_().primitiveFieldRef();
        scalarField& weight = wallLawWeightPtr_().primitiveFieldRef();
        tauW = 0;
        weight = 0;
        for (const label celli : wallCells_)
        {
            vector slip(U[celli]);
            if (wallVelocity_)
            {
                slip -= UbPtr_()[celli];
            }
            const vector n(Sw_[celli]/max(Awall_[celli], VSMALL));
            tauW[celli] =
                nuW[celli]*mag(slip - (slip & n)*n)/max(dWall_[celli], VSMALL);
            weight[celli] = w[celli];
        }
        wallShearPtr_().correctBoundaryConditions();
        wallLawWeightPtr_().correctBoundaryConditions();
    }

    tnuW.ref().correctBoundaryConditions();
    return tnuW;
}
