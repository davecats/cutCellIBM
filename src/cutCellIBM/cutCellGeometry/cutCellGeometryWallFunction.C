/*---------------------------------------------------------------------------*\
    cutCellIBM: wall shear of the immersed no-slip walls in turbulent flow
\*---------------------------------------------------------------------------*/

#include "cutCellGeometry.H"
#include "DynamicField.H"
#include "Map.H"
#include "ListOps.H"
#include "cyclicPolyPatch.H"

namespace
{

// Spalding's law, y+ = u+ + (exp(k u+) - 1 - k u+ - (k u+)^2/2 - (k u+)^3/6)/E,
// written with Re_y = u+ y+ = |u| y/nu as g(u+) = u+ + f(u+) - Re_y/u+ = 0.
// g is increasing, negative for small u+ and positive at u+ = sqrt(Re_y) (the
// viscous-sublayer root), so bisection on (0, sqrt(Re_y)] always converges.
Foam::scalar spaldingUPlus(const Foam::scalar Rey)
{
    using namespace Foam;
    const scalar kappa = 0.41, E = 9.8;
    if (Rey <= 0) return 0;

    auto g = [&](const scalar up)
    {
        const scalar ku = kappa*up;
        if (ku > 50) return GREAT;   // exp would overflow: far beyond the root
        return up + (::exp(ku) - 1 - ku - 0.5*sqr(ku) - pow3(ku)/6)/E - Rey/up;
    };

    scalar lo = 0, hi = ::sqrt(Rey);
    if (g(hi) <= 0) return hi;
    for (label it = 0; it < 100 && hi - lo > 1e-12*hi; ++it)
    {
        const scalar mid = 0.5*(lo + hi);
        if (g(mid) > 0) hi = mid; else lo = mid;
    }
    return 0.5*(lo + hi);
}


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


// u_tau from the SA law averaged over the wall distances y[0..n) (equal
// volumes): |u_t| = u_tau <u+(u_tau y/nu)>. The left side increases with
// u_tau; u+ <= y+ makes the viscous-sublayer value sqrt(|u_t| nu/<y>) a
// lower bound, and the upper bound is found by doubling
Foam::scalar saUTauAverage
(
    const Foam::scalar magUt,
    const Foam::scalar nu,
    const Foam::scalar* y,
    const Foam::label n
)
{
    using namespace Foam;
    scalar yMean = 0;
    for (label k = 0; k < n; ++k) yMean += y[k];
    yMean /= n;
    if (magUt <= 0 || yMean <= 0) return 0;

    auto f = [&](const scalar ut)
    {
        scalar s = 0;
        for (label k = 0; k < n; ++k) s += saUPlusOfYPlus(ut*y[k]/nu);
        return ut*s/n - magUt;
    };

    scalar lo = ::sqrt(magUt*nu/yMean), hi = 2*lo;
    for (label it = 0; it < 60 && f(hi) < 0; ++it) hi *= 2;
    for (label it = 0; it < 100 && hi - lo > 1e-12*hi; ++it)
    {
        const scalar mid = 0.5*(lo + hi);
        if (f(mid) > 0) hi = mid; else lo = mid;
    }
    return 0.5*(lo + hi);
}

} // End anonymous namespace


Foam::scalar Foam::cutCellGeometry::wallFunctionUTau
(
    const volVectorField& U,
    const volScalarField& nu,
    const label celli,
    const label samplei,
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

    // Tangential slip of the sample along the wall of the cell (a linked
    // cell has the value of the merged cell and the distance of its
    // centroid)
    const vector n(Sw_[celli]/max(Awall_[celli], VSMALL));
    vector slip(U[samplei]);
    if (moving_)
    {
        slip -= UbPtr_()[celli];
    }
    const scalar magUt = mag(slip - (slip & n)*n);
    const scalar Rey = magUt*y/nu[celli];
    const scalar up =
        (bodies_[bodyi].wallFunction() == ibmBody::SPALARTALLMARAS)
      ? saUPlus(Rey)
      : spaldingUPlus(Rey);
    return (up > SMALL) ? magUt/up : 0;
}


void Foam::cutCellGeometry::calcWallSamples()
{
    const label nCells = mesh_.nCells();
    wallCellSize_.setSize(nCells);
    wallCellSize_ = 0;
    sampleCell_.setSize(nCells);
    sampleCell_ = -1;
    sampleDist_.setSize(nCells);
    sampleDist_ = 0;

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
            WarningInFunction
                << "wallSampleDistance needs level-set bodies (the distance "
                << "of the sample from the wall): the wall cells sample "
                << "themselves" << endl;
            return;
        }
    }

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
    // group (or in the solid, or on another processor) leaves the wall
    // cell sampling itself.
    auto group = [&](const label c)
    {
        return (linkMaster_[c] >= 0) ? linkMaster_[c] : c;
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

    label nMissed = 0;
    DynamicList<label> found;
    DynamicField<point> pts;
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
        if (s < 0 || alpha_[s] <= 0)
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

    Info<< "    wall-law samples at " << wallSampleDistance_ << " h: "
        << returnReduce(found.size(), sumOp<label>()) << " wall cells";
    reduce(nMissed, sumOp<label>());
    if (nMissed > 0)
    {
        Info<< ", " << nMissed << " sample themselves (sample in the solid "
            << "or on another processor)";
    }
    Info<< endl;
}


void Foam::cutCellGeometry::calcWallQuantiles()
{
    wallIndex_.setSize(mesh_.nCells());
    wallIndex_ = -1;
    wallQuantiles_.clear();
    if (!wallLawCellAverage_)
    {
        return;
    }
    forAll(bodies_, bodyi)
    {
        if (!bodies_[bodyi].hasLevelSet())
        {
            WarningInFunction
                << "wallLawDistance cellAverage needs level-set bodies: the wall "
                << "law is taken at the centroid distance" << endl;
            return;
        }
    }

    // Members of every merged group
    auto group = [&](const label c)
    {
        return (linkMaster_[c] >= 0) ? linkMaster_[c] : c;
    };
    Map<DynamicList<label>> members;
    forAll(linkMaster_, celli)
    {
        if (linkMaster_[celli] >= 0)
        {
            members(linkMaster_[celli]).append(celli);
        }
    }

    // Sub-samples of the fluid of every group with a wall cell, weighted by
    // the volume they stand for; 12 per direction is enough for 16
    // quantiles
    const label nSub = min(nSubSamples_, label(12));
    const scalarField& V = mesh_.V();
    label nWall = 0;
    for (const label celli : wallCells_)
    {
        wallIndex_[celli] = nWall++;
    }
    DynamicField<point> pts;
    DynamicList<scalar> wts;
    labelList start(nWall + 1, 0);
    DynamicList<point> cellPts;
    for (const label celli : wallCells_)
    {
        const label m = group(celli);
        DynamicList<label> cells(1, m);
        const auto iter = members.cfind(m);
        if (iter.found()) cells.append(iter.val());
        for (const label c : cells)
        {
            if (alpha_[c] <= 0) continue;
            cellSamplePoints(c, nSub, cellPts);
            for (const point& p : cellPts)
            {
                pts.append(p);
                wts.append(V[c]/max(cellPts.size(), label(1)));
            }
        }
        start[wallIndex_[celli] + 1] = pts.size();
    }
    scalarField d;
    wallDistance(pts, d);

    wallQuantiles_.setSize(nWall*nWallQuantiles_, 0);
    for (const label celli : wallCells_)
    {
        const label wi = wallIndex_[celli];
        DynamicList<scalar> dk;
        DynamicList<scalar> wk;
        for (label i = start[wi]; i < start[wi + 1]; ++i)
        {
            if (d[i] > 0)
            {
                dk.append(d[i]);
                wk.append(wts[i]);
            }
        }
        if (dk.empty())
        {
            for (label q = 0; q < nWallQuantiles_; ++q)
            {
                wallQuantiles_[wi*nWallQuantiles_ + q] = dWall_[celli];
            }
            continue;
        }
        labelList order(sortedOrder(dk));
        scalar total = 0;
        for (const scalar w : wk) total += w;

        // Distance at the midpoint of each of nWallQuantiles_ equal volumes
        scalar cum = 0;
        label q = 0;
        forAll(order, j)
        {
            const scalar w = wk[order[j]];
            while
            (
                q < nWallQuantiles_
             && cum + w >= (q + 0.5)/nWallQuantiles_*total
            )
            {
                wallQuantiles_[wi*nWallQuantiles_ + q] = dk[order[j]];
                ++q;
            }
            cum += w;
        }
        for (; q < nWallQuantiles_; ++q)
        {
            wallQuantiles_[wi*nWallQuantiles_ + q] = dk[order.last()];
        }
    }

    Info<< "    wall law averaged over the fluid of the wall cells ("
        << nWallQuantiles_ << " wall-distance quantiles)" << endl;
}


Foam::scalar Foam::cutCellGeometry::cellAverageUTau
(
    const volVectorField& U,
    const volScalarField& nu,
    const label celli
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

    const vector n(Sw_[celli]/max(Awall_[celli], VSMALL));
    vector slip(U[celli]);
    if (moving_)
    {
        slip -= UbPtr_()[celli];
    }
    const scalar magUt = mag(slip - (slip & n)*n);
    return saUTauAverage
    (
        magUt,
        nu[celli],
        &wallQuantiles_[wallIndex_[celli]*nWallQuantiles_],
        nWallQuantiles_
    );
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

    for (const label celli : wallCells_)
    {
        const label s = sampleCell_[celli];
        const scalar utCell =
            wallQuantiles_.size()
          ? cellAverageUTau(U, nu, celli)
          : wallFunctionUTau
            (
                U, nu, celli, celli,
                (wallLawDistance_ == "halfCell")
              ? 0.5*wallCellSize_[celli]
              : (wallLawDistance_ == "groupCentroid" && linkWallDist_[celli] > 0)
              ? linkWallDist_[celli]
              : dWall_[celli]
            );
        const scalar utSample =
            (s < 0)
          ? utCell
          : wallFunctionUTau(U, nu, celli, s, sampleDist_[celli]);
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
            (s < 0 || w <= 0)
          ? utCell
          : ::sqrt((1 - w)*sqr(utCell) + w*sqr(utSample));
    }
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
        // value Spalding's law gives nu_w >= nu (u+ <= y+); with a sample
        // further out the wall cell may be faster than the linear law, and a
        // floor nu would add to the wall shear of the wall law. A sample
        // also leaves nu_w unbounded as the wall cell's own slip vanishes:
        // the implicit wall drag then holds the cell at rest, a second
        // steady state. On the wall law nu_w = u_tau d/u+(d+) stays below
        // nu + kappa u_tau d (u+ >= ln(1 + kappa d+)/kappa), which caps it
        // with the cell size for d
        vector slip(U[celli]);
        if (moving_)
        {
            slip -= UbPtr_()[celli];
        }
        const vector n(Sw_[celli]/max(Awall_[celli], VSMALL));
        const scalar magUt = max(mag(slip - (slip & n)*n), VSMALL);
        const scalar nuLaw = sqr(uTau[celli])*dWall_[celli]/magUt;
        nuW[celli] =
            (sampleCell_[celli] < 0)
          ? max(nu[celli], nuLaw)
          : min(nuLaw, nu[celli] + 0.41*uTau[celli]*wallCellSize_[celli]);
    }

    // Diagnostics at write times: the kinematic wall shear the wall cells
    // apply, nu_w |u_t|/dWall, and the weight of the wall law
    if (mesh_.time().writeTime())
    {
        volScalarField tauW
        (
            IOobject("ibmWallShear", mesh_.time().timeName(), mesh_),
            mesh_,
            dimensionedScalar(sqr(dimVelocity), Zero),
            "zeroGradient"
        );
        volScalarField weight
        (
            IOobject("ibmWallLawWeight", mesh_.time().timeName(), mesh_),
            mesh_,
            dimensionedScalar(dimless, Zero),
            "zeroGradient"
        );
        for (const label celli : wallCells_)
        {
            vector slip(U[celli]);
            if (moving_)
            {
                slip -= UbPtr_()[celli];
            }
            const vector n(Sw_[celli]/max(Awall_[celli], VSMALL));
            tauW[celli] =
                nuW[celli]*mag(slip - (slip & n)*n)/max(dWall_[celli], VSMALL);
            weight[celli] = w[celli];
        }
        tauW.write();
        weight.write();
    }

    tnuW.ref().correctBoundaryConditions();
    return tnuW;
}
