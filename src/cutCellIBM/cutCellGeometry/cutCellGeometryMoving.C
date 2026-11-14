/*---------------------------------------------------------------------------*\
    cutCellIBM: moving-body part of cutCellGeometry (doc/moving_ibm.pdf)

    Nothing here is executed for static bodies: update() returns at once
    when no body moves, and the moving-only fields are never allocated.
    Invariant bodies (static geometry with a wall velocity) only use
    calcWallVelocity and the invariance check.
\*---------------------------------------------------------------------------*/

#include "cutCellGeometry.H"
#include "fvMesh.H"
#include "Time.H"
#include "boundBox.H"
#include "OStringStream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cutCellGeometry::update()
{
    if (!wallVelocity_)
    {
        return;
    }

    const scalar t = mesh_.time().value();

    if (!moving_)
    {
        // Invariant bodies only: the geometry stays, the wall velocity
        // follows the motion
        forAll(bodies_, bodyi)
        {
            bodies_[bodyi].moveTo(t);
        }
        checkInvariance(0);
        calcWallVelocity();
        return;
    }

    Info<< "Cut-cell immersed boundary: moving the bodies to t = " << t << endl;

    forAll(bodies_, bodyi)
    {
        bodies_[bodyi].moveTo(t);
        if (bodies_[bodyi].moving())
        {
            Info<< "    body " << bodies_[bodyi].name()
                << ": centre " << bodies_[bodyi].currentCentre()
                << ", velocity " << bodies_[bodyi].velocityOfCentre() << endl;
        }
    }

    // Keep the face fractions of the geometry the carried flux belongs to
    if (!thetaOldPtr_)
    {
        thetaOldPtr_.reset(new surfaceScalarField("ibmThetaOld", theta_));
    }
    else
    {
        thetaOldPtr_() = theta_;
    }

    solidOld_ = solidCells_;

    // Geometry at the new time level. The old geometry is never needed: the
    // fluid volume the cells had at the start of the step is recovered from
    // the wall flux (calcMovingTerms), which is what makes the discrete
    // geometric conservation law exact.
    calcGeometry();
    checkInvariance(0);
    calcMovingTerms();

    Info<< "    cells: "
        << returnReduce(solidCells_.count(), sumOp<label>()) << " solid, "
        << returnReduce(cutCells_.count(), sumOp<label>()) << " cut; "
        << "smallest old fluid fraction " << alphaOldMin_
        << ", cells swept through by more than they hold: " << nSwept_
        << endl;
}


void Foam::cutCellGeometry::updateMesh()
{
    forAll(bodies_, bodyi)
    {
        if (!bodies_[bodyi].hasLevelSet())
        {
            FatalErrorInFunction
                << "Body " << bodies_[bodyi].name() << " has no level set: "
                << "the fluid fraction of an alpha-field body cannot follow "
                << "a change of the mesh topology. Use a level-set body "
                << "with adaptive refinement." << exit(FatalError);
        }
    }

    Info<< "Cut-cell immersed boundary: rebuilding the geometry on the "
        << "changed mesh" << endl;

    // The registered fields have been mapped to the new mesh (and are about
    // to be overwritten); the plain lists are resized here
    const label nCells = mesh_.nCells();
    cellBody_.setSize(nCells);
    nearestBody_.setSize(nCells);
    solidCells_.resize(nCells);
    cutCells_.resize(nCells);
    wallCells_.resize(nCells);

    // Mapped face values are meaningless on the faces the change created
    dCentroid_ = mesh_.delta();
    cutDeltaCoeffs_ = mesh_.nonOrthDeltaCoeffs();
    faceValueOffset_ = dimensionedVector(dimLength, Zero);

    calcGeometry();
    report();

    // Static bodies: the old time levels of alpha are the current one (a
    // refined cut cell inherited the fraction of its parent). Moving
    // bodies get theirs from the next update(). Invariant bodies need the
    // wall velocity on the changed cells now.
    if (!moving_)
    {
        if (wallVelocity_)
        {
            calcWallVelocity();
        }

        volScalarField* a0 = &alpha_;
        while (a0->nOldTimes())
        {
            a0 = &a0->oldTime();
            a0->primitiveFieldRef() = alpha_.primitiveField();
            a0->correctBoundaryConditions();
        }
    }
}


void Foam::cutCellGeometry::calcWallVelocity()
{
    const vectorField& Sw = Sw_.primitiveField();

    // Body velocity for the cells a body owns (cut, solid and wall cells),
    // zero elsewhere, for the no-slip wall traction and the blanking. For a
    // wall cell it is evaluated at the foot of the wall normal through the
    // fluid centroid, x_w = x_c + dWall n_w, the point the two-point wall
    // gradient differences against. Solid cells use their centre. The mass
    // source does not use it for rotating bodies (see below).
    //
    // Without the second-order treatment dWall = alpha V/(2 Awall) is an
    // effective distance (with it, the centroid distance): in a nearly full
    // cell with a tiny wall segment it exceeds the cell size by orders of
    // magnitude and x_w would lie far outside the cell (in a periodic domain
    // possibly beyond the minimum image). The wall of such a cell lies inside
    // the cell, so the distance is limited to where the normal leaves the
    // cell's bounding box; the traction coefficient keeps dWall.
    vectorField& Ub = UbPtr_().primitiveFieldRef();
    Ub = Zero;

    // (second order: from the value point, the group centroid of a linked
    // cell, which dWall is measured from)
    pointField wallPoints
    (
        mesh_.cellCentres()
      + (secondOrder_ ? valueOffset_ : centroidOffset_).primitiveField()
    );
    label nClipped = 0;
    for (const label celli : wallCells_)
    {
        const scalar dExit = cellExitDistance
        (
            celli, wallPoints[celli], Sw[celli]/max(Awall_[celli], VSMALL)
        );
        if (dWall_[celli] <= dExit)
        {
            wallPoints[celli] +=
                dWall_[celli]*Sw[celli]/max(Awall_[celli], VSMALL);
        }
        else
        {
            wallPoints[celli] += dExit*Sw[celli]/max(Awall_[celli], VSMALL);
            ++nClipped;
        }
    }

    List<vectorField> bodyU(bodies_.size());
    forAll(bodies_, bodyi)
    {
        if (bodies_[bodyi].hasWallVelocity())
        {
            bodies_[bodyi].velocity(wallPoints, bodyU[bodyi]);
        }
    }
    forAll(cellBody_, celli)
    {
        const label bodyi = cellBody_[celli];
        if (bodyi >= 0 && bodies_[bodyi].hasWallVelocity())
        {
            Ub[celli] = bodyU[bodyi][celli];
        }
    }
    UbPtr_().correctBoundaryConditions();
    reduce(nClipped, sumOp<label>());
    if (nClipped > 0)
    {
        Info<< "    wall velocity: evaluation point limited to the cell in "
            << nClipped << " wall cells (tiny wall segments)" << endl;
    }
}


void Foam::cutCellGeometry::calcMovingTerms()
{
    const scalar dt = mesh_.time().deltaTValue();
    const scalarField& V = mesh_.V();
    const vectorField& Sw = Sw_.primitiveField();

    calcWallVelocity();
    const vectorField& Ub = UbPtr_().primitiveField();

    // Wall flux S (Sw points from the fluid into the body): the volume the
    // wall sweeps per unit time, taken as PRIMARY. Continuity becomes
    // sum_f theta_f phi_f = -S, and the fluid volume at the start of the
    // step is integrated from it, alpha^n V = alpha^{n+1} V - dt S, so that
    // unsteady and convective terms cancel identically for a fluid moving
    // rigidly with the body (discrete geometric conservation law).
    //
    // Translation: S = U . Sw, exact and uniform (it includes the faces
    // closed by sliver absorption, which carry the swept volume of the
    // absorbed cell). Rotation: no single evaluation point of u_b gives the
    // flux through the wall (u_b(x) . Sw depends on the offset of x across
    // the wall normal), so it is computed as the flux of omega x (x - c)
    // through the wall faces, S_rot (rotationalWallFlux). It vanishes for a
    // body of revolution turning about its axis away from removed cells,
    // and sums to zero over the domain for any rigid motion. Bodies without
    // rotation, and alpha-field bodies, keep S = u_b . Sw. Invariant bodies
    // have no wall flux (S = 0): their wall velocity is tangential to the
    // exact wall, and u_b . Sw of the discrete wall would only be its error.
    scalarField Srot;
    bool anyRotating = false;
    forAll(bodies_, bodyi)
    {
        if (bodies_[bodyi].rotating()) anyRotating = true;
    }
    if (anyRotating)
    {
        rotationalWallFlux(Srot);
    }

    scalarField& S = massSourcePtr_().primitiveFieldRef();
    scalarField alphaOld(alpha_.primitiveField());

    alphaOldMin_ = GREAT;
    nSwept_ = 0;
    scalar sumS = 0;
    scalar maxS = 0;
    scalar maxSrot = 0;
    scalar maxDalpha = 0;

    forAll(S, celli)
    {
        const label bodyi = cellBody_[celli];
        scalar Sint;
        if (bodyi >= 0 && bodies_[bodyi].invariant())
        {
            Sint = 0;
        }
        else if (anyRotating && bodyi >= 0 && bodies_[bodyi].rotating())
        {
            Sint =
                (bodies_[bodyi].velocityOfCentre() & Sw[celli])
              + Srot[celli];
        }
        else if (anyRotating)
        {
            // S_rot is zero here unless the cell has no body of its own but
            // touches the cut of a rotating one
            Sint = (Ub[celli] & Sw[celli]) + Srot[celli];
        }
        else
        {
            Sint = Ub[celli] & Sw[celli];
        }
        if (anyRotating) maxSrot = max(maxSrot, mag(Srot[celli]));
        S[celli] = Sint/V[celli];
        alphaOld[celli] = alpha_[celli] - dt*Sint/V[celli];
        sumS += Sint;
        maxS = max(maxS, mag(Sint));
        maxDalpha = max(maxDalpha, mag(alpha_[celli] - alphaOld[celli]));

        if (!solidCells_.test(celli))
        {
            alphaOldMin_ = min(alphaOldMin_, alphaOld[celli]);
            if (alphaOld[celli] < 0) ++nSwept_;
        }
    }
    reduce(alphaOldMin_, minOp<scalar>());
    reduce(nSwept_, sumOp<label>());
    reduce(sumS, sumOp<scalar>());
    reduce(maxS, maxOp<scalar>());
    reduce(maxDalpha, maxOp<scalar>());
    massSourcePtr_().correctBoundaryConditions();

    if (anyRotating)
    {
        reduce(maxSrot, maxOp<scalar>());
        Info<< "    rotational wall flux: max |S_rot| " << maxSrot << endl;
    }

    // The swept volume is what fvm::ddt(alpha, U) must see as the old
    // fluid fraction
    volScalarField& alpha0 = alpha_.oldTime();
    alpha0.primitiveFieldRef() = alphaOld;
    alpha0.correctBoundaryConditions();

    Info<< "    mass source: sum over the domain " << sumS
        << " (zero to roundoff by closure), max |S| " << maxS
        << ", max |alpha^{n+1} - alpha^n| " << maxDalpha << endl;
}


Foam::scalar Foam::cutCellGeometry::cellExitDistance
(
    const label celli,
    const point& x0,
    const vector& n
) const
{
    // Distance from x0 (inside the cell) along n to the boundary of the
    // cell's bounding box, over the solved directions
    const boundBox bb
    (
        mesh_.cells()[celli].points(mesh_.faces(), mesh_.points()),
        false
    );
    const Vector<label>& solD = mesh_.solutionD();
    scalar t = GREAT;
    for (direction d = 0; d < 3; ++d)
    {
        if (solD[d] < 0 || mag(n[d]) < SMALL) continue;
        const scalar face = (n[d] > 0) ? bb.max()[d] : bb.min()[d];
        t = min(t, max((face - x0[d])/n[d], scalar(0)));
    }
    return t;
}


void Foam::cutCellGeometry::rotationalWallFlux(scalarField& Srot) const
{
    // S_rot,P = -sum_f theta_f (omega x (x'_f - c)) . Sf over the OPEN faces
    // of P (Sf out of P, x'_f the wet-face centroid of the geometric cut):
    // the flux of the (divergence-free, linear) rotational velocity through
    // the wall of the post-processed fluid cell, i.e. through the chord of
    // the geometric cut and through the geometric wet parts of the faces
    // closed by sliver absorption, which move with the body. This is the
    // rotational counterpart of U . Sw, which also includes the closed
    // faces. Every open face shared by two cells of the same body enters
    // both with opposite sign, and the two copies of a periodic face are
    // mapped to the same minimum image about the centre, hence
    // sum_P S_rot,P = 0 for any motion and mesh (except on faces between
    // cells of two different rotating bodies). Cells without fluid have all
    // faces closed and no flux.
    //
    // The geometric form over all wet faces with the flux of removed cells
    // handed to their neighbours (doc/openIssues) is exact for a body
    // turning in place, but a cell next to a removed one then keeps the flux
    // through its closed face, O(|omega| R h), which the discrete continuity
    // cannot carry: the Taylor-Couette pressure errors grew threefold. Here
    // a body turning in place has S = 0 except in the cells next to removed
    // ones (doc/secondOrderStudy.md).
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const vectorField& Sf = mesh_.faceAreas();
    const cellList& cells = mesh_.cells();

    Srot.setSize(mesh_.nCells());
    Srot = 0;

    // Body whose rotation a cell sees: its own body, or for a cell without
    // one that still touches the geometric cut, the nearest body. Cells away
    // from the bodies are skipped (their flux is zero to round-off).
    auto rotatingBody = [&](const label celli) -> label
    {
        label bodyi = cellBody_[celli];
        if (bodyi < 0)
        {
            for (const label facei : cells[celli])
            {
                if (thetaGeo_[facei] < 1)
                {
                    bodyi = nearestBody_[celli];
                    break;
                }
            }
        }
        return (bodyi >= 0 && bodies_[bodyi].rotating()) ? bodyi : -1;
    };

    // Flux of the rotation of the body of cell celli through the wet part
    // of face facei, with Sf
    auto faceFlux = [&](const label celli, const label facei, const scalar t)
    {
        const label bodyi = rotatingBody(celli);
        if (bodyi < 0) return scalar(0);
        const ibmBody& b = bodies_[bodyi];
        const vector u(b.omega() ^ b.relativePosition(wetCentreGeo_[facei]));
        return t*(u & Sf[facei]);
    };

    forAll(nei, facei)
    {
        const scalar t = theta_[facei];
        if (t <= 0) continue;
        Srot[own[facei]] -= faceFlux(own[facei], facei, t);
        Srot[nei[facei]] += faceFlux(nei[facei], facei, t);
    }
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];
        const fvsPatchScalarField& tp = theta_.boundaryField()[patchi];
        forAll(fvp, i)
        {
            if (tp[i] <= 0) continue;
            const label facei = fvp.start() + i;
            Srot[own[facei]] -= faceFlux(own[facei], facei, tp[i]);
        }
    }
}


namespace
{
    // Edge crossings of the cut that belong to body bodyi: those at which
    // its level set is the smallest in magnitude among the level-set
    // bodies, with the edge lengths
    void bodyCrossings
    (
        const Foam::fvMesh& mesh,
        const Foam::PtrList<Foam::ibmBody>& bodies,
        const Foam::pointField& edgeCross,
        const Foam::label bodyi,
        Foam::pointField& pts,
        Foam::scalarField& h
    )
    {
        using namespace Foam;
        DynamicList<point> x;
        DynamicList<scalar> l;
        const edgeList& edges = mesh.edges();
        forAll(edgeCross, edgei)
        {
            if (edgeCross[edgei] == point::max) continue;
            x.append(edgeCross[edgei]);
            l.append(edges[edgei].mag(mesh.points()));
        }
        pointField xc(std::move(x));
        scalarField lc(std::move(l));

        scalarField phiMin(xc.size(), GREAT);
        scalarField phiOwn(xc.size(), GREAT);
        scalarField phib;
        forAll(bodies, bj)
        {
            if (!bodies[bj].hasLevelSet()) continue;
            bodies[bj].evaluate(xc, phib);
            phiMin = min(phiMin, mag(phib));
            if (bj == bodyi) phiOwn = mag(phib);
        }

        DynamicList<point> xo;
        DynamicList<scalar> lo;
        forAll(xc, i)
        {
            if (phiOwn[i] <= phiMin[i])
            {
                xo.append(xc[i]);
                lo.append(lc[i]);
            }
        }
        pts = std::move(xo);
        h = std::move(lo);
    }
}


Foam::scalar Foam::cutCellGeometry::invarianceTime() const
{
    // tau = (smallest edge length at the cut) / (largest wall speed), over
    // the invariant bodies: the virtual motion moves the fastest wall point
    // by one cell, so a body that is not invariant changes its level set
    // by a fraction of the edge length of the order of the normal part of
    // its velocity relative to the largest one
    scalar hMin = GREAT;
    scalar uMax = 0;
    forAll(bodies_, bodyi)
    {
        if (!bodies_[bodyi].invariant()) continue;
        pointField pts;
        scalarField h;
        bodyCrossings(mesh_, bodies_, edgeCrossGeo_, bodyi, pts, h);
        vectorField u;
        bodies_[bodyi].velocity(pts, u);
        if (h.size())
        {
            hMin = min(hMin, min(h));
            uMax = max(uMax, max(mag(u)));
        }
    }
    reduce(hMin, minOp<scalar>());
    reduce(uMax, maxOp<scalar>());
    return (uMax > 0 && hMin < GREAT) ? hMin/uMax : 0;
}


void Foam::cutCellGeometry::checkInvariance(const scalar tau) const
{
    // The level set of an invariant body at the edge crossings of its cut,
    // with the body at rest (the geometry) and in its motion state advanced
    // by tau, must agree: |d phi| <= invarianceTolerance h (h the edge
    // length; the level sets are distances, so d phi is the normal
    // displacement of the wall)
    forAll(bodies_, bodyi)
    {
        const ibmBody& b = bodies_[bodyi];
        if (!b.invariant()) continue;

        pointField pts;
        scalarField h;
        bodyCrossings(mesh_, bodies_, edgeCrossGeo_, bodyi, pts, h);

        scalar err = 0;
        label nPts = pts.size();
        if (pts.size())
        {
            scalarField dphi;
            b.levelSetChange(pts, tau, dphi);
            err = max(mag(dphi)/max(h, VSMALL));
        }
        reduce(err, maxOp<scalar>());
        reduce(nPts, sumOp<label>());

        if (tau > 0)
        {
            Info<< "    body " << b.name() << ": invariance under a virtual "
                << "motion over " << tau << " s (wall moved by one cell): "
                << "max |d phi|/h " << err << " at " << nPts
                << " edge crossings" << endl;
        }
        if (err > invarianceTolerance_)
        {
            OStringStream when;
            if (tau > 0)
            {
                when<< "under a virtual motion over " << tau << " s";
            }
            else
            {
                when<< "at t = " << mesh_.time().value();
            }
            FatalErrorInFunction
                << "Body " << b.name() << " is declared invariant, but its "
                << "motion changes its shape: the level set at the cut moves "
                << "by " << err << " edge lengths (invarianceTolerance "
                << invarianceTolerance_ << ") "
                << when.str().c_str()
                << ". Only a motion that maps the body onto itself (rotation "
                << "about an axis of symmetry, sliding along a direction of "
                << "translational symmetry) can be invariant."
                << exit(FatalError);
        }
    }
}


const Foam::volVectorField& Foam::cutCellGeometry::Ub() const
{
    if (!UbPtr_)
    {
        FatalErrorInFunction
            << "No moving or invariant body: the body velocity field does "
            << "not exist"
            << exit(FatalError);
    }
    return UbPtr_();
}


const Foam::volScalarField& Foam::cutCellGeometry::massSource() const
{
    if (!massSourcePtr_)
    {
        FatalErrorInFunction
            << "No moving body: the mass source field does not exist"
            << exit(FatalError);
    }
    return massSourcePtr_();
}


Foam::tmp<Foam::volVectorField> Foam::cutCellGeometry::movingWallSource
(
    const volScalarField& nu
) const
{
    // The u_b half of the no-slip wall traction nu (u_b - u_P) Awall/d, and
    // the blanking source that drives the solid interior to u_b. Both
    // vanish for a body at rest.
    tmp<volVectorField> tsrc(nu*(noSlipCoeff_ + blankCoeff_)*Ub());
    tsrc.ref().rename("ibmMovingWallSource");
    return tsrc;
}


Foam::tmp<Foam::volVectorField> Foam::cutCellGeometry::movingWallSource
(
    const volScalarField& nuWall,
    const volScalarField& nu
) const
{
    tmp<volVectorField> tsrc
    (
        (nuWall*noSlipCoeff_ + nu*blankCoeff_)*Ub()
    );
    tsrc.ref().rename("ibmMovingWallSource");
    return tsrc;
}


Foam::label Foam::cutCellGeometry::pressureRefCell
(
    const label userRefCell
) const
{
    bool userOk = (userRefCell >= 0) && !solidCells_.test(userRefCell);
    reduce(userOk, orOp<bool>());
    if (userOk)
    {
        return userRefCell;
    }

    // The fluid cell farthest from the bodies (largest level set), so that
    // a moving body does not swallow the reference; without a level set,
    // the first fluid cell
    scalar best = -GREAT;
    label bestCell = -1;
    const bool haveLevelSet = (cellLevelSet_.size() == mesh_.nCells());
    forAll(alpha_, celli)
    {
        if (solidCells_.test(celli)) continue;
        const scalar score = haveLevelSet ? cellLevelSet_[celli] : -scalar(celli);
        if (score > best)
        {
            best = score;
            bestCell = celli;
        }
    }

    const scalar gBest = returnReduce(best, maxOp<scalar>());
    label owner = (bestCell >= 0 && best == gBest) ? Pstream::myProcNo() : labelMax;
    owner = returnReduce(owner, minOp<label>());

    return (owner == Pstream::myProcNo()) ? bestCell : -1;
}
