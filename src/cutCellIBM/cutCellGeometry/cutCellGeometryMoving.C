/*---------------------------------------------------------------------------*\
    cutCellIBM: moving-body part of cutCellGeometry (doc/moving_ibm.pdf)

    Nothing here is executed for static bodies: update() returns at once
    when no body moves, and the moving-only fields are never allocated.
\*---------------------------------------------------------------------------*/

#include "cutCellGeometry.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cutCellGeometry::update()
{
    if (!moving_)
    {
        return;
    }

    const scalar t = mesh_.time().value();
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

    // Geometry at the new time level. The old geometry is never needed: the
    // fluid volume the cells had at the start of the step is recovered from
    // the wall flux (calcMovingTerms), which is what makes the discrete
    // geometric conservation law exact.
    calcGeometry();
    calcMovingTerms();

    Info<< "    cells: "
        << returnReduce(solidCells_.count(), sumOp<label>()) << " solid, "
        << returnReduce(cutCells_.count(), sumOp<label>()) << " cut; "
        << "smallest old fluid fraction " << alphaOldMin_
        << ", cells swept through by more than they hold: " << nSwept_
        << endl;
}


void Foam::cutCellGeometry::calcMovingTerms()
{
    const scalar dt = mesh_.time().deltaTValue();
    const scalarField& V = mesh_.V();
    const vectorField& Sw = Sw_.primitiveField();

    // Body velocity for the cells a body owns (cut, solid and wall cells),
    // zero elsewhere, for the no-slip wall traction and the blanking. For a
    // wall cell it is evaluated at the foot of the wall normal through the
    // cell centre, x_w = x_P + dWall n_w, the point the two-point wall
    // gradient differences against. Solid cells use their centre. The mass
    // source does not use it for rotating bodies (see below).
    vectorField& Ub = UbPtr_().primitiveFieldRef();
    Ub = Zero;

    pointField wallPoints(mesh_.cellCentres());
    for (const label celli : wallCells_)
    {
        wallPoints[celli] += dWall_[celli]*Sw[celli]/max(Awall_[celli], VSMALL);
    }

    List<vectorField> bodyU(bodies_.size());
    forAll(bodies_, bodyi)
    {
        if (bodies_[bodyi].moving())
        {
            bodies_[bodyi].velocity(wallPoints, bodyU[bodyi]);
        }
    }
    forAll(cellBody_, celli)
    {
        const label bodyi = cellBody_[celli];
        if (bodyi >= 0 && bodies_[bodyi].moving())
        {
            Ub[celli] = bodyU[bodyi][celli];
        }
    }
    UbPtr_().correctBoundaryConditions();

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
    // rotation, and alpha-field bodies, keep S = u_b . Sw.
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
        if (anyRotating && bodyi >= 0 && bodies_[bodyi].rotating())
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


const Foam::volVectorField& Foam::cutCellGeometry::Ub() const
{
    if (!UbPtr_)
    {
        FatalErrorInFunction
            << "No moving body: the body velocity field does not exist"
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
