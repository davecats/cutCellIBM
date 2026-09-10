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
    // zero elsewhere. For a wall cell it is evaluated at its wall segment,
    // x_w = x_c + dWall n_w: the wall flux u_b.Sw of a rotating body must
    // vanish for a tangential motion, which the cell-centre value does not
    // give (an O(h) spurious mass source alternating in sign around the
    // body). Solid cells use their centre.
    vectorField& Ub = UbPtr_().primitiveFieldRef();
    Ub = Zero;

    pointField wallPoints(mesh_.cellCentres() + centroidOffset_.primitiveField());
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

    // Wall flux S = u_b . Sw (Sw points from the fluid into the body): the
    // volume the wall sweeps per unit time, taken as PRIMARY. Continuity
    // becomes sum_f theta_f phi_f = -S, and the fluid volume at the start of
    // the step is integrated from it, alpha^n V = alpha^{n+1} V - dt S, so
    // that unsteady and convective terms cancel identically for a fluid
    // moving rigidly with the body (discrete geometric conservation law).
    scalarField& S = massSourcePtr_().primitiveFieldRef();
    scalarField alphaOld(alpha_.primitiveField());

    alphaOldMin_ = GREAT;
    nSwept_ = 0;
    scalar sumS = 0;

    forAll(S, celli)
    {
        const scalar Sint = Ub[celli] & Sw[celli];
        S[celli] = Sint/V[celli];
        alphaOld[celli] = alpha_[celli] - dt*Sint/V[celli];
        sumS += Sint;

        if (!solidCells_.test(celli))
        {
            alphaOldMin_ = min(alphaOldMin_, alphaOld[celli]);
            if (alphaOld[celli] < 0) ++nSwept_;
        }
    }
    reduce(alphaOldMin_, minOp<scalar>());
    reduce(nSwept_, sumOp<label>());
    reduce(sumS, sumOp<scalar>());
    massSourcePtr_().correctBoundaryConditions();

    // The swept volume is what fvm::ddt(alpha, U) must see as the old
    // fluid fraction
    volScalarField& alpha0 = alpha_.oldTime();
    alpha0.primitiveFieldRef() = alphaOld;
    alpha0.correctBoundaryConditions();

    Info<< "    mass source: sum over the domain " << sumS
        << " (zero to roundoff by closure)" << endl;
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
