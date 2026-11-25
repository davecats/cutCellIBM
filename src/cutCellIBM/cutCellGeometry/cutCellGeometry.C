#include "cutCellGeometry.H"
#include "ibmBodies.H"
#include "fvMesh.H"
#include "Time.H"
#include "fvc.H"
#include "fvmSup.H"
#include "syncTools.H"
#include "volPointInterpolation.H"
#include "pointFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "coupledFvPatchField.H"
#include "OStringStream.H"
#include "boundBox.H"
#include "globalIndex.H"
#include "coupledPolyPatch.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "clockTime.H"

namespace Foam
{
    defineTypeNameAndDebug(cutCellGeometry, 0);
}

namespace
{

// Lexicographic order of points, so that both sides of a coupled boundary
// bisect a shared edge in the same direction and obtain the same crossing.
inline bool lessPoint(const Foam::point& a, const Foam::point& b)
{
    if (a.x() != b.x()) return a.x() < b.x();
    if (a.y() != b.y()) return a.y() < b.y();
    return a.z() < b.z();
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cutCellGeometry::cutCellGeometry(const fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "ibmProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    bodies_(),
    // Linked small cells need no absorption beyond the cells without
    // fluid; without linking slivers are absorbed below 2 %
    alphaMin_
    (
        getOrDefault<scalar>
        (
            "alphaMin",
            getOrDefault<scalar>("mergeAlpha", 0.5) > 0 ? 0 : 0.02
        )
    ),
    nSubSamples_(getOrDefault<label>("nSubSamples", 32)),
    nBisect_(getOrDefault<label>("nBisect", 40)),
    writeFields_(getOrDefault<bool>("writeFields", true)),
    writeForces_(getOrDefault<bool>("writeForces", true)),
    alpha_
    (
        IOobject("ibmAlpha", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar(dimless, 1.0), "zeroGradient"
    ),
    theta_
    (
        IOobject("ibmTheta", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar(dimless, 1.0)
    ),
    Sw_
    (
        IOobject("ibmSw", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedVector(dimArea, Zero), "zeroGradient"
    ),
    Awall_
    (
        IOobject("ibmAwall", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar(dimArea, Zero), "zeroGradient"
    ),
    dWall_
    (
        IOobject("ibmWallDistance", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar(dimLength, GREAT), "zeroGradient"
    ),
    wallCoeff_
    (
        IOobject("ibmWallCoeff", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar(dimless/dimArea, Zero), "zeroGradient"
    ),
    noSlipCoeff_
    (
        IOobject("ibmNoSlipCoeff", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar(dimless/dimArea, Zero), "zeroGradient"
    ),
    blankCoeff_
    (
        IOobject("ibmBlankCoeff", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar(dimless/dimArea, Zero), "zeroGradient"
    ),
    bodyIndex_
    (
        IOobject("ibmBodyIndex", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedScalar(dimless, -1), "zeroGradient"
    ),
    secondOrder_(getOrDefault<bool>("secondOrder", true)),
    ibmConvection_(getOrDefault<bool>("ibmConvection", secondOrder_)),
    ibmConvectionWallLaw_(getOrDefault<bool>("ibmConvectionWallLaw", false)),
    zeroShear_(getOrDefault<bool>("zeroShear", true)),
    zeroShearTangential_(getOrDefault<bool>("zeroShearTangential", true)),
    centroidOffset_
    (
        IOobject("ibmCentroidOffset", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedVector(dimLength, Zero), "zeroGradient"
    ),
    valueOffset_
    (
        IOobject("ibmValueOffset", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedVector(dimLength, Zero), "zeroGradient"
    ),
    wetFaceOffset_
    (
        IOobject("ibmWetFaceOffset", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedVector(dimLength, Zero)
    ),
    dCentroid_
    (
        IOobject("ibmCentroidDelta", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh.delta()
    ),
    cutDeltaCoeffs_
    (
        IOobject("ibmDeltaCoeffs", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh.nonOrthDeltaCoeffs()
    ),
    centroidLevelSet_(),
    faceValueOffset_
    (
        IOobject("ibmFaceValueOffset", mesh.time().timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedVector(dimLength, Zero)
    ),
    pointPhiGeo_(),
    edgeCrossGeo_(),
    wallMoment_(),
    wallPieceSum_(),
    wallTorqueMoment_(),
    wetOffsetAll_(),
    cellBody_(mesh.nCells(), -1),
    nearestBody_(mesh.nCells(), -1),
    solidCells_(mesh.nCells()),
    cutCells_(mesh.nCells()),
    wallCells_(mesh.nCells()),
    fluidVolume_(0),
    upToDate_(false),
    forcesFilePtr_(nullptr),
    momentsFilePtr_(nullptr),
    moving_(false),
    wallVelocity_(false),
    invarianceTolerance_(getOrDefault<scalar>("invarianceTolerance", 1e-6)),
    periodicLengths_(getOrDefault<vector>("periodicLengths", Zero)),
    cellLevelSet_(),
    cacheStatic_(getOrDefault<bool>("cacheStaticGeometry", true)),
    hierSampling_(getOrDefault<bool>("hierarchicalSampling", true)),
    staticCacheValid_(false),
    staticPointPhi_(),
    staticFacePhi_(),
    staticCellPhi_(),
    staticEdgeCross_(),
    staticAlpha_(),
    staticCentroid_(),
    staticBodyPhi_(),
    cellRadius_(),
    UbPtr_(nullptr),
    massSourcePtr_(nullptr),
    thetaOldPtr_(nullptr),
    solidOld_(),
    alphaOldMin_(1),
    nSwept_(0),
    hPlusResolved_(getOrDefault<scalar>("hPlusResolved", 1)),
    hPlusModelled_(getOrDefault<scalar>("hPlusModelled", 10)),
    wallSampleDistance_(getOrDefault<scalar>("wallSampleDistance", 1.5)),
    wallCellSize_(),
    sampleCell_(),
    sampleDist_(),
    sampleSlot_(),
    nSampleSlots_(0),
    sampleSendCells_(),
    sampleRecvSlots_(),
    wallSamplesDone_(false),
    lsqDone_(false),
    mergeAlpha_(getOrDefault<scalar>("mergeAlpha", 0.5)),
    mergePenalty_(getOrDefault<scalar>("mergePenalty", 1e4)),
    linkFace_(),
    linkGroup_(),
    nbrLinkGroup_(),
    boundaryLink_(),
    linkRel_(),
    linkWallDist_(),
    groupSlot_(),
    slotGroup_(),
    sendSlots_(),
    recvSlots_(),
    splitGroupCells_(),
    nSplitGroups_(0)
{
    Info<< "Cut-cell immersed boundary: reading " << name() << nl
        << "    mergeAlpha " << mergeAlpha_
        << ", alphaMin " << alphaMin_
        << ", nSubSamples " << nSubSamples_
        << (hierSampling_ ? " (hierarchical)" : "")
        << ", nBisect " << nBisect_
        << (secondOrder_ ? ", second-order wall treatment" : "")
        << (ibmConvection_ ? ", convective face values from the value points" : "")
        << (ibmConvection_ && ibmConvectionWallLaw_ ? " (wall law in wall-function cells)" : "")
        << (zeroShear_ ? ", zero-shear free slip" : "")
        << endl;
    if (secondOrder_)
    {
        centroidOffset_.writeOpt(IOobject::AUTO_WRITE);
        if (mergeAlpha_ > 0) valueOffset_.writeOpt(IOobject::AUTO_WRITE);
    }
    else
    {
        WarningInFunction
            << "secondOrder false: the first-order wall treatment takes the "
            << "flux of a cut face at the face centre, so the cut-cell "
            << "pressure has an O(1) error that does not converge "
            << "(doc/secondOrderStudy.md). Velocity and forces converge at "
            << "first order." << endl;
    }
    checkSchemes();

    // Bodies
    const dictionary& bdict = subDict("bodies");
    bodies_.setSize(bdict.size());
    label nBodies = 0;
    label nLevelSet = 0;
    for (const entry& e : bdict)
    {
        if (!e.isDict())
        {
            FatalIOErrorInFunction(bdict)
                << "Entry " << e.keyword() << " is not a dictionary"
                << exit(FatalIOError);
        }
        bodies_.set(nBodies, ibmBody::New(e.keyword(), mesh_, e.dict()));
        if (bodies_[nBodies].hasLevelSet()) ++nLevelSet;
        ++nBodies;
    }
    bodies_.setSize(nBodies);

    if (nBodies == 0)
    {
        WarningInFunction << "No immersed bodies defined" << endl;
    }

    // Periodicity of the domain, for bodies crossing cyclic boundaries
    forAll(bodies_, bodyi)
    {
        const ibmBody& b = bodies_[bodyi];
        bodies_[bodyi].setPeriodicLengths(periodicLengths_);
        if (b.motion().invariant() && !b.hasLevelSet())
        {
            FatalErrorInFunction
                << "Body " << b.name() << ": an invariant motion needs a "
                << "level set, to check the invariance; an alpha-field body "
                << "has none" << exit(FatalError);
        }
        if (b.moving()) moving_ = true;
        if (b.hasWallVelocity()) wallVelocity_ = true;
    }
    reduce(moving_, orOp<bool>());
    reduce(wallVelocity_, orOp<bool>());

    if (moving_)
    {
        Info<< "    moving bodies: geometry rebuilt every time step"
            << (cmptMax(periodicLengths_) > 0
                ? ", periodic lengths " + Foam::name(periodicLengths_) : "")
            << endl;
    }
    else if (wallVelocity_)
    {
        Info<< "    invariant bodies only: static geometry, wall velocity "
            << "updated every time step, invariance tolerance "
            << invarianceTolerance_ << endl;
    }

    // Free bodies: mass properties and the state of a restart, before the
    // bodies are put in place
    initFreeBodies();

    if (wallVelocity_)
    {
        forAll(bodies_, bodyi)
        {
            bodies_[bodyi].moveTo(mesh_.time().value());
        }
        UbPtr_.reset
        (
            new volVectorField
            (
                IOobject("ibmUb", mesh.time().timeName(), mesh,
                         IOobject::NO_READ, IOobject::AUTO_WRITE),
                mesh, dimensionedVector(dimVelocity, Zero), "zeroGradient"
            )
        );
    }

    if (moving_)
    {
        massSourcePtr_.reset
        (
            new volScalarField
            (
                IOobject("ibmMassSource", mesh.time().timeName(), mesh,
                         IOobject::NO_READ, IOobject::NO_WRITE),
                mesh, dimensionedScalar(dimless/dimTime, Zero), "zeroGradient"
            )
        );
        // geometry written with the solution
        alpha_.writeOpt(IOobject::AUTO_WRITE);
        theta_.writeOpt(IOobject::AUTO_WRITE);
        Sw_.writeOpt(IOobject::AUTO_WRITE);
        bodyIndex_.writeOpt(IOobject::AUTO_WRITE);
    }

    calcGeometry();
    report();

    if (wallVelocity_)
    {
        // Invariant bodies: under a virtual motion that moves the wall by
        // one cell, then in the actual state at the start time
        const scalar tau = invarianceTime();
        if (tau > 0)
        {
            checkInvariance(tau);
            checkInvariance(0);
        }
    }

    if (moving_)
    {
        calcMovingTerms();
    }
    else if (wallVelocity_)
    {
        calcWallVelocity();
    }

    if (writeFields_)
    {
        writeGeometryFields();
    }
}


void Foam::cutCellGeometry::calcGeometry()
{
    // Wall-clock times of the stages (reported for moving geometry)
    clockTime clock;
    FixedList<scalar, 8> times(Zero);

    // Fluid fractions: union of all bodies. Level-set bodies are combined
    // exactly through the minimum of their level sets; alpha-field bodies
    // and the level-set group are combined through the minimum of their
    // fractions, exact wherever they do not share a cell.
    scalarField alpha(mesh_.nCells(), 1.0);
    scalarField theta(mesh_.nFaces(), 1.0);
    pointField wetCentre(mesh_.faceCentres());
    cellBody_ = -1;
    nearestBody_ = -1;
    vectorField centroidOffset(mesh_.nCells(), Zero);
    wetOffsetAll_.setSize(mesh_.nFaces());
    wetOffsetAll_ = Zero;

    // Level set at the points and edge crossings of the union: a point is
    // fluid if it is fluid for every body; on an edge cut by several, the
    // crossing nearest to the fluid end bounds the fluid
    pointPhiGeo_.setSize(mesh_.nPoints());
    pointPhiGeo_ = GREAT;
    edgeCrossGeo_.setSize(mesh_.nEdges());
    edgeCrossGeo_ = point::max;
    auto mergeCut = [&](const scalarField& phi, const pointField& ec)
    {
        const edgeList& edges = mesh_.edges();
        const pointField& points = mesh_.points();
        forAll(edges, edgei)
        {
            if (ec[edgei] == point::max) continue;
            const edge& e = edges[edgei];
            if (edgeCrossGeo_[edgei] == point::max)
            {
                edgeCrossGeo_[edgei] = ec[edgei];
            }
            else
            {
                const point& xf =
                    (phi[e[0]] > 0) ? points[e[0]] : points[e[1]];
                if (magSqr(ec[edgei] - xf) < magSqr(edgeCrossGeo_[edgei] - xf))
                {
                    edgeCrossGeo_[edgei] = ec[edgei];
                }
            }
        }
        pointPhiGeo_ = min(pointPhiGeo_, phi);
    };

    auto merge = [&](const scalarField& a, const scalarField& t,
                     const pointField& xw, const labelList& body)
    {
        forAll(alpha, celli)
        {
            if (a[celli] < alpha[celli])
            {
                alpha[celli] = a[celli];
                cellBody_[celli] = body[celli];
            }
        }
        forAll(theta, facei)
        {
            if (t[facei] < theta[facei])
            {
                theta[facei] = t[facei];
                wetCentre[facei] = xw[facei];
            }
        }
    };

    label nLevelSet = 0;
    forAll(bodies_, bodyi)
    {
        if (bodies_[bodyi].hasLevelSet()) ++nLevelSet;
    }

    if (nLevelSet > 0)
    {
        scalarField a, t, phi;
        pointField xw, ec;
        labelList body;
        calcFromLevelSet(a, t, xw, body, centroidOffset, phi, ec);
        merge(a, t, xw, body);
        mergeCut(phi, ec);
        cellCentreLevelSet(cellLevelSet_);
    }
    times[0] = clock.timeIncrement();

    forAll(bodies_, bodyi)
    {
        if (!bodies_[bodyi].hasLevelSet())
        {
            scalarField a, t, phi;
            pointField xw, ec;
            calcFromAlphaField
            (
                refCast<const alphaFieldBody>(bodies_[bodyi]),
                a, t, xw, phi, ec
            );
            merge(a, t, xw, labelList(mesh_.nCells(), bodyi));
            mergeCut(phi, ec);

            // Second order: fluid centroids of the cut cells of this body
            // (the value points; the cell centre otherwise)
            if (secondOrder_)
            {
                label nNoCentroid = 0;
                forAll(alpha, celli)
                {
                    if
                    (
                        cellBody_[celli] != bodyi
                     || alpha[celli] <= 0 || alpha[celli] >= 1
                    )
                    {
                        continue;
                    }
                    point c;
                    if (sampleInterpolatedCentroid(celli, phi, c))
                    {
                        centroidOffset[celli] = c - mesh_.cellCentres()[celli];
                    }
                    else
                    {
                        ++nNoCentroid;
                    }
                }
                reduce(nNoCentroid, sumOp<label>());
                if (nNoCentroid)
                {
                    WarningInFunction
                        << "alpha-field body " << bodies_[bodyi].name() << ": "
                        << nNoCentroid << " cut cells that are not axis-aligned"
                        << " hexahedra keep the cell centre as value point"
                        << endl;
                }
            }
        }
    }

    // The geometric cut as computed, before any post-processing: the exact
    // wet faces of the fluid polygon of every cell (rotational wall flux)
    thetaGeo_ = theta;
    wetCentreGeo_ = wetCentre;

    // Absorb slivers (C1): a deactivated cell also has its faces closed, or
    // the discrete divergence would see a flux into a cell without equation.
    // Cells without fluid always (alphaMin 0: only those).
    const cellList& cells = mesh_.cells();
    label nAbsorbed = 0;
    forAll(alpha, celli)
    {
        if (alpha[celli] <= 0 || alpha[celli] < alphaMin_)
        {
            if (alpha[celli] > 0) ++nAbsorbed;
            alpha[celli] = 0;
            for (const label facei : cells[celli])
            {
                theta[facei] = 0;
            }
        }
    }
    reduce(nAbsorbed, sumOp<label>());
    if (nAbsorbed > 0 || !moving_)
    {
        Info<< "    absorbed " << nAbsorbed << " sliver cells with alpha < "
            << alphaMin_ << endl;
    }

    // Coupled faces (C2): the two sides of a cyclic or processor face are the
    // same face and must carry the same fraction. Take the minimum, which
    // also propagates the closing of absorbed cells across the boundary.
    // The geometric fractions (before the absorption, which legitimately
    // closes one side first) must agree on both sides.
    {
        scalarField thetaGeoSynced(thetaGeo_);
        syncTools::syncFaceList(mesh_, thetaGeoSynced, minEqOp<scalar>());
        syncTools::syncFaceList(mesh_, theta, minEqOp<scalar>());
        scalar maxDiff = 0;
        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& fvp = mesh_.boundary()[patchi];
            if (fvp.coupled())
            {
                forAll(fvp, i)
                {
                    const label facei = fvp.start() + i;
                    maxDiff = max
                    (
                        maxDiff,
                        mag(thetaGeoSynced[facei] - thetaGeo_[facei])
                    );
                }
            }
        }
        reduce(maxDiff, maxOp<scalar>());
        if (maxDiff > 1e-10)
        {
            WarningInFunction
                << "The fluid fraction differs across coupled faces by up to "
                << maxDiff << ". The body is not periodic/consistent across "
                << "the coupled boundary; the minimum has been taken." << endl;
        }
    }

    // Fluid cells without an open face to another cell (all their internal
    // and coupled faces closed, e.g. by absorbed neighbours) would have an
    // empty pressure row: blank them like slivers, whatever their alpha.
    // Their internal and coupled faces are closed already, so the
    // synchronisation above stays valid.
    {
        label nClosed = 0;
        forAll(alpha, celli)
        {
            if (alpha[celli] <= 0) continue;
            bool open = false;
            for (const label facei : cells[celli])
            {
                if (theta[facei] <= 0) continue;
                if (facei < mesh_.nInternalFaces())
                {
                    open = true;
                    break;
                }
                const label patchi =
                    mesh_.boundaryMesh().whichPatch(facei);
                if (mesh_.boundary()[patchi].coupled())
                {
                    open = true;
                    break;
                }
            }
            if (!open)
            {
                ++nClosed;
                alpha[celli] = 0;
                for (const label facei : cells[celli])
                {
                    theta[facei] = 0;
                }
            }
        }
        reduce(nClosed, sumOp<label>());
        if (nClosed > 0)
        {
            Info<< "    blanked " << nClosed
                << " fluid cells without an open face" << endl;
        }
    }

    alpha_.primitiveFieldRef() = alpha;
    alpha_.correctBoundaryConditions();

    // Fluid centroids (cell centre in fluid cells); solid cells keep zero
    forAll(alpha, celli)
    {
        if (alpha[celli] <= 0) centroidOffset[celli] = Zero;
    }
    centroidOffset_.primitiveFieldRef() = centroidOffset;
    centroidOffset_.correctBoundaryConditions();
    wetFaceOffset_.primitiveFieldRef() =
        SubField<vector>(wetOffsetAll_, mesh_.nInternalFaces());
    forAll(wetFaceOffset_.boundaryField(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];
        wetFaceOffset_.boundaryFieldRef()[patchi] =
            SubField<vector>(wetOffsetAll_, fvp.size(), fvp.start());
    }

    theta_.primitiveFieldRef() = SubField<scalar>(theta, mesh_.nInternalFaces());
    forAll(theta_.boundaryField(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];
        theta_.boundaryFieldRef()[patchi] =
            SubField<scalar>(theta, fvp.size(), fvp.start());
    }

    times[1] = clock.timeIncrement();
    calcDerived();
    times[2] = clock.timeIncrement();
    calcLinks();
    times[3] = clock.timeIncrement();
    if (secondOrder_)
    {
        calcValueOffset();
        times[4] = clock.timeIncrement();
        calcCentroidGeometry();
        times[5] = clock.timeIncrement();
        calcWallPieces();
        times[6] = clock.timeIncrement();
    }
    if (zeroShear_)
    {
        calcZeroShear();
        times[7] = clock.timeIncrement();
    }
    if (moving_)
    {
        for (scalar& x : times)
        {
            reduce(x, maxOp<scalar>());
        }
        Info<< "    geometry times [s, max over ranks]: level set "
            << times[0] << ", post " << times[1] << ", derived " << times[2]
            << ", links " << times[3] << ", value offset " << times[4]
            << ", centroids " << times[5] << ", wall pieces " << times[6]
            << ", zero shear " << times[7] << endl;
    }
    wallSamplesDone_ = false;
    lsqDone_ = false;

    upToDate_ = true;
}


Foam::label Foam::cutCellGeometry::outwardFace
(
    const label celli,
    const scalarField& nbrAlpha
) const
{
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const vectorField& Sf = mesh_.faceAreas();
    const scalarField& alpha = alpha_.primitiveField();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const label nInternal = mesh_.nInternalFaces();
    const vector nw(Sw_[celli]/max(Awall_[celli], VSMALL));

    label bestFace = -1;
    scalar best = -GREAT;
    for (const label facei : mesh_.cells()[celli])
    {
        // Wet fraction of the face and fluid fraction across it: internal
        // faces, and coupled faces whose transform is a translation
        // (processor, translational cyclic)
        scalar t, alphaOther;
        vector sOut(Sf[facei]);
        if (facei < nInternal)
        {
            t = theta_[facei];
            const bool isOwn = (own[facei] == celli);
            alphaOther = alpha[isOwn ? nei[facei] : own[facei]];
            if (!isOwn) sOut = -sOut;
        }
        else
        {
            const label patchi = pbm.whichPatch(facei);
            const coupledPolyPatch* cpp = isA<coupledPolyPatch>(pbm[patchi]);
            if (!cpp || !cpp->parallel()) continue;
            t = theta_.boundaryField()[patchi][facei - pbm[patchi].start()];
            alphaOther = nbrAlpha[facei - nInternal];
        }
        if (t <= 0) continue;
        if (alphaOther <= alpha[celli]) continue;
        const scalar wet = t*mag(Sf[facei]);
        const scalar score =
            wet*(1 + max(scalar(0), -(nw & sOut)/max(mag(sOut), VSMALL)))
           *alphaOther;
        // Equal scores (symmetric geometry) are decided by the face centre,
        // not by the face order, which differs between decompositions
        const scalar tol = 1e-10*max(score, best);
        bool better = bestFace < 0 || score > best + tol;
        if (!better && score >= best - tol)
        {
            const point& a = mesh_.faceCentres()[facei];
            const point& b = mesh_.faceCentres()[bestFace];
            for (direction d = 0; d < vector::nComponents; ++d)
            {
                const scalar tolX = 1e-10*(mag(a[d]) + mag(b[d]) + SMALL);
                if (mag(a[d] - b[d]) > tolX)
                {
                    better = a[d] < b[d];
                    break;
                }
            }
        }
        if (better)
        {
            best = score;
            bestFace = facei;
        }
    }
    return bestFace;
}


void Foam::cutCellGeometry::calcLinks()
{
    const label nCells = mesh_.nCells();
    const label nInternal = mesh_.nInternalFaces();
    const label nBoundary = mesh_.nBoundaryFaces();
    linkFace_.setSize(nCells);
    linkFace_ = -1;
    linkGroup_.setSize(nCells);
    linkGroup_ = -1;
    nbrLinkGroup_.setSize(nBoundary);
    nbrLinkGroup_ = -1;
    boundaryLink_.setSize(nBoundary);
    boundaryLink_ = false;
    linkRel_.setSize(nCells);
    linkRel_ = Zero;
    linkWallDist_.setSize(nCells);
    linkWallDist_ = 0;
    groupSlot_.setSize(nCells);
    groupSlot_ = -1;
    slotGroup_.clear();
    sendSlots_.setSize(Pstream::nProcs());
    recvSlots_.setSize(Pstream::nProcs());
    forAll(sendSlots_, proci)
    {
        sendSlots_[proci].clear();
        recvSlots_[proci].clear();
    }
    splitGroupCells_.clear();
    nSplitGroups_ = 0;

    if (mergeAlpha_ <= 0)
    {
        return;
    }

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const labelUList& faceOwner = mesh_.faceOwner();
    const scalarField& alpha = alpha_.primitiveField();
    const scalarField& V = mesh_.V();
    const vectorField& C = mesh_.C().primitiveField();
    const vectorField& off = centroidOffset_.primitiveField();
    const globalIndex cellIndex(nCells);

    // The coupled faces of the other boundary faces get no group
    boolList coupledFace(nBoundary, false);
    for (const polyPatch& pp : mesh_.boundaryMesh())
    {
        if (pp.coupled())
        {
            SubList<bool>(coupledFace, pp.size(), pp.offset()) = true;
        }
    }

    auto isSlave = [&](const label celli)
    {
        return alpha[celli] > 0 && alpha[celli] < mergeAlpha_;
    };

    // Link face of every small cell: the open internal or coupled face to a
    // fluid neighbour with more fluid, preferring the direction away from
    // the wall (-n_w) and large wet areas
    scalarField nbrAlpha;
    syncTools::swapBoundaryCellList(mesh_, alpha, nbrAlpha);
    label nOrphan = 0;
    forAll(alpha, celli)
    {
        if (!isSlave(celli)) continue;
        linkFace_[celli] = outwardFace(celli, nbrAlpha);
        if (linkFace_[celli] < 0) ++nOrphan;
    }

    // Groups: every fluid cell that is not linked starts as the master of
    // its own group, the global index of the cell. The links pass the group
    // down the chains (alpha increases along a link: no cycles), also
    // across coupled faces, with the number of hops and the fluid centroid
    // relative to the master's. A chain longer than maxHops stays without
    // a master: independent of the numbering and the decomposition.
    const label maxHops = 10;
    labelList hops(nCells, 0);
    forAll(alpha, celli)
    {
        if (alpha[celli] > 0 && linkFace_[celli] < 0)
        {
            linkGroup_[celli] = cellIndex.toGlobal(celli);
        }
    }
    const surfaceVectorField delta(mesh_.delta());
    while (true)
    {
        bool changed = false;
        bool passChanged = true;
        while (passChanged)
        {
            passChanged = false;
            forAll(linkFace_, celli)
            {
                const label f = linkFace_[celli];
                if (f < 0 || f >= nInternal || linkGroup_[celli] >= 0) continue;
                const label o = (own[f] == celli) ? nei[f] : own[f];
                if (linkGroup_[o] < 0 || hops[o] >= maxHops) continue;
                linkGroup_[celli] = linkGroup_[o];
                hops[celli] = hops[o] + 1;
                linkRel_[celli] = linkRel_[o] + off[celli] - off[o] - (C[o] - C[celli]);
                passChanged = true;
                changed = true;
            }
        }

        labelList nbrGroup, nbrHops;
        vectorField nbrRel, nbrOff;
        syncTools::swapBoundaryCellList(mesh_, linkGroup_, nbrGroup);
        syncTools::swapBoundaryCellList(mesh_, hops, nbrHops);
        syncTools::swapBoundaryCellList(mesh_, linkRel_, nbrRel);
        syncTools::swapBoundaryCellList(mesh_, off, nbrOff);
        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& fvp = mesh_.boundary()[patchi];
            if (!fvp.coupled()) continue;
            const fvsPatchVectorField& dp = delta.boundaryField()[patchi];
            forAll(fvp, i)
            {
                const label f = fvp.start() + i;
                const label celli = fvp.faceCells()[i];
                if (linkFace_[celli] != f || linkGroup_[celli] >= 0) continue;
                const label bf = f - nInternal;
                if (nbrGroup[bf] < 0 || nbrHops[bf] >= maxHops) continue;
                linkGroup_[celli] = nbrGroup[bf];
                hops[celli] = nbrHops[bf] + 1;
                linkRel_[celli] = nbrRel[bf] + off[celli] - nbrOff[bf] - dp[i];
                changed = true;
            }
        }
        if (!returnReduce(changed, orOp<bool>())) break;
    }

    // Masters: cells that are not linked but linked into
    boolList linkedInto(nCells, false);
    label nLinked = 0;
    forAll(linkFace_, celli)
    {
        const label f = linkFace_[celli];
        if (f < 0) continue;
        if (linkGroup_[celli] < 0)
        {
            linkFace_[celli] = -1;
            ++nOrphan;
            continue;
        }
        ++nLinked;
        if (f < nInternal)
        {
            linkedInto[(own[f] == celli) ? nei[f] : own[f]] = true;
        }
        else
        {
            boundaryLink_[f - nInternal] = true;
        }
    }
    boolList linkIn(boundaryLink_);
    syncTools::swapBoundaryFaceList(mesh_, linkIn);
    forAll(linkIn, bf)
    {
        if (coupledFace[bf] && linkIn[bf])
        {
            boundaryLink_[bf] = true;
            linkedInto[faceOwner[nInternal + bf]] = true;
        }
    }
    forAll(linkGroup_, celli)
    {
        if (linkFace_[celli] < 0 && !linkedInto[celli])
        {
            linkGroup_[celli] = -1;
            linkRel_[celli] = Zero;
        }
    }
    syncTools::swapBoundaryCellList(mesh_, linkGroup_, nbrLinkGroup_);
    forAll(nbrLinkGroup_, bf)
    {
        if (!coupledFace[bf]) nbrLinkGroup_[bf] = -1;
    }

    // Slots of the groups on this processor, and the exchange of the group
    // sums with the processors of the masters
    {
        Map<label> slotOf;
        DynamicList<label> groups;
        forAll(linkGroup_, celli)
        {
            const label g = linkGroup_[celli];
            if (g < 0) continue;
            const auto iter = slotOf.cfind(g);
            if (iter.found())
            {
                groupSlot_[celli] = *iter;
            }
            else
            {
                groupSlot_[celli] = groups.size();
                slotOf.insert(g, groups.size());
                groups.append(g);
            }
        }
        slotGroup_.transfer(groups);

        if (Pstream::parRun())
        {
            List<DynamicList<label>> send(Pstream::nProcs());
            forAll(slotGroup_, s)
            {
                const label proci = cellIndex.whichProcID(slotGroup_[s]);
                if (proci != Pstream::myProcNo()) send[proci].append(s);
            }
            PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
            forAll(send, proci)
            {
                if (send[proci].size())
                {
                    UOPstream os(proci, pBufs);
                    os << labelList(UIndirectList<label>(slotGroup_, send[proci]));
                }
                sendSlots_[proci].transfer(send[proci]);
            }
            labelList recvSizes;
            pBufs.finishedSends(recvSizes);
            forAll(recvSizes, proci)
            {
                if (proci == Pstream::myProcNo() || !recvSizes[proci]) continue;
                UIPstream is(proci, pBufs);
                const labelList groupsOf(is);
                labelList& slots = recvSlots_[proci];
                slots.setSize(groupsOf.size());
                forAll(groupsOf, k)
                {
                    const auto iter = slotOf.cfind(groupsOf[k]);
                    if (!iter.found())
                    {
                        FatalErrorInFunction
                            << "Group " << groupsOf[k] << " of processor "
                            << proci << " has no master here"
                            << abort(FatalError);
                    }
                    slots[k] = *iter;
                }
            }
        }
    }

    // Split groups: those with a link over a coupled face, counted on the
    // processor of the master
    {
        labelList nCoupledLinks(slotGroup_.size(), 0);
        forAll(linkFace_, celli)
        {
            if (linkFace_[celli] >= nInternal)
            {
                ++nCoupledLinks[groupSlot_[celli]];
            }
        }
        groupSlotSum(nCoupledLinks);
        DynamicList<label> cells;
        forAll(groupSlot_, celli)
        {
            const label s = groupSlot_[celli];
            if (s < 0 || nCoupledLinks[s] == 0) continue;
            cells.append(celli);
            if (linkFace_[celli] < 0) ++nSplitGroups_;
        }
        splitGroupCells_.transfer(cells);
        reduce(nSplitGroups_, sumOp<label>());
    }

    // Wall distance of the centroid of each merged group (master and all its
    // slaves), from the level sets at the fluid centroids
    bool levelSetOnly = true;
    forAll(bodies_, bodyi)
    {
        if (!bodies_[bodyi].hasLevelSet()) levelSetOnly = false;
    }
    scalarField volSum(alpha*V);
    groupSum(volSum);
    scalarField yGroupSum(nCells, 0), weightSum(nCells, 0);
    if (levelSetOnly)
    {
        // The level set at the fluid centroids, fluid-volume weighted
        const pointField xc(fluidCentroids());
        scalarField wd;
        wallDistance(xc, wd);
        yGroupSum = alpha*V*wd;
        weightSum = volSum;
    }
    else
    {
        // Without a level set only the wall cells know their distance, the
        // first-order alpha V/(2 Awall) of their fluid centroid: each of a
        // group extends its wall plane to the group centroid, weighted by
        // its wall area (the level-set result for a planar wall)
        vectorField momSum(alpha*V*linkRel_);
        groupSum(momSum);
        for (const label celli : wallCells_)
        {
            if (linkGroup_[celli] < 0) continue;
            const vector n(Sw_[celli]/max(Awall_[celli], VSMALL));
            const vector toCentroid
            (
                momSum[celli]/max(volSum[celli], VSMALL) - linkRel_[celli]
            );
            yGroupSum[celli] =
                Awall_[celli]*(dWall_[celli] - (n & toCentroid));
            weightSum[celli] = Awall_[celli];
        }
    }
    groupSum(yGroupSum);
    if (!levelSetOnly) groupSum(weightSum);
    scalarField& dWall = dWall_.primitiveFieldRef();
    scalarField& wallCoeff = wallCoeff_.primitiveFieldRef();
    scalarField& noSlip = noSlipCoeff_.primitiveFieldRef();
    forAll(alpha, celli)
    {
        if (linkGroup_[celli] < 0 || weightSum[celli] <= 0) continue;
        const scalar yGroup = yGroupSum[celli]/weightSum[celli];
        linkWallDist_[celli] = yGroup;

        // The wall of a slave belongs to the merged cell, whose value sits at
        // the merged centroid: wall distance and wall coefficients from it
        if (linkFace_[celli] >= 0 && wallCells_.test(celli))
        {
            dWall[celli] = yGroup;
            const scalar c = Awall_[celli]/(yGroup*V[celli]);
            if (wallCoeff[celli] > 0) wallCoeff[celli] = c;
            if (noSlip[celli] > 0) noSlip[celli] = c;
        }
    }
    dWall_.correctBoundaryConditions();
    wallCoeff_.correctBoundaryConditions();
    noSlipCoeff_.correctBoundaryConditions();

    label nCoupled = 0;
    forAll(linkFace_, celli)
    {
        if (linkFace_[celli] >= nInternal) ++nCoupled;
    }
    reduce(nLinked, sumOp<label>());
    reduce(nCoupled, sumOp<label>());
    reduce(nOrphan, sumOp<label>());
    Info<< "    cell linking: " << nLinked << " cells with alpha < "
        << mergeAlpha_ << " linked to a master";
    if (nCoupled > 0)
    {
        Info<< " (" << nCoupled << " over coupled faces)";
    }
    if (nOrphan > 0)
    {
        Info<< ", " << nOrphan << " without a master (left unlinked)";
    }
    Info<< endl;
}


// * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

void Foam::cutCellGeometry::unionLevelSet
(
    const pointField& pts,
    scalarField& phi
) const
{
    phi.setSize(pts.size());
    phi = GREAT;
    scalarField phib;
    forAll(bodies_, bodyi)
    {
        if (bodies_[bodyi].hasLevelSet())
        {
            bodies_[bodyi].evaluate(pts, phib);
            phi = min(phi, phib);
        }
    }
}


bool Foam::cutCellGeometry::inGroup
(
    const label bodyi,
    const bodyGroup group
) const
{
    const ibmBody& b = bodies_[bodyi];
    if (!b.hasLevelSet()) return false;
    switch (group)
    {
        case STATICBODIES: return !b.moving();
        case MOVINGBODIES: return b.moving();
        default: return true;
    }
}


void Foam::cutCellGeometry::unionLevelSet
(
    const pointField& pts,
    scalarField& phi,
    const bodyGroup group
) const
{
    phi.setSize(pts.size());
    phi = GREAT;
    scalarField phib;
    forAll(bodies_, bodyi)
    {
        if (inGroup(bodyi, group))
        {
            bodies_[bodyi].evaluate(pts, phib);
            phi = min(phi, phib);
        }
    }
}


bool Foam::cutCellGeometry::staticCacheActive() const
{
    if (!cacheStatic_ || !moving_) return false;
    bool hasStatic = false, hasMoving = false;
    forAll(bodies_, bodyi)
    {
        if (inGroup(bodyi, STATICBODIES)) hasStatic = true;
        if (inGroup(bodyi, MOVINGBODIES)) hasMoving = true;
    }
    return hasStatic && hasMoving;
}


Foam::label Foam::cutCellGeometry::cellSide
(
    const label celli,
    const scalarField& pointPhi,
    const scalarField& facePhi,
    const scalarField& cellPhi
) const
{
    const faceList& faces = mesh_.faces();
    label nPos = 0, nNeg = 0;
    (cellPhi[celli] > 0 ? nPos : nNeg)++;
    for (const label facei : mesh_.cells()[celli])
    {
        (facePhi[facei] > 0 ? nPos : nNeg)++;
        const face& f = faces[facei];
        forAll(f, i)
        {
            (pointPhi[f[i]] > 0 ? nPos : nNeg)++;
        }
    }
    if (nPos > 0 && nNeg > 0) return 0;
    return (nPos > 0) ? 1 : -1;
}


void Foam::cutCellGeometry::bisectEdges
(
    const scalarField& pointPhi,
    const labelUList& cutEdges,
    const bodyGroup group,
    pointField& edgeCross
) const
{
    // All edges are bisected together so that the level set is evaluated
    // in batches (matters for STL surfaces)
    const pointField& points = mesh_.points();
    const edgeList& edges = mesh_.edges();

    pointField A(cutEdges.size()), B(cutEdges.size());
    boolList fluidA(cutEdges.size());
    forAll(cutEdges, i)
    {
        const edge& e = edges[cutEdges[i]];
        label a = e[0], b = e[1];
        if (lessPoint(points[b], points[a])) Swap(a, b);
        A[i] = points[a];
        B[i] = points[b];
        fluidA[i] = pointPhi[a] > 0;
    }

    scalarField lo(cutEdges.size(), 0.0), hi(cutEdges.size(), 1.0);
    pointField mid(cutEdges.size());
    scalarField phiMid;
    for (label it = 0; it < nBisect_; ++it)
    {
        forAll(mid, i)
        {
            mid[i] = A[i] + 0.5*(lo[i] + hi[i])*(B[i] - A[i]);
        }
        unionLevelSet(mid, phiMid, group);
        forAll(mid, i)
        {
            if ((phiMid[i] > 0) == fluidA[i])
            {
                lo[i] = 0.5*(lo[i] + hi[i]);
            }
            else
            {
                hi[i] = 0.5*(lo[i] + hi[i]);
            }
        }
    }

    forAll(cutEdges, i)
    {
        edgeCross[cutEdges[i]] = A[i] + 0.5*(lo[i] + hi[i])*(B[i] - A[i]);
    }
}


void Foam::cutCellGeometry::buildStaticCache() const
{
    const pointField& points = mesh_.points();
    const edgeList& edges = mesh_.edges();
    const cellList& cells = mesh_.cells();

    unionLevelSet(points, staticPointPhi_, STATICBODIES);
    unionLevelSet(mesh_.faceCentres(), staticFacePhi_, STATICBODIES);
    unionLevelSet(mesh_.cellCentres(), staticCellPhi_, STATICBODIES);

    DynamicList<label> cutEdges(edges.size()/8);
    forAll(edges, edgei)
    {
        const edge& e = edges[edgei];
        if ((staticPointPhi_[e[0]] > 0) != (staticPointPhi_[e[1]] > 0))
        {
            cutEdges.append(edgei);
        }
    }
    staticEdgeCross_.setSize(edges.size());
    staticEdgeCross_ = point::max;
    bisectEdges(staticPointPhi_, cutEdges, STATICBODIES, staticEdgeCross_);

    staticAlpha_.setSize(cells.size());
    staticAlpha_ = -1;
    staticCentroid_ = mesh_.cellCentres();
    label nCut = 0;
    forAll(cells, celli)
    {
        if
        (
            cellSide(celli, staticPointPhi_, staticFacePhi_, staticCellPhi_)
         == 0
        )
        {
            staticAlpha_[celli] =
                sampleCellFraction(celli, staticCentroid_[celli], STATICBODIES);
            ++nCut;
        }
    }

    staticBodyPhi_.setSize(bodies_.size());
    forAll(bodies_, bodyi)
    {
        staticBodyPhi_[bodyi].clear();
        if (inGroup(bodyi, STATICBODIES))
        {
            bodies_[bodyi].evaluate(mesh_.cellCentres(), staticBodyPhi_[bodyi]);
        }
    }

    // Circumradius of every cell about its centre: of its bounding box,
    // which holds its points, face centres and sub-samples
    cellRadius_.setSize(cells.size());
    forAll(cells, celli)
    {
        const boundBox bb
        (
            cells[celli].points(mesh_.faces(), points),
            false
        );
        const point& C = mesh_.cellCentres()[celli];
        const tmp<pointField> tcorners(bb.points());
        scalar r2 = 0;
        for (const point& p : tcorners())
        {
            r2 = max(r2, magSqr(p - C));
        }
        cellRadius_[celli] = Foam::sqrt(r2);
    }

    staticCacheValid_ = true;
    Info<< "    static bodies cached: "
        << returnReduce(nCut, sumOp<label>()) << " sub-sampled cells, "
        << returnReduce(cutEdges.size(), sumOp<label>()) << " cut edges"
        << endl;
}


void Foam::cutCellGeometry::cellCentreLevelSet(scalarField& phi) const
{
    if (staticCacheActive() && staticCacheValid_)
    {
        scalarField phiM;
        unionLevelSet(mesh_.cellCentres(), phiM, MOVINGBODIES);
        phi = min(staticCellPhi_, phiM);
    }
    else
    {
        unionLevelSet(mesh_.cellCentres(), phi);
    }
}


void Foam::cutCellGeometry::wallDistance
(
    const pointField& pts,
    scalarField& d
) const
{
    forAll(bodies_, bodyi)
    {
        if (!bodies_[bodyi].hasLevelSet())
        {
            FatalErrorInFunction
                << "Body " << bodies_[bodyi].name() << " has no level set: "
                << "the distance to the immersed wall is only available for "
                << "level-set bodies" << exit(FatalError);
        }
    }
    unionLevelSet(pts, d);
}


void Foam::cutCellGeometry::checkSchemes() const
{
    // The second-order viscous fluxes live in the cutCellCorrected schemes:
    // with the switch on and another laplacian or snGrad scheme the wall
    // treatment is only partly second order (V1/F0, still first order), and
    // with the switch off cutCellCorrected is the uncorrected scheme
    const dictionary& schemes = mesh_.schemesDict();
    for (const word sub : {"laplacianSchemes", "snGradSchemes"})
    {
        const dictionary* dictPtr = schemes.findDict(sub);
        if (!dictPtr) continue;
        for (const entry& e : *dictPtr)
        {
            if (e.isDict()) continue;
            OStringStream os;
            os << e.stream();
            const bool cut = os.str().find("cutCellCorrected") != string::npos;
            const bool none = os.str().find("none") != string::npos;
            if (none || cut == secondOrder_) continue;
            WarningInFunction
                << sub << ' ' << e.keyword() << " is " << os.str()
                << (
                       secondOrder_
                     ? ": the second-order wall treatment needs cutCellCorrected"
                     : ": cutCellCorrected without secondOrder is the"
                       " uncorrected scheme"
                   )
                << endl;
        }
    }
}


Foam::tmp<Foam::pointField> Foam::cutCellGeometry::fluidCentroids() const
{
    return mesh_.cellCentres() + centroidOffset_.primitiveField();
}


void Foam::cutCellGeometry::faceFractions
(
    const scalarField& pointPhi,
    const pointField& edgeCross,
    scalarField& theta,
    pointField& wetCentre
) const
{
    const faceList& faces = mesh_.faces();
    const pointField& points = mesh_.points();
    const labelListList& faceEdges = mesh_.faceEdges();
    const edgeList& edges = mesh_.edges();
    const vectorField& Sf = mesh_.faceAreas();

    theta.setSize(faces.size());
    wetCentre = mesh_.faceCentres();

    DynamicList<point> poly(16);

    forAll(faces, facei)
    {
        const face& f = faces[facei];
        const labelList& fe = faceEdges[facei];

        label nFluid = 0;
        forAll(f, i)
        {
            if (pointPhi[f[i]] > 0) ++nFluid;
        }

        if (nFluid == f.size())
        {
            theta[facei] = 1;
            continue;
        }
        if (nFluid == 0)
        {
            theta[facei] = 0;
            continue;
        }

        // Fluid polygon: fluid vertices and the crossing points of cut edges
        poly.clear();
        forAll(f, i)
        {
            const label pi = f[i];
            const label pn = f[f.fcIndex(i)];
            const bool fluidi = pointPhi[pi] > 0;
            const bool fluidn = pointPhi[pn] > 0;

            if (fluidi)
            {
                poly.append(points[pi]);
            }
            if (fluidi != fluidn)
            {
                const label edgei = fe[i];
                if (debug && edges[edgei] != edge(pi, pn))
                {
                    FatalErrorInFunction
                        << "faceEdges ordering assumption violated"
                        << abort(FatalError);
                }
                poly.append(edgeCross[edgei]);
            }
        }

        vector area(Zero);
        vector moment(Zero);
        const point& p0 = poly[0];
        for (label i = 1; i + 1 < poly.size(); ++i)
        {
            const vector a(0.5*((poly[i] - p0) ^ (poly[i+1] - p0)));
            area += a;
            moment += mag(a)*(p0 + poly[i] + poly[i+1])/3.0;
        }

        // Centroid of the wet polygon by fan triangulation, the triangles
        // weighted by their area projected on the face normal
        {
            const vector n(Sf[facei]/max(mag(Sf[facei]), VSMALL));
            scalar sumA = 0;
            vector sumAx(Zero);
            for (label i = 1; i + 1 < poly.size(); ++i)
            {
                const scalar a =
                    0.5*(((poly[i] - p0) ^ (poly[i+1] - p0)) & n);
                sumA += a;
                sumAx += a*(p0 + poly[i] + poly[i+1])/3.0;
            }
            if (mag(sumA) > VSMALL)
            {
                wetCentre[facei] = sumAx/sumA;
            }
        }

        const scalar magSf = mag(Sf[facei]);
        theta[facei] =
            (magSf > VSMALL) ? min(max(mag(area)/magSf, 0.0), 1.0) : 0.0;
        if (mag(area) > VSMALL && wetOffsetAll_.size() == faces.size())
        {
            wetOffsetAll_[facei] = moment/mag(area) - mesh_.faceCentres()[facei];
        }
    }
}


bool Foam::cutCellGeometry::sampleInterpolatedCentroid
(
    const label celli,
    const scalarField& pointPhi,
    point& centroid
) const
{
    const labelList& cp = mesh_.cellPoints()[celli];
    if (cp.size() != 8 || mesh_.cells()[celli].size() != 6) return false;
    const pointField pts(mesh_.points(), cp);
    const boundBox bb(pts, false);
    const vector span(bb.span());
    if
    (
        mag(span.x()*span.y()*span.z() - mesh_.cellVolumes()[celli])
      > 1e-8*mesh_.cellVolumes()[celli]
    )
    {
        return false;
    }

    // Corner values of the trilinear interpolant, indexed by the corner bits
    scalar phiC[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    forAll(pts, i)
    {
        label k = 0;
        for (direction d = 0; d < 3; ++d)
        {
            if (pts[i][d] - bb.min()[d] > 0.5*span[d]) k |= (1 << d);
        }
        phiC[k] = pointPhi[cp[i]];
    }

    const Vector<label>& solD = mesh_.solutionD();
    const label n = nSubSamples_;
    Vector<label> nDir(n, n, n);
    for (direction d = 0; d < 3; ++d)
    {
        if (solD[d] < 0) nDir[d] = 1;
    }

    vector sum(Zero);
    label nFluid = 0;
    for (label i = 0; i < nDir.x(); ++i)
    {
        for (label j = 0; j < nDir.y(); ++j)
        {
            for (label k = 0; k < nDir.z(); ++k)
            {
                const vector t
                (
                    (i + 0.5)/nDir.x(), (j + 0.5)/nDir.y(), (k + 0.5)/nDir.z()
                );
                scalar phi = 0;
                for (label c = 0; c < 8; ++c)
                {
                    scalar w = 1;
                    for (direction d = 0; d < 3; ++d)
                    {
                        w *= (c & (1 << d)) ? t[d] : 1 - t[d];
                    }
                    phi += w*phiC[c];
                }
                if (phi > 0)
                {
                    sum += cmptMultiply(t, span);
                    ++nFluid;
                }
            }
        }
    }
    if (nFluid == 0) return false;

    centroid = bb.min() + sum/nFluid;
    for (direction d = 0; d < 3; ++d)
    {
        if (solD[d] < 0) centroid[d] = mesh_.cellCentres()[celli][d];
    }
    return true;
}


bool Foam::cutCellGeometry::distanceGroup(const bodyGroup group) const
{
    forAll(bodies_, bodyi)
    {
        if (inGroup(bodyi, group) && !bodies_[bodyi].distanceLevelSet())
        {
            return false;
        }
    }
    return true;
}


namespace
{
    // Coordinate d of the sub-sample with index i: the arithmetic of the
    // full sub-sampling, so the hierarchical one evaluates the same points
    inline Foam::scalar sampleCoord
    (
        const Foam::boundBox& bb,
        const Foam::Vector<Foam::label>& nDir,
        const Foam::Vector<Foam::label>& solD,
        const Foam::point& C,
        const Foam::direction d,
        const Foam::label i
    )
    {
        return
            (solD[d] < 0)
          ? C[d]
          : bb.min()[d] + (i + 0.5)/nDir[d]*bb.span()[d];
    }

    struct sampleBox
    {
        Foam::label lo[3];
        Foam::label hi[3];   // exclusive
    };
}


void Foam::cutCellGeometry::classifySamples
(
    const label celli,
    const boundBox& bb,
    const Vector<label>& nDir,
    const bodyGroup group,
    List<signed char>& side
) const
{
    // Margin on the half-diagonal for the round-off of the body transform
    static const scalar REACH = 1 + 1e-6;

    const Vector<label>& solD = mesh_.solutionD();
    const point& C = mesh_.cellCentres()[celli];
    const label ny = nDir.y(), nz = nDir.z();
    side.setSize(nDir.x()*ny*nz);
    side = 0;

    auto fill = [&](const sampleBox& b, const signed char s)
    {
        for (label i = b.lo[0]; i < b.hi[0]; ++i)
        {
            for (label j = b.lo[1]; j < b.hi[1]; ++j)
            {
                for (label k = b.lo[2]; k < b.hi[2]; ++k)
                {
                    side[(i*ny + j)*nz + k] = s;
                }
            }
        }
    };

    DynamicList<sampleBox> boxes(1), next, multi;
    {
        sampleBox all;
        for (direction d = 0; d < 3; ++d)
        {
            all.lo[d] = 0;
            all.hi[d] = nDir[d];
        }
        boxes.append(all);
    }
    DynamicList<point> centres;
    DynamicList<scalar> radii;
    DynamicList<label> leafIdx;
    DynamicList<point> leafPts;
    scalarField phi;

    while (boxes.size())
    {
        centres.clear();
        radii.clear();
        multi.clear();
        next.clear();
        for (const sampleBox& b : boxes)
        {
            if
            (
                b.hi[0] - b.lo[0] == 1
             && b.hi[1] - b.lo[1] == 1
             && b.hi[2] - b.lo[2] == 1
            )
            {
                point p;
                for (direction d = 0; d < 3; ++d)
                {
                    p[d] = sampleCoord(bb, nDir, solD, C, d, b.lo[d]);
                }
                leafIdx.append((b.lo[0]*ny + b.lo[1])*nz + b.lo[2]);
                leafPts.append(p);
                continue;
            }
            point c;
            vector half;
            for (direction d = 0; d < 3; ++d)
            {
                const scalar a = sampleCoord(bb, nDir, solD, C, d, b.lo[d]);
                const scalar z = sampleCoord(bb, nDir, solD, C, d, b.hi[d] - 1);
                c[d] = 0.5*(a + z);
                half[d] = 0.5*(z - a);
            }
            multi.append(b);
            centres.append(c);
            radii.append(mag(half));
        }
        if (centres.size())
        {
            unionLevelSet(pointField(centres), phi, group);
        }
        forAll(multi, bi)
        {
            const sampleBox& b = multi[bi];
            const scalar r = REACH*radii[bi];
            if (phi[bi] > r)
            {
                fill(b, 1);
            }
            else if (phi[bi] < -r)
            {
                fill(b, -1);
            }
            else
            {
                // Split every direction with more than one sample in two
                label m[3][3];
                label nParts[3];
                for (direction d = 0; d < 3; ++d)
                {
                    if (b.hi[d] - b.lo[d] > 1)
                    {
                        m[d][0] = b.lo[d];
                        m[d][1] = (b.lo[d] + b.hi[d])/2;
                        m[d][2] = b.hi[d];
                        nParts[d] = 2;
                    }
                    else
                    {
                        m[d][0] = b.lo[d];
                        m[d][1] = b.hi[d];
                        nParts[d] = 1;
                    }
                }
                for (label pi = 0; pi < nParts[0]; ++pi)
                {
                    for (label pj = 0; pj < nParts[1]; ++pj)
                    {
                        for (label pk = 0; pk < nParts[2]; ++pk)
                        {
                            sampleBox s;
                            s.lo[0] = m[0][pi]; s.hi[0] = m[0][pi + 1];
                            s.lo[1] = m[1][pj]; s.hi[1] = m[1][pj + 1];
                            s.lo[2] = m[2][pk]; s.hi[2] = m[2][pk + 1];
                            next.append(s);
                        }
                    }
                }
            }
        }
        boxes.transfer(next);
    }

    // The single samples, evaluated at the coordinates of the full
    // sub-sampling: the same values, the same signs
    if (leafPts.size())
    {
        unionLevelSet(pointField(leafPts), phi, group);
        forAll(leafIdx, i)
        {
            side[leafIdx[i]] = (phi[i] > 0) ? 1 : -1;
        }
    }
}


Foam::scalar Foam::cutCellGeometry::sampleCellFraction
(
    const label celli,
    point& centroid,
    const bodyGroup group
) const
{
    const pointField cellPts
    (
        mesh_.cells()[celli].points(mesh_.faces(), mesh_.points())
    );
    const boundBox bb(cellPts, false);
    const Vector<label>& solD = mesh_.solutionD();
    const point& C = mesh_.cellCentres()[celli];

    // Axis-aligned hexahedron: every point of the bounding box is inside
    bool aligned = false;
    if (cellPts.size() == 8 && mesh_.cells()[celli].size() == 6)
    {
        const vector span(bb.span());
        const scalar bbVol = span.x()*span.y()*span.z();
        aligned = mag(bbVol - mesh_.cellVolumes()[celli])
                < 1e-8*mesh_.cellVolumes()[celli];
    }

    const label n = nSubSamples_;
    Vector<label> nDir(n, n, n);
    for (direction d = 0; d < 3; ++d)
    {
        if (solD[d] < 0) nDir[d] = 1;   // empty direction: one sample
    }

    centroid = C;

    if (hierSampling_ && distanceGroup(group))
    {
        // The signs of all sub-samples by hierarchical classification, then
        // the sums of the full sub-sampling in its order
        List<signed char> side;
        classifySamples(celli, bb, nDir, group, side);
        label nIn = 0, nFluid = 0;
        vector sum(Zero);
        for (label i = 0; i < nDir.x(); ++i)
        {
            for (label j = 0; j < nDir.y(); ++j)
            {
                for (label k = 0; k < nDir.z(); ++k)
                {
                    point p
                    (
                        bb.min().x() + (i + 0.5)/nDir.x()*bb.span().x(),
                        bb.min().y() + (j + 0.5)/nDir.y()*bb.span().y(),
                        bb.min().z() + (k + 0.5)/nDir.z()*bb.span().z()
                    );
                    for (direction d = 0; d < 3; ++d)
                    {
                        if (solD[d] < 0) p[d] = C[d];
                    }
                    if
                    (
                        aligned
                     || mesh_.pointInCell(p, celli, polyMesh::FACE_PLANES)
                    )
                    {
                        ++nIn;
                        if (side[(i*nDir.y() + j)*nDir.z() + k] > 0)
                        {
                            ++nFluid;
                            sum += p;
                        }
                    }
                }
            }
        }
        if (nIn == 0)
        {
            return 1.0;
        }
        if (nFluid > 0)
        {
            centroid = sum/nFluid;
        }
        return scalar(nFluid)/nIn;
    }

    DynamicList<point> samples(nDir.x()*nDir.y()*nDir.z());
    for (label i = 0; i < nDir.x(); ++i)
    {
        for (label j = 0; j < nDir.y(); ++j)
        {
            for (label k = 0; k < nDir.z(); ++k)
            {
                point p
                (
                    bb.min().x() + (i + 0.5)/nDir.x()*bb.span().x(),
                    bb.min().y() + (j + 0.5)/nDir.y()*bb.span().y(),
                    bb.min().z() + (k + 0.5)/nDir.z()*bb.span().z()
                );
                for (direction d = 0; d < 3; ++d)
                {
                    if (solD[d] < 0) p[d] = C[d];
                }
                if (aligned || mesh_.pointInCell(p, celli, polyMesh::FACE_PLANES))
                {
                    samples.append(p);
                }
            }
        }
    }

    if (samples.empty())
    {
        return 1.0;
    }

    scalarField phi;
    unionLevelSet(pointField(samples), phi, group);

    label nFluid = 0;
    vector sum(Zero);
    forAll(phi, i)
    {
        if (phi[i] > 0)
        {
            ++nFluid;
            sum += samples[i];
        }
    }
    if (nFluid > 0)
    {
        centroid = sum/nFluid;
    }
    return scalar(nFluid)/samples.size();
}


void Foam::cutCellGeometry::calcFromLevelSet
(
    scalarField& alpha,
    scalarField& theta,
    pointField& wetCentre,
    labelList& cellBody,
    vectorField& centroidOffset,
    scalarField& pointPhi,
    pointField& edgeCross
) const
{
    const pointField& points = mesh_.points();
    const edgeList& edges = mesh_.edges();
    const cellList& cells = mesh_.cells();

    // Wall-clock times of the sub-stages (reported for moving geometry)
    clockTime clock;
    FixedList<scalar, 5> times(Zero);

    // Static-body cache: the static part once, the moving bodies here
    const bool useCache = staticCacheActive();
    if (useCache && !staticCacheValid_)
    {
        buildStaticCache();
    }

    // Level set at the mesh points, face centres and cell centres
    scalarField facePhi, cellPhi;
    scalarField movingPointPhi, movingCellPhi;
    if (useCache)
    {
        scalarField movingFacePhi;
        unionLevelSet(points, movingPointPhi, MOVINGBODIES);
        unionLevelSet(mesh_.faceCentres(), movingFacePhi, MOVINGBODIES);
        unionLevelSet(mesh_.cellCentres(), movingCellPhi, MOVINGBODIES);
        pointPhi = min(staticPointPhi_, movingPointPhi);
        facePhi = min(staticFacePhi_, movingFacePhi);
        cellPhi = min(staticCellPhi_, movingCellPhi);
    }
    else
    {
        unionLevelSet(points, pointPhi);
        unionLevelSet(mesh_.faceCentres(), facePhi);
        unionLevelSet(mesh_.cellCentres(), cellPhi);
    }

    times[0] = clock.timeIncrement();

    // Crossing point on every cut edge, by bisection of the exact level set.
    // With the cache, an edge the moving bodies cannot reach (their level
    // set above the edge length at an end; distances, so positive along the
    // whole edge) is cut by the static bodies alone, at their crossing.
    // A small relative margin covers the round-off of the body transform.
    static const scalar REACH = 1 + 1e-6;
    edgeCross.setSize(edges.size());
    edgeCross = point::max;
    DynamicList<label> cutEdges(edges.size()/8);
    label nCutEdges = 0, nReused = 0;
    forAll(edges, edgei)
    {
        const edge& e = edges[edgei];
        if ((pointPhi[e[0]] > 0) == (pointPhi[e[1]] > 0)) continue;
        ++nCutEdges;
        if (useCache)
        {
            const scalar L = REACH*e.mag(points);
            if
            (
                (movingPointPhi[e[0]] > L || movingPointPhi[e[1]] > L)
             && staticEdgeCross_[edgei] != point::max
            )
            {
                edgeCross[edgei] = staticEdgeCross_[edgei];
                ++nReused;
                continue;
            }
        }
        cutEdges.append(edgei);
    }
    bisectEdges(pointPhi, cutEdges, ALLBODIES, edgeCross);
    times[1] = clock.timeIncrement();

    // Face fractions
    faceFractions(pointPhi, edgeCross, theta, wetCentre);
    times[2] = clock.timeIncrement();

    // Cell fractions. A cell is a candidate for cutting if its vertices, face
    // centres or centre do not all lie on the same side; otherwise it is
    // wholly fluid or wholly solid (a body smaller than a cell is not
    // resolved, as in the reference implementation).
    alpha.setSize(cells.size());
    cellBody.setSize(cells.size(), -1);

    // Samples per cell (an empty direction has one)
    label nSamples = 1;
    for (direction d = 0; d < 3; ++d)
    {
        if (mesh_.solutionD()[d] > 0) nSamples *= nSubSamples_;
    }
    const vectorField& faceAreas = mesh_.faceAreas();

    label nCandidates = 0, nMissed = 0, nCached = 0;
    forAll(cells, celli)
    {
        const cell& c = cells[celli];
        const label side = cellSide(celli, pointPhi, facePhi, cellPhi);

        if (side == 0)
        {
            point xc;
            ++nCandidates;
            // With the cache, a cell the moving bodies cannot reach (their
            // level set at its centre above its circumradius) has the
            // fraction and centroid of the static bodies
            if
            (
                useCache
             && movingCellPhi[celli] > REACH*cellRadius_[celli]
             && staticAlpha_[celli] >= 0
            )
            {
                alpha[celli] = staticAlpha_[celli];
                xc = staticCentroid_[celli];
                ++nCached;
            }
            else
            {
                alpha[celli] = sampleCellFraction(celli, xc);
            }

            // A cell with a wet face holds fluid (the face fractions are
            // exact up to the bisection), also when its fluid is smaller than
            // a sample and no sample falls into it. Without fluid it would be
            // closed below and the wall area and fluid volume would depend on
            // the sub-sampling (G1, finding 4). Its fraction is then below
            // one sample: half of one; its fluid centroid the mean of the wet
            // face centres, weighted by their wet areas.
            if (alpha[celli] <= 0)
            {
                vector sum(Zero);
                scalar wSum = 0;
                for (const label facei : c)
                {
                    const scalar w = theta[facei]*mag(faceAreas[facei]);
                    if (w > 0)
                    {
                        sum += w*wetCentre[facei];
                        wSum += w;
                    }
                }
                if (wSum > 0)
                {
                    alpha[celli] = 0.5/nSamples;
                    xc = sum/wSum;
                    ++nMissed;
                }
            }
            centroidOffset[celli] = xc - mesh_.cellCentres()[celli];
        }
        else
        {
            alpha[celli] = (side > 0) ? 1.0 : 0.0;
        }
    }

    times[3] = clock.timeIncrement();

    // Owning body of every cut or solid cell: the body closest to its centre
    List<scalarField> bodyPhi(bodies_.size());
    forAll(bodies_, bodyi)
    {
        if (useCache && inGroup(bodyi, STATICBODIES))
        {
            bodyPhi[bodyi] = staticBodyPhi_[bodyi];
        }
        else if (bodies_[bodyi].hasLevelSet())
        {
            bodies_[bodyi].evaluate(mesh_.cellCentres(), bodyPhi[bodyi]);
        }
    }
    forAll(cells, celli)
    {
        scalar phiMin = GREAT;
        label nearest = -1;
        forAll(bodies_, bodyi)
        {
            if
            (
                bodies_[bodyi].hasLevelSet()
             && bodyPhi[bodyi][celli] < phiMin
            )
            {
                phiMin = bodyPhi[bodyi][celli];
                nearest = bodyi;
            }
        }
        nearestBody_[celli] = nearest;
        if (alpha[celli] < 1.0)
        {
            cellBody[celli] = nearest;
        }
    }

    times[4] = clock.timeIncrement();

    reduce(nCandidates, sumOp<label>());
    reduce(nMissed, sumOp<label>());
    reduce(nCutEdges, sumOp<label>());
    Info<< "    level-set bodies: " << nCutEdges
        << " cut edges, " << nCandidates << " sub-sampled cells, " << nMissed
        << " with fluid smaller than a sample";
    if (useCache)
    {
        Info<< " (static-body cache: " << returnReduce(nCached, sumOp<label>())
            << " cells, " << returnReduce(nReused, sumOp<label>())
            << " edges reused)";
    }
    Info<< endl;

    if (moving_)
    {
        // Max and mean over the ranks: their ratio is the load imbalance
        // (the cut cells of a moving body sit on a few ranks)
        static const char* names[5] =
            {"level set at points", "edges", "faces", "cells", "owners"};
        Info<< "    level-set stage [s, max/mean over ranks]:";
        for (label k = 0; k < 5; ++k)
        {
            const scalar tMax = returnReduce(times[k], maxOp<scalar>());
            const scalar tMean =
                returnReduce(times[k], sumOp<scalar>())/Pstream::nProcs();
            Info<< (k ? ", " : " ") << names[k] << ' ' << tMax << '/' << tMean;
        }
        Info<< endl;
    }
}


void Foam::cutCellGeometry::calcFromAlphaField
(
    const alphaFieldBody& body,
    scalarField& alpha,
    scalarField& theta,
    pointField& wetCentre,
    scalarField& pointPhi,
    pointField& edgeCross
) const
{
    const tmp<volScalarField> talpha(body.alpha());
    alpha = talpha().primitiveField();

    // Face fractions from the 0.5 iso-surface of the point-interpolated
    // fraction, with linear interpolation along the edges
    const pointScalarField pAlpha
    (
        volPointInterpolation::New(mesh_).interpolate(talpha())
    );
    pointPhi = pAlpha.primitiveField() - 0.5;

    const edgeList& edges = mesh_.edges();
    const pointField& points = mesh_.points();
    edgeCross.setSize(edges.size());
    edgeCross = point::max;
    forAll(edges, edgei)
    {
        const edge& e = edges[edgei];
        const scalar fa = pointPhi[e[0]];
        const scalar fb = pointPhi[e[1]];
        if ((fa > 0) != (fb > 0))
        {
            const scalar t = fa/(fa - fb);
            edgeCross[edgei] = points[e[0]] + t*(points[e[1]] - points[e[0]]);
        }
    }

    faceFractions(pointPhi, edgeCross, theta, wetCentre);
}


void Foam::cutCellGeometry::calcDerived()
{
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    const surfaceScalarField& deltaCoeffs = mesh_.deltaCoeffs();
    const scalarField& V = mesh_.V();

    // Wall area vector from closure, Sw = -sum_f theta_f Sf
    vectorField& Sw = Sw_.primitiveFieldRef();
    scalarField& blank = blankCoeff_.primitiveFieldRef();
    Sw = Zero;
    blank = 0;

    forAll(own, facei)
    {
        const vector tSf(theta_[facei]*Sf[facei]);
        Sw[own[facei]] -= tSf;
        Sw[nei[facei]] += tSf;
        const scalar b = magSf[facei]*deltaCoeffs[facei];
        blank[own[facei]] += b;
        blank[nei[facei]] += b;
    }
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];
        const labelUList& fc = fvp.faceCells();
        const fvsPatchScalarField& tp = theta_.boundaryField()[patchi];
        const fvsPatchVectorField& Sfp = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& magSfp = magSf.boundaryField()[patchi];
        const fvsPatchScalarField& dcp = deltaCoeffs.boundaryField()[patchi];
        forAll(fvp, i)
        {
            Sw[fc[i]] -= tp[i]*Sfp[i];
            blank[fc[i]] += magSfp[i]*dcp[i];
        }
    }

    // Cell sets
    solidCells_.reset();
    cutCells_.reset();
    wallCells_.reset();

    const scalarField& alpha = alpha_.primitiveField();
    scalarField& Awall = Awall_.primitiveFieldRef();
    scalarField& dWall = dWall_.primitiveFieldRef();
    scalarField& noSlip = noSlipCoeff_.primitiveFieldRef();

    // Level set at the fluid centroids: the wall distance of the second-order
    // treatment (level-set bodies; alpha-field bodies use their edge crossings)
    bool haveCentroidLevelSet = false;
    if (secondOrder_)
    {
        label nLevelSet = 0;
        forAll(bodies_, bodyi)
        {
            if (bodies_[bodyi].hasLevelSet()) ++nLevelSet;
        }
        if (nLevelSet > 0)
        {
            unionLevelSet
            (
                mesh_.cellCentres() + centroidOffset_.primitiveField(),
                centroidLevelSet_
            );
            haveCentroidLevelSet = true;
        }
    }

    fluidVolume_ = 0;
    forAll(alpha, celli)
    {
        Awall[celli] = mag(Sw[celli]);
        fluidVolume_ += alpha[celli]*V[celli];

        if (alpha[celli] <= 0)
        {
            solidCells_.set(celli);
            blank[celli] /= V[celli];
        }
        else
        {
            blank[celli] = 0;
            if (alpha[celli] < 1) cutCells_.set(celli);
        }

        // A cell is a wall cell if its wall segment is not negligible against
        // its faces: cut cells, and fluid cells facing a closed face
        const scalar Aref = pow(V[celli], 2.0/3.0);
        if (alpha[celli] > 0 && Awall[celli] > 1e-10*Aref)
        {
            wallCells_.set(celli);
            dWall[celli] = alpha[celli]*V[celli]/(2*Awall[celli]);
            const label bodyi = cellBody_[celli];
            const bool alphaFieldCell =
                bodyi >= 0 && !bodies_[bodyi].hasLevelSet();
            if (haveCentroidLevelSet && !alphaFieldCell)
            {
                // distance of the fluid centroid to the wall, floored at a
                // small fraction of the cell size
                dWall[celli] =
                    max(centroidLevelSet_[celli], 1e-3*cbrt(V[celli]));
            }
            else if (secondOrder_ && alphaFieldCell)
            {
                // alpha-field body: distance of the fluid centroid to the
                // wall plane through the mean of the cell's edge crossings
                point xw(Zero);
                label nCross = 0;
                for (const label edgei : mesh_.cellEdges()[celli])
                {
                    if (edgeCrossGeo_[edgei] == point::max) continue;
                    xw += edgeCrossGeo_[edgei];
                    ++nCross;
                }
                if (nCross > 0)
                {
                    xw /= nCross;
                    const point xc
                    (
                        mesh_.cellCentres()[celli]
                      + centroidOffset_.primitiveField()[celli]
                    );
                    dWall[celli] = max
                    (
                        (Sw[celli]/Awall[celli]) & (xw - xc),
                        1e-3*cbrt(V[celli])
                    );
                }
            }
        }
        else
        {
            dWall[celli] = GREAT;
        }
    }
    reduce(fluidVolume_, sumOp<scalar>());

    // Owning body of wall cells that are not cut (fluid cells behind a
    // closed face): inherit it from the solid/cut neighbour across the face
    bodyIndex_.primitiveFieldRef() = scalarField(cellBody_.size(), -1);
    forAll(cellBody_, celli)
    {
        bodyIndex_[celli] = cellBody_[celli];
    }
    bodyIndex_.correctBoundaryConditions();

    forAll(own, facei)
    {
        if (theta_[facei] <= 0)
        {
            const label o = own[facei], n = nei[facei];
            if (cellBody_[o] < 0 && cellBody_[n] >= 0) cellBody_[o] = cellBody_[n];
            if (cellBody_[n] < 0 && cellBody_[o] >= 0) cellBody_[n] = cellBody_[o];
        }
    }
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];
        if (!fvp.coupled()) continue;
        const labelUList& fc = fvp.faceCells();
        const fvsPatchScalarField& tp = theta_.boundaryField()[patchi];
        const scalarField nbrBody
        (
            bodyIndex_.boundaryField()[patchi].patchNeighbourField()
        );
        forAll(fvp, i)
        {
            if (tp[i] <= 0 && cellBody_[fc[i]] < 0 && nbrBody[i] >= 0)
            {
                cellBody_[fc[i]] = label(nbrBody[i] + 0.5);
            }
        }
    }
    // Wall cells still without a body (a cell whose fluid fraction rounds
    // to one while one of its faces is partly cut): the nearest body, or
    // the only body if there is no level set
    for (const label celli : wallCells_)
    {
        if (cellBody_[celli] < 0)
        {
            cellBody_[celli] =
                (nearestBody_[celli] >= 0) ? nearestBody_[celli]
              : (bodies_.size() == 1 ? 0 : -1);
        }
    }

    forAll(cellBody_, celli)
    {
        bodyIndex_[celli] = cellBody_[celli];
    }
    bodyIndex_.correctBoundaryConditions();

    // Dirichlet wall coefficient, Awall/(dWall V) = 2 Awall^2/(alpha V^2), and
    // its no-slip mask
    scalarField& wallCoeff = wallCoeff_.primitiveFieldRef();
    wallCoeff = 0;
    noSlip = 0;
    for (const label celli : wallCells_)
    {
        wallCoeff[celli] = Awall[celli]/(dWall[celli]*V[celli]);
        const label bodyi = cellBody_[celli];
        const bool isNoSlip =
            (bodyi < 0) || (bodies_[bodyi].wall() == ibmBody::NOSLIP);
        if (isNoSlip)
        {
            noSlip[celli] = wallCoeff[celli];
        }
    }

    Sw_.correctBoundaryConditions();
    Awall_.correctBoundaryConditions();
    dWall_.correctBoundaryConditions();
    wallCoeff_.correctBoundaryConditions();
    noSlipCoeff_.correctBoundaryConditions();
    blankCoeff_.correctBoundaryConditions();
}


void Foam::cutCellGeometry::calcValueOffset()
{
    valueOffset_ == centroidOffset_;
    if (mergeAlpha_ <= 0) return;

    // Centroid of every merged group: the shared value of a master and its
    // slaves belongs there, so the fluxes through the outer faces of the
    // group see it at one point and stay exact for a linear field. From
    // the fluid centroids relative to the master's (linkRel_): the same
    // across processors and cyclics.
    const scalarField& alpha = alpha_.primitiveField();
    const scalarField& V = mesh_.V();
    const vectorField& offFluid = centroidOffset_.primitiveField();
    scalarField volSum(alpha*V);
    vectorField momSum(alpha*V*linkRel_);
    groupSum(volSum);
    groupSum(momSum);
    vectorField& off = valueOffset_.primitiveFieldRef();
    forAll(linkGroup_, celli)
    {
        if (linkGroup_[celli] < 0) continue;
        off[celli] =
            offFluid[celli] - linkRel_[celli] + momSum[celli]/volSum[celli];
    }
    valueOffset_.correctBoundaryConditions();

    // The wall of a master that is a cut cell itself also sees the value at
    // the group centroid (calcLinks moves the walls of the slaves only)
    scalarField& dWall = dWall_.primitiveFieldRef();
    scalarField& wallCoeff = wallCoeff_.primitiveFieldRef();
    scalarField& noSlip = noSlipCoeff_.primitiveFieldRef();
    forAll(linkGroup_, celli)
    {
        if (linkGroup_[celli] < 0 || linkFace_[celli] >= 0) continue;
        if (!wallCells_.test(celli)) continue;
        dWall[celli] = linkWallDist_[celli];
        const scalar c = Awall_[celli]/(dWall[celli]*V[celli]);
        if (wallCoeff[celli] > 0) wallCoeff[celli] = c;
        if (noSlip[celli] > 0) noSlip[celli] = c;
    }
    dWall_.correctBoundaryConditions();
    wallCoeff_.correctBoundaryConditions();
    noSlipCoeff_.correctBoundaryConditions();
}


bool Foam::cutCellGeometry::insideGroup(const label facei) const
{
    if (mergeAlpha_ <= 0) return false;
    if (facei < mesh_.nInternalFaces())
    {
        const label g = linkGroup_[mesh_.owner()[facei]];
        return g >= 0 && g == linkGroup_[mesh_.neighbour()[facei]];
    }
    const label g = linkGroup_[mesh_.faceOwner()[facei]];
    return g >= 0 && g == nbrLinkGroup_[facei - mesh_.nInternalFaces()];
}


void Foam::cutCellGeometry::calcCentroidGeometry()
{
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    const surfaceVectorField delta(mesh_.delta());
    const surfaceScalarField& meshDC = mesh_.nonOrthDeltaCoeffs();
    const vectorField& off = valueOffset_.primitiveField();
    const vectorField& offOwn = centroidOffset_.primitiveField();
    label nGroupFaces = 0, nClamped = 0;
    scalar minNd = GREAT, maxSkew = 0;

    dCentroid_ = delta;
    cutDeltaCoeffs_ = meshDC;
    faceValueOffset_ = dimensionedVector(dimLength, Zero);

    // The linear interpolate w p_O + (1-w) p_N of values at the centroids is
    // the value at m_f = x_c,O + (1-w) d_c, not at the wet-face centroid
    // x'_f: the gradient corrects it by g_f.(x'_f - m_f) on every open face
    // of a cell with a centroid offset, and on every partially wet face
    const surfaceScalarField& weights = mesh_.weights();
    const volVectorField& C = mesh_.C();
    const surfaceVectorField& Cf = mesh_.Cf();
    auto needsCorrection = [&](const scalar t, const vector& oO, const vector& oN)
    {
        return t > 0 && (t < 1 || magSqr(oO) > 0 || magSqr(oN) > 0);
    };

    forAll(own, facei)
    {
        const label o = own[facei], n = nei[facei];
        if (theta_[facei] <= 0 || solidCells_.test(o) || solidCells_.test(n))
        {
            continue;
        }
        // Inside a group both values sit at the group centroid; the fluxes
        // there cancel in the merged equation, keep the fluid centroids
        const bool inside = insideGroup(facei);
        const vector& oO = inside ? offOwn[o] : off[o];
        const vector& oN = inside ? offOwn[n] : off[n];
        const vector d(delta[facei] + oN - oO);
        const scalar nd = (Sf[facei] & d)/magSf[facei];
        dCentroid_[facei] = d;
        cutDeltaCoeffs_[facei] = 1.0/max(nd, 0.1/meshDC[facei]);
        if (needsCorrection(theta_[facei], oO, oN))
        {
            const point m(C[o] + oO + (1 - weights[facei])*d);
            faceValueOffset_[facei] = Cf[facei] + wetFaceOffset_[facei] - m;
        }
        if (!inside && (oO != offOwn[o] || oN != offOwn[n]))
        {
            ++nGroupFaces;
            minNd = min(minNd, nd*meshDC[facei]);
            maxSkew = max(maxSkew, mag(d)/max(nd, VSMALL));
            if (nd < 0.1/meshDC[facei]) ++nClamped;
        }
    }
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];
        if (!fvp.coupled()) continue;
        const labelUList& fc = fvp.faceCells();
        const vectorField offN
        (
            valueOffset_.boundaryField()[patchi].patchNeighbourField()
        );
        const vectorField offOwnN
        (
            centroidOffset_.boundaryField()[patchi].patchNeighbourField()
        );
        const fvsPatchScalarField& tp = theta_.boundaryField()[patchi];
        const fvsPatchVectorField& deltap = delta.boundaryField()[patchi];
        const fvsPatchVectorField& Sfp = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& magSfp = magSf.boundaryField()[patchi];
        const fvsPatchScalarField& dcp = meshDC.boundaryField()[patchi];
        fvsPatchVectorField& dp = dCentroid_.boundaryFieldRef()[patchi];
        fvsPatchScalarField& cp = cutDeltaCoeffs_.boundaryFieldRef()[patchi];
        const fvsPatchScalarField& wp = weights.boundaryField()[patchi];
        const fvsPatchVectorField& Cfp = Cf.boundaryField()[patchi];
        const fvsPatchVectorField& wop = wetFaceOffset_.boundaryField()[patchi];
        fvsPatchVectorField& fvop = faceValueOffset_.boundaryFieldRef()[patchi];
        forAll(fvp, i)
        {
            if (tp[i] <= 0 || solidCells_.test(fc[i])) continue;
            const bool inside = insideGroup(fvp.start() + i);
            const vector& oO = inside ? offOwn[fc[i]] : off[fc[i]];
            const vector& oN = inside ? offOwnN[i] : offN[i];
            const vector d(deltap[i] + oN - oO);
            const scalar nd = (Sfp[i] & d)/magSfp[i];
            dp[i] = d;
            cp[i] = 1.0/max(nd, 0.1/dcp[i]);
            if (needsCorrection(tp[i], oO, oN))
            {
                const point m(C[fc[i]] + oO + (1 - wp[i])*d);
                fvop[i] = Cfp[i] + wop[i] - m;
            }
        }
    }
    if (mergeAlpha_ > 0)
    {
        reduce(nGroupFaces, sumOp<label>());
        reduce(nClamped, sumOp<label>());
        reduce(minNd, minOp<scalar>());
        reduce(maxSkew, maxOp<scalar>());
        Info<< "    group faces: " << nGroupFaces << ", min n.d/h " << minNd
            << ", max |d|/n.d " << maxSkew << ", clamped " << nClamped << endl;
    }
}


void Foam::cutCellGeometry::calcWallPieces()
{
    // The wall of a cell is split into planar pieces k with area vector S_k
    // (out of the fluid) and centroid x_k, so that the wall pressure force
    // sum_k p(x_k) S_k is exact for a linear pressure:
    //  - the wall polygon(s) of the geometric cut: the trace of the wall on
    //    every cut face of the cell runs from one edge crossing to the next,
    //    and the traces close into loops; each loop is fan-triangulated about
    //    its mean point (in 2-D the loop is the chord times the depth);
    //  - the geometric wet part of every face closed by the post-processing
    //    (absorption of the neighbour, coupled-face sync), area vector
    //    (thetaGeo - theta) Sf at the wet-face centroid.
    // By closure of the geometric fluid polyhedron, sum_k S_k = Sw.
    const faceList& faces = mesh_.faces();
    const cellList& cells = mesh_.cells();
    const labelList& own = mesh_.faceOwner();
    const labelListList& faceEdges = mesh_.faceEdges();
    const vectorField& Sf = mesh_.faceAreas();
    const vectorField& Sw = Sw_.primitiveField();
    const vectorField xc(mesh_.cellCentres() + valueOffset_.primitiveField());

    // Face fractions after post-processing over all faces; -1 on faces of
    // empty patches, which are not part of the finite-volume boundary
    scalarField thetaAll(mesh_.nFaces(), -1);
    SubField<scalar>(thetaAll, mesh_.nInternalFaces()) = theta_.primitiveField();
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];
        SubField<scalar>(thetaAll, fvp.size(), fvp.start()) =
            theta_.boundaryField()[patchi];
    }

    wallMoment_.setSize(mesh_.nCells());
    wallMoment_ = Zero;
    wallPieceSum_.setSize(mesh_.nCells());
    wallPieceSum_ = Zero;
    wallTorqueMoment_.setSize(mesh_.nCells());
    wallTorqueMoment_ = Zero;

    label nFallback = 0;
    scalar maxMismatch = 0;
    Map<label> next(16);
    DynamicList<point> loop(16);
    DynamicList<label> crossings(8);
    DynamicList<bool> isExit(8);

    for (const label celli : wallCells_)
    {
        tensor T(Zero);
        tensor Q(Zero);
        vector Ssum(Zero);
        bool broken = false;

        // Directed wall edges between edge crossings. The fluid polygon of a
        // face, in face order (anticlockwise about Sf), runs through the
        // wall trace from the crossing X where the boundary leaves the fluid
        // to the next crossing Y; the wall polygon, anticlockwise about its
        // own normal (into the body), runs through it the other way for the
        // owner, and the same way for the neighbour (Sf points into it).
        next.clear();
        for (const label facei : cells[celli])
        {
            const face& f = faces[facei];
            const labelList& fe = faceEdges[facei];
            const bool isOwner = (own[facei] == celli);

            // Crossings in face order, flagged as exits (fluid to solid);
            // states alternate, so each exit is followed by an entry
            crossings.clear();
            isExit.clear();
            forAll(f, i)
            {
                const bool fluidi = pointPhiGeo_[f[i]] > 0;
                const bool fluidn = pointPhiGeo_[f[f.fcIndex(i)]] > 0;
                if (fluidi != fluidn)
                {
                    crossings.append(fe[i]);
                    isExit.append(fluidi);
                }
            }
            forAll(crossings, k)
            {
                if (!isExit[k]) continue;
                const label k1 = crossings.fcIndex(k);
                if (isExit[k1]) { broken = true; continue; }
                const label X = crossings[k], Y = crossings[k1];
                if (isOwner) next.set(Y, X); else next.set(X, Y);
            }

            // Geometric wet part closed by the post-processing
            if (thetaAll[facei] >= 0)
            {
                const scalar dTheta = thetaGeo_[facei] - thetaAll[facei];
                if (dTheta > SMALL)
                {
                    const vector Sk((isOwner ? 1 : -1)*dTheta*Sf[facei]);
                    const vector rk(wetCentreGeo_[facei] - xc[celli]);
                    Ssum += Sk;
                    T += rk*Sk;
                    Q += rk*(rk ^ Sk);
                }
            }
        }

        // Close the directed edges into loops
        labelHashSet used(2*next.size());
        forAllConstIters(next, iter)
        {
            if (broken) break;
            const label start = iter.key();
            if (used.found(start)) continue;
            loop.clear();
            label e = start;
            do
            {
                used.insert(e);
                loop.append(edgeCrossGeo_[e]);
                const auto it = next.cfind(e);
                if (!it.found()) { broken = true; break; }
                e = it.val();
            } while (e != start && loop.size() <= next.size());
            if (broken || e != start) { broken = true; break; }

            point m(Zero);
            for (const point& x : loop) m += x;
            m /= loop.size();
            forAll(loop, i)
            {
                const point& a = loop[i];
                const point& b = loop[loop.fcIndex(i)];
                const vector A(0.5*((a - m) ^ (b - m)));
                const vector rk((m + a + b)/3.0 - xc[celli]);
                Ssum += A;
                T += rk*A;
                Q += rk*(rk ^ A);
            }
        }

        const scalar mismatch = mag(Ssum - Sw[celli])/max(Awall_[celli], VSMALL);
        if (broken || mismatch > 1e-6)
        {
            // Fall back on the wall value at the foot of the centroid normal
            ++nFallback;
            T = (dWall_[celli]*Sw[celli]/max(Awall_[celli], VSMALL))*Sw[celli];
            Q = Zero;
            Ssum = Sw[celli];
        }
        else
        {
            maxMismatch = max(maxMismatch, mismatch);
        }
        wallMoment_[celli] = T;
        wallPieceSum_[celli] = Ssum;
        wallTorqueMoment_[celli] = Q;
    }

    reduce(nFallback, sumOp<label>());
    reduce(maxMismatch, maxOp<scalar>());
    Info<< "    wall pieces: max |sum_k S_k - Sw|/|Sw| " << maxMismatch
        << ", cells on the fallback (centroid normal foot) " << nFallback
        << endl;
}


namespace
{
    // A face-gradient correction g_f.(n - d/(n.d)) with |n - d/(n.d)| below
    // this is round-off of an orthogonal face; the correction is dropped there
    constexpr Foam::scalar snGradTol = 1e-10;
}


void Foam::cutCellGeometry::calcZeroShear()
{
    // Zero tangential stress on a free-slip wall (doc/freeSlipPlan.md). For
    // the velocity w = u - u_b relative to a rigid body motion u_b = U +
    // omega x r (whose stress vanishes), w.n = 0 on the wall and the stress
    // is zero if P dw/dn = H.w_t, H = grad n the curvature tensor of the
    // wall (n into the fluid; the cylinder: dw_theta/dr = w_theta/R). With
    // du_b/dn = omega x n, the wall flux of the Laplacian form of the
    // momentum equation is nu (H.w_t + omega x n) per unit area.
    //
    // H from the Hessian of the body's level set at the fluid centroid of
    // the cell (second differences, step 1e-3 of the cell size), projected
    // on the tangent plane, K = P grad(grad phi) P/|grad phi|, and moved to
    // the wall along the normal per principal direction,
    // kappa_w = kappa/(1 - phi kappa) (exact for a distance function).
    // The tangential wall value of w follows from the value at the value
    // point, a distance d = dWall from the wall, by the same condition:
    // w_t,P = w_t,w + d H_w.w_t,w, so w_t,w = (I + d H_w)^-1 w_t,P. Under-
    // resolved curvature (|phi kappa| or |d kappa_w| beyond 1/2) is limited
    // so that the friction stays below the no-slip coefficient Awall/(d V).
    // The implicit part is the largest |eigenvalue| of the friction tensor:
    // on a concave wall (negative friction) a fully explicit term grows
    // slowly in time (Taylor-Couette with a free-slip outer wall), with
    // +|f| implicit and -2|f| explicit the lagged part is damped.
    zsIndex_.setSize(mesh_.nCells());
    zsIndex_ = -1;
    DynamicList<label> cells;
    DynamicList<symmTensor> friction;
    DynamicList<scalar> implicitCoeff;
    DynamicList<symmTensor> wallMap;
    DynamicList<vector> offset;

    const scalarField& V = mesh_.V();
    const vectorField& C = mesh_.cellCentres();
    const vectorField& offC = centroidOffset_.primitiveField();
    label nFallback = 0;
    label nLimited = 0;
    label nFlat = 0;
    scalar kMin = GREAT;
    scalar kMax = -GREAT;

    // second-difference stencil: centre, +-e_i, and +-e_a +-e_b (a < b)
    const label nSt = 19;

    forAll(bodies_, bodyi)
    {
        const ibmBody& body = bodies_[bodyi];
        if (body.wall() != ibmBody::FREESLIP) continue;

        DynamicList<label> bodyCells;
        for (const label celli : wallCells_)
        {
            if (cellBody_[celli] == bodyi) bodyCells.append(celli);
        }
        if (!body.smoothLevelSet())
        {
            nFallback += bodyCells.size();
            continue;
        }

        pointField pts(nSt*bodyCells.size());
        scalarField eps(bodyCells.size());
        forAll(bodyCells, i)
        {
            const label celli = bodyCells[i];
            const point x(C[celli] + offC[celli]);
            const scalar e = 1e-3*cbrt(V[celli]);
            eps[i] = e;
            label m = nSt*i;
            pts[m++] = x;
            for (direction d = 0; d < vector::nComponents; ++d)
            {
                vector ed(Zero);
                ed[d] = e;
                pts[m++] = x + ed;
                pts[m++] = x - ed;
            }
            for (direction a = 0; a < vector::nComponents; ++a)
            {
                for (direction b = a + 1; b < vector::nComponents; ++b)
                {
                    vector ea(Zero);
                    vector eb(Zero);
                    ea[a] = e;
                    eb[b] = e;
                    pts[m++] = x + ea + eb;
                    pts[m++] = x + ea - eb;
                    pts[m++] = x - ea + eb;
                    pts[m++] = x - ea - eb;
                }
            }
        }
        scalarField phi;
        body.evaluate(pts, phi);

        forAll(bodyCells, i)
        {
            const label celli = bodyCells[i];
            const label m = nSt*i;
            const scalar e = eps[i];
            const scalar phi0 = phi[m];

            vector g;
            tensor hess(Zero);
            for (direction d = 0; d < vector::nComponents; ++d)
            {
                const scalar pp = phi[m + 1 + 2*d];
                const scalar pm = phi[m + 2 + 2*d];
                g[d] = (pp - pm)/(2*e);
                hess(d, d) = (pp + pm - 2*phi0)/sqr(e);
            }
            label mm = m + 7;
            for (direction a = 0; a < vector::nComponents; ++a)
            {
                for (direction b = a + 1; b < vector::nComponents; ++b)
                {
                    const scalar h =
                        (phi[mm] - phi[mm + 1] - phi[mm + 2] + phi[mm + 3])
                       /(4*sqr(e));
                    hess(a, b) = h;
                    hess(b, a) = h;
                    mm += 4;
                }
            }

            const scalar magG = mag(g);
            if (magG < SMALL)
            {
                // no normal (a ridge of the distance function): no friction
                ++nFallback;
                continue;
            }
            const vector n(g/magG);

            // orthonormal tangent basis and the 2x2 curvature matrix in it
            const direction kMinCmpt =
                (mag(n.x()) <= mag(n.y()) && mag(n.x()) <= mag(n.z())) ? 0
              : (mag(n.y()) <= mag(n.z()) ? 1 : 2);
            vector ek(Zero);
            ek[kMinCmpt] = 1;
            vector t1(n ^ ek);
            t1 /= mag(t1);
            const vector t2(n ^ t1);
            const scalar k11 = (t1 & (hess & t1))/magG;
            const scalar k22 = (t2 & (hess & t2))/magG;
            const scalar k12 = 0.5*((t1 & (hess & t2)) + (t2 & (hess & t1)))/magG;

            // principal curvatures and directions (Jacobi rotation)
            const scalar ang = 0.5*atan2(2*k12, k11 - k22);
            const scalar cs = cos(ang);
            const scalar sn = sin(ang);
            const scalar kappa[2] =
            {
                k11*sqr(cs) + 2*k12*sn*cs + k22*sqr(sn),
                k11*sqr(sn) - 2*k12*sn*cs + k22*sqr(cs)
            };
            const vector dir[2] = {cs*t1 + sn*t2, -sn*t1 + cs*t2};

            const scalar d = dWall_[celli];
            const scalar a = Awall_[celli]/V[celli];
            symmTensor F(Zero);
            symmTensor W(Zero);
            scalar cImp = 0;
            bool limited = false;
            bool flat = true;
            for (label j = 0; j < 2; ++j)
            {
                scalar r = 1 - phi0*kappa[j];
                if (r < 0.5)
                {
                    r = 0.5;
                    limited = true;
                }
                const scalar kw = kappa[j]/r;
                scalar s = 1 + d*kw;
                if (s < 0.5)
                {
                    s = 0.5;
                    limited = true;
                }
                F += (a*kw/s)*sqr(dir[j]);
                W += (1/s)*sqr(dir[j]);
                cImp = max(cImp, mag(a*kw/s));
                kMin = min(kMin, kw);
                kMax = max(kMax, kw);
                if (mag(kw)*cbrt(V[celli]) > 1e-8) flat = false;
            }
            if (limited) ++nLimited;
            if (flat) ++nFlat;

            // tangential offset of the wall-piece centroid from the value
            // point (second order: the wall pieces about the value point)
            vector delta(Zero);
            if (secondOrder_ && zeroShearTangential_)
            {
                const vector& Swc = Sw_[celli];
                const vector rw
                (
                    (wallMoment_[celli] & Swc)/max(magSqr(Swc), VSMALL)
                );
                delta = rw - (rw & n)*n;
            }

            zsIndex_[celli] = cells.size();
            cells.append(celli);
            friction.append(F);
            implicitCoeff.append(cImp);
            wallMap.append(W);
            offset.append(delta);
        }
    }

    zsCells_.transfer(cells);
    zsFriction_.transfer(friction);
    zsImplicit_.transfer(implicitCoeff);
    zsWallMap_.transfer(wallMap);
    zsOffset_.transfer(offset);

    if (!zsReported_)
    {
        zsReported_ = true;
        const label nTotal = returnReduce(zsCells_.size(), sumOp<label>());
        reduce(nFallback, sumOp<label>());
        reduce(nLimited, sumOp<label>());
        reduce(nFlat, sumOp<label>());
        reduce(kMin, minOp<scalar>());
        reduce(kMax, maxOp<scalar>());
        Info<< "    zero-shear free slip: " << nTotal << " wall cells";
        if (nTotal > 0)
        {
            Info<< ", wall curvature " << kMin << " .. " << kMax
                << " (" << nFlat << " flat), " << nLimited
                << " limited (radius of curvature under twice the distance)";
        }
        Info<< endl;
        if (nFallback > 0)
        {
            WarningInFunction
                << nFallback << " free-slip wall cells of bodies without a "
                << "smooth level set (triSurfaceMesh, alpha field) keep "
                << "du/dn = 0: their curvature is not computed" << endl;
        }
    }
}


void Foam::cutCellGeometry::zeroShearRelative
(
    const volVectorField& U,
    vectorField& wP,
    vectorField& ubP,
    scalarField& dunP
) const
{
    wP.setSize(zsCells_.size());
    ubP.setSize(zsCells_.size());
    dunP.setSize(zsCells_.size());
    ubP = Zero;
    dunP = Zero;
    forAll(zsCells_, k)
    {
        wP[k] = U[zsCells_[k]];
    }

    // Second order, from the least-squares gradient over the fluid (a
    // linked group: one gradient at its centroid), G_ij = d_i u_j:
    // - the normal derivative of the normal velocity from continuity,
    //   du_n/dn = -tr(P G P) = n.G.n - tr(G) (n into the fluid). Equal to
    //   n.G.n for a solenoidal field, but n.G.n depends on the cell's own
    //   normal velocity through a negative self-coefficient, which made
    //   the explicit flux unstable (SIMPLE without linking, PIMPLE with
    //   it); the tangential divergence hardly does.
    // - the velocity at the foot of the wall-piece centroid, u_P + delta.G
    if (secondOrder_)
    {
        PtrList<volVectorField> g(vector::nComponents);
        for (direction c = 0; c < vector::nComponents; ++c)
        {
            g.set(c, lsqGrad(U.component(c)(), true).ptr());
        }
        forAll(zsCells_, k)
        {
            const label celli = zsCells_[k];
            const vector nf(-Sw_[celli]/max(Awall_[celli], VSMALL));
            const vector& delta = zsOffset_[k];
            for (direction c = 0; c < vector::nComponents; ++c)
            {
                dunP[k] += nf[c]*(nf & g[c][celli]) - g[c][celli][c];
                if (zeroShearTangential_)
                {
                    wP[k][c] += delta & g[c][celli];
                }
            }
        }
    }

    if (!wallVelocity_) return;

    // body velocity at the value points, the point u_P belongs to
    const vectorField& off =
        (secondOrder_ ? valueOffset_ : centroidOffset_).primitiveField();
    pointField x(zsCells_.size());
    forAll(zsCells_, k)
    {
        x[k] = mesh_.cellCentres()[zsCells_[k]] + off[zsCells_[k]];
    }
    forAll(bodies_, bodyi)
    {
        if (!bodies_[bodyi].hasWallVelocity()) continue;
        vectorField ub;
        bodies_[bodyi].velocity(x, ub);
        forAll(zsCells_, k)
        {
            if (cellBody_[zsCells_[k]] == bodyi)
            {
                // rigid velocity at the same point as wP
                ubP[k] = ub[k] + (bodies_[bodyi].omega() ^ zsOffset_[k]);
            }
        }
    }
    wP -= ubP;
}


Foam::tmp<Foam::fvVectorMatrix> Foam::cutCellGeometry::zeroShearSource
(
    const volScalarField& nu,
    const volVectorField& U
) const
{
    volScalarField::Internal coeff
    (
        IOobject("ibmZeroShearCoeff", mesh_.time().timeName(), mesh_,
                 IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh_,
        dimensionedScalar(nu.dimensions()/dimArea, Zero)
    );
    volVectorField::Internal su
    (
        IOobject("ibmZeroShearSource", mesh_.time().timeName(), mesh_,
                 IOobject::NO_READ, IOobject::NO_WRITE, false),
        mesh_,
        dimensionedVector(nu.dimensions()*U.dimensions()/dimArea, Zero)
    );

    vectorField wP;
    vectorField ubP;
    scalarField dunP;
    zeroShearRelative(U, wP, ubP, dunP);

    forAll(zsCells_, k)
    {
        const label celli = zsCells_[k];
        const ibmBody& body = bodies_[cellBody_[celli]];
        const scalar nuc = nu[celli];
        // tangential: the curvature friction; normal: d u_n/dn (explicit),
        // the normal viscous flux, which without it is carried by the
        // cut-cell pressure as an O(1) layer that does not converge and
        // leaves the velocity first order (doc/freeSlipPlan.md)
        vector flux(zsFriction_[k] & wP[k]);
        flux -= dunP[k]*Sw_[celli]/mesh_.V()[celli];
        if (body.hasWallVelocity())
        {
            // omega x n, n = -Sw/Awall into the fluid, times Awall/V
            flux -= (body.omega() ^ Sw_[celli])/mesh_.V()[celli];
        }
        coeff[celli] = nuc*zsImplicit_[k];
        su[celli] = nuc*(flux - zsImplicit_[k]*U[celli]);
    }

    tmp<fvVectorMatrix> tm(fvm::Sp(coeff, U));
    tm.ref() += su;
    return tm;
}


void Foam::cutCellGeometry::calcLsq() const
{
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const label nCells = mesh_.nCells();
    const label nInternal = mesh_.nInternalFaces();
    const fvBoundaryMesh& bm = mesh_.boundary();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    // Normal matrices A = sum_f w d d over the open faces to fluid
    // neighbours, w = 1/|d|^2; a linked group has one value at its
    // centroid and gets one gradient, fitted over the outer faces of the
    // whole group (the sums of its cells)
    symmTensorField A(nCells, Zero);
    labelList nUsed(nCells, 0);
    auto accumulate = [&](const label celli, const vector& d)
    {
        const scalar w = 1.0/max(magSqr(d), VSMALL);
        A[celli] += w*symm(d*d);
        ++nUsed[celli];
        return vector(w*d);
    };

    DynamicList<label> faces(nInternal);
    DynamicList<vector> wd(nInternal);
    forAll(own, facei)
    {
        const label o = own[facei], n = nei[facei];
        if
        (
            theta_[facei] <= 0 || solidCells_.test(o) || solidCells_.test(n)
         || insideGroup(facei)      // one shared value, no information
        )
        {
            continue;
        }
        const vector& d = dCentroid_[facei];
        faces.append(facei);
        wd.append(accumulate(o, d));
        accumulate(n, d);
    }
    lsqFaces_.transfer(faces);
    lsqWd_.transfer(wd);

    lsqPatchFaces_.setSize(bm.size());
    lsqPatchWd_.setSize(bm.size());
    forAll(bm, patchi)
    {
        lsqPatchFaces_[patchi].clear();
        lsqPatchWd_[patchi].clear();
        const fvPatch& fvp = bm[patchi];
        if (!fvp.coupled()) continue;
        const labelUList& fc = fvp.faceCells();
        const fvsPatchScalarField& tp = theta_.boundaryField()[patchi];
        const fvsPatchVectorField& dp = dCentroid_.boundaryField()[patchi];
        DynamicList<label> pf(fvp.size());
        DynamicList<vector> pw(fvp.size());
        forAll(fvp, i)
        {
            if
            (
                tp[i] <= 0 || solidCells_.test(fc[i])
             || insideGroup(fvp.start() + i)
            )
            {
                continue;
            }
            pf.append(i);
            pw.append(accumulate(fc[i], dp[i]));
        }
        lsqPatchFaces_[patchi].transfer(pf);
        lsqPatchWd_[patchi].transfer(pw);
    }
    groupSum(A);
    groupSum(nUsed);

    // Regularise the empty directions, then invert
    const Vector<label>& solD = mesh_.solutionD();
    label nD = 0;
    for (direction d = 0; d < 3; ++d) if (solD[d] > 0) ++nD;

    lsqInvA_.setSize(nCells);
    lsqInvA_ = Zero;
    lsqOk_.reset();
    lsqOk_.resize(nCells);
    forAll(A, celli)
    {
        if (nUsed[celli] < nD || solidCells_.test(celli)) continue;
        symmTensor Ai(A[celli]);
        for (direction d = 0; d < 3; ++d)
        {
            if (solD[d] < 0)
            {
                Ai[symmTensor::XX + 3*d - (d*(d - 1))/2] += 1.0;   // diagonal: XX, YY, ZZ
            }
        }
        if (det(Ai) > SMALL*pow3(tr(Ai)/3))
        {
            lsqInvA_[celli] = inv(Ai);
            lsqOk_.set(celli);
        }
    }

    // Fluid cells with too few fluid neighbours for a least-squares fit
    // (slivers, mostly in 3-D) take the mean gradient of their fluid face
    // neighbours that have one
    boolList ok(nCells, false);
    for (const label celli : lsqOk_) ok[celli] = true;
    boolList okN;
    syncTools::swapBoundaryCellList(mesh_, ok, okN);

    labelList slot(nCells, -1);
    DynamicList<label> slotCells;
    DynamicList<label> slotCount;
    auto slotOf = [&](const label celli)
    {
        if (slot[celli] < 0)
        {
            slot[celli] = slotCells.size();
            slotCells.append(celli);
            slotCount.append(0);
        }
        ++slotCount[slot[celli]];
        return slot[celli];
    };
    DynamicList<labelPair> fallback;
    forAll(own, facei)
    {
        const label o = own[facei], n = nei[facei];
        if (theta_[facei] <= 0) continue;
        if (ok[n] && !ok[o] && !solidCells_.test(o))
        {
            fallback.append(labelPair(slotOf(o), n));
        }
        if (ok[o] && !ok[n] && !solidCells_.test(n))
        {
            fallback.append(labelPair(slotOf(n), o));
        }
    }
    lsqFallback_.transfer(fallback);

    // Coupled faces whose cell on this side is a fallback cell: the cell on
    // the other side contributes, so its gradient is needed there
    boolList fallbackFace(mesh_.nBoundaryFaces(), false);
    lsqPatchFallback_.setSize(bm.size());
    forAll(bm, patchi)
    {
        lsqPatchFallback_[patchi].clear();
        const fvPatch& fvp = bm[patchi];
        if (!fvp.coupled()) continue;
        const labelUList& fc = fvp.faceCells();
        const fvsPatchScalarField& tp = theta_.boundaryField()[patchi];
        DynamicList<labelPair> pf;
        forAll(fvp, i)
        {
            const label bfi = fvp.start() + i - nInternal;
            if (tp[i] > 0 && okN[bfi] && !ok[fc[i]] && !solidCells_.test(fc[i]))
            {
                pf.append(labelPair(slotOf(fc[i]), i));
                fallbackFace[bfi] = true;
            }
        }
        lsqPatchFallback_[patchi].transfer(pf);
    }
    lsqFallbackCells_.transfer(slotCells);
    lsqFallbackCount_.transfer(slotCount);


    // Cells whose gradient the second-order corrections read: the cells with
    // a value offset (cut cells, linked groups), the wall cells, the cells
    // of the faces with a correction (a face-value offset, a partly wet
    // face or a face-gradient correction) and the fallback cells with the
    // neighbours they take the gradient from
    bitSet near(nCells);
    const vectorField& off = valueOffset_.primitiveField();
    forAll(off, celli)
    {
        if (solidCells_.test(celli)) continue;
        if
        (
            magSqr(off[celli]) > 0
         || (linkGroup_.size() == nCells && linkGroup_[celli] >= 0)
        )
        {
            near.set(celli);
        }
    }
    for (const label celli : wallCells_) near.set(celli);
    for (const label celli : lsqFallbackCells_) near.set(celli);
    for (const labelPair& p : lsqFallback_) near.set(p.second());

    auto needsSnGrad = [](const vector& n, const vector& d, const scalar dc)
    {
        return mag(n - d*dc) > snGradTol;
    };

    DynamicList<label> sgFaces;
    forAll(own, facei)
    {
        const label o = own[facei], n = nei[facei];
        if (theta_[facei] <= 0 || solidCells_.test(o) || solidCells_.test(n))
        {
            continue;
        }
        const bool sg = needsSnGrad
        (
            Sf[facei]/magSf[facei], dCentroid_[facei], cutDeltaCoeffs_[facei]
        );
        if (sg) sgFaces.append(facei);
        if (sg || theta_[facei] < 1 || magSqr(faceValueOffset_[facei]) > 0)
        {
            near.set(o);
            near.set(n);
        }
    }
    snGradFaces_.transfer(sgFaces);

    // Coupled faces: a correction on either side (the face-gradient
    // correction must be the same on both sides, conservation), and the
    // fallback neighbours on the other side
    boolList sgB(mesh_.nBoundaryFaces(), false);
    boolList nearB(mesh_.nBoundaryFaces(), false);
    forAll(bm, patchi)
    {
        const fvPatch& fvp = bm[patchi];
        const labelUList& fc = fvp.faceCells();
        const fvsPatchScalarField& tp = theta_.boundaryField()[patchi];
        forAll(fvp, i)
        {
            if (tp[i] <= 0 || solidCells_.test(fc[i])) continue;
            if (!fvp.coupled())
            {
                if (tp[i] < 1) near.set(fc[i]);   // fixed value on a cut face
                continue;
            }
            const label bfi = fvp.start() + i - nInternal;
            sgB[bfi] = needsSnGrad
            (
                Sf.boundaryField()[patchi][i]/magSf.boundaryField()[patchi][i],
                dCentroid_.boundaryField()[patchi][i],
                cutDeltaCoeffs_.boundaryField()[patchi][i]
            );
            nearB[bfi] =
                sgB[bfi] || tp[i] < 1
             || magSqr(faceValueOffset_.boundaryField()[patchi][i]) > 0
             || fallbackFace[bfi];
        }
    }
    boolList sgN(sgB);
    boolList nearN(nearB);
    syncTools::swapBoundaryFaceList(mesh_, sgN);
    syncTools::swapBoundaryFaceList(mesh_, nearN);

    snGradPatchFaces_.setSize(bm.size());
    forAll(bm, patchi)
    {
        snGradPatchFaces_[patchi].clear();
        const fvPatch& fvp = bm[patchi];
        if (!fvp.coupled()) continue;
        const labelUList& fc = fvp.faceCells();
        const fvsPatchScalarField& tp = theta_.boundaryField()[patchi];
        DynamicList<label> pf;
        forAll(fvp, i)
        {
            if (tp[i] <= 0 || solidCells_.test(fc[i])) continue;
            const label bfi = fvp.start() + i - nInternal;
            if (sgB[bfi] || sgN[bfi]) pf.append(i);
            if (sgB[bfi] || sgN[bfi] || nearB[bfi] || nearN[bfi])
            {
                near.set(fc[i]);
            }
        }
        snGradPatchFaces_[patchi].transfer(pf);
    }

    lsqNearCells_ = near.sortedToc();
    DynamicList<label> nf;
    forAll(lsqFaces_, k)
    {
        const label facei = lsqFaces_[k];
        if (near.test(own[facei]) || near.test(nei[facei])) nf.append(k);
    }
    lsqNearFaces_.transfer(nf);
    lsqNearPatchFaces_.setSize(bm.size());
    forAll(bm, patchi)
    {
        lsqNearPatchFaces_[patchi].clear();
        const labelUList& fc = bm[patchi].faceCells();
        const labelList& pf = lsqPatchFaces_[patchi];
        DynamicList<label> npf;
        forAll(pf, k)
        {
            if (near.test(fc[pf[k]])) npf.append(k);
        }
        lsqNearPatchFaces_[patchi].transfer(npf);
    }

    if (debug)
    {
        Info<< "    least-squares gradient: "
            << returnReduce(lsqNearCells_.size(), sumOp<label>())
            << " cells read by the corrections, "
            << returnReduce(snGradFaces_.size(), sumOp<label>())
            << " internal faces with a face-gradient correction" << endl;
    }

    lsqDone_ = true;
}


Foam::tmp<Foam::volVectorField> Foam::cutCellGeometry::lsqGrad
(
    const volScalarField& f
) const
{
    return lsqGrad(f, false);
}


Foam::tmp<Foam::volVectorField> Foam::cutCellGeometry::lsqGrad
(
    const volScalarField& f,
    const bool near
) const
{
    if (!lsqDone_)
    {
        calcLsq();
    }

    tmp<volVectorField> tg
    (
        volVectorField::New
        (
            "lsqGrad(" + f.name() + ')',
            mesh_,
            f.dimensions()/dimLength,
            extrapolatedCalculatedFvPatchField<vector>::typeName
        )
    );
    volVectorField& g = tg.ref();
    vectorField& gI = g.primitiveFieldRef();
    gI = Zero;

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const scalarField& fI = f.primitiveField();

    // r = sum_f w d (f_N - f_P), in face order: the near cells get the sums
    // of all their faces, as on the whole mesh
    vectorField r(mesh_.nCells(), Zero);
    auto addFace = [&](const label k)
    {
        const label facei = lsqFaces_[k];
        const label o = own[facei], n = nei[facei];
        const scalar df = fI[n] - fI[o];
        r[o] += lsqWd_[k]*df;
        r[n] += lsqWd_[k]*df;       // (-d)(-df) = d df
    };
    if (near)
    {
        for (const label k : lsqNearFaces_) addFace(k);
    }
    else
    {
        forAll(lsqFaces_, k) addFace(k);
    }
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];
        if (!fvp.coupled()) continue;
        const labelUList& fc = fvp.faceCells();
        const scalarField fN(f.boundaryField()[patchi].patchNeighbourField());
        const labelList& pf = lsqPatchFaces_[patchi];
        const vectorField& pw = lsqPatchWd_[patchi];
        auto addPatchFace = [&](const label k)
        {
            const label i = pf[k];
            r[fc[i]] += pw[k]*(fN[i] - fI[fc[i]]);
        };
        if (near)
        {
            for (const label k : lsqNearPatchFaces_[patchi]) addPatchFace(k);
        }
        else
        {
            forAll(pf, k) addPatchFace(k);
        }
    }
    groupSum(r);

    if (near)
    {
        for (const label celli : lsqNearCells_)
        {
            if (lsqOk_.test(celli)) gI[celli] = lsqInvA_[celli] & r[celli];
        }
    }
    else
    {
        for (const label celli : lsqOk_)
        {
            gI[celli] = lsqInvA_[celli] & r[celli];
        }
    }
    g.correctBoundaryConditions();

    // Fallback cells: the mean gradient of their fluid face neighbours with
    // a fit; still exact for a linear field
    {
        vectorField sumG(lsqFallbackCells_.size(), Zero);
        for (const labelPair& p : lsqFallback_)
        {
            sumG[p.first()] += gI[p.second()];
        }
        forAll(mesh_.boundary(), patchi)
        {
            const labelPairList& pf = lsqPatchFallback_[patchi];
            if (pf.empty()) continue;
            const vectorField gN(g.boundaryField()[patchi].patchNeighbourField());
            for (const labelPair& p : pf)
            {
                sumG[p.first()] += gN[p.second()];
            }
        }
        forAll(lsqFallbackCells_, s)
        {
            gI[lsqFallbackCells_[s]] = sumG[s]/lsqFallbackCount_[s];
        }
    }

    g.correctBoundaryConditions();

    // Check (DebugSwitches cutCellGeometry 2): the near cells as on the whole mesh
    if (near && debug > 1)
    {
        const tmp<volVectorField> tfull(lsqGrad(f, false));
        label nDiff = 0;
        for (const label celli : lsqNearCells_)
        {
            if (tfull()[celli] != gI[celli]) ++nDiff;
        }
        reduce(nDiff, sumOp<label>());
        Info<< "    lsqGrad(" << f.name() << ") near cells "
            << returnReduce(lsqNearCells_.size(), sumOp<label>())
            << ", differing from the whole-mesh gradient " << nDiff << endl;
        if (nDiff)
        {
            FatalErrorInFunction
                << "near-cell gradient differs on " << nDiff << " cells"
                << exit(FatalError);
        }
    }

    return tg;
}


void Foam::cutCellGeometry::report() const
{
    const label nSolid = returnReduce(solidCells_.count(), sumOp<label>());
    const label nCut = returnReduce(cutCells_.count(), sumOp<label>());
    const label nWall = returnReduce(wallCells_.count(), sumOp<label>());
    const label nCells = mesh_.globalData().nTotalCells();

    vector SwSum(gSum(Sw_.primitiveField()));
    scalar AwallSum(gSum(Awall_.primitiveField()));
    scalar minAlphaCut = GREAT;
    for (const label celli : cutCells_)
    {
        minAlphaCut = min(minAlphaCut, alpha_[celli]);
    }
    reduce(minAlphaCut, minOp<scalar>());

    Info<< "    cells: " << nCells << " total, " << nSolid << " solid, "
        << nCut << " cut, " << nWall << " with a wall segment" << nl
        << "    fluid volume " << fluidVolume_
        << " (domain " << gSum(mesh_.V()) << ")" << nl
        << "    total wall area " << AwallSum
        << ", net wall area vector " << SwSum
        << " (zero for closed bodies inside the domain)" << nl
        << "    smallest retained cell fraction "
        << (nCut > 0 ? minAlphaCut : 1.0) << nl << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::cutCellGeometry::grad
(
    const volScalarField& p
) const
{
    tmp<volVectorField> tg
    (
        volVectorField::New
        (
            "cutCellGrad(" + p.name() + ')',
            mesh_,
            p.dimensions()/dimLength,
            extrapolatedCalculatedFvPatchField<vector>::typeName
        )
    );
    volVectorField& g = tg.ref();
    vectorField& gI = g.primitiveFieldRef();
    gI = Zero;

    const surfaceScalarField pf(fvc::interpolate(p));
    const surfaceVectorField& Sf = mesh_.Sf();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    forAll(own, facei)
    {
        const vector v(theta_[facei]*pf[facei]*Sf[facei]);
        gI[own[facei]] += v;
        gI[nei[facei]] -= v;
    }
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];
        const labelUList& fc = fvp.faceCells();
        const fvsPatchScalarField& tp = theta_.boundaryField()[patchi];
        const fvsPatchScalarField& pfp = pf.boundaryField()[patchi];
        const fvsPatchVectorField& Sfp = Sf.boundaryField()[patchi];
        forAll(fvp, i)
        {
            gI[fc[i]] += tp[i]*pfp[i]*Sfp[i];
        }
    }

    // Wall pressure force, p_P Sw, through closure; the sum is then
    // sum_f theta_f (p_f - p_P) Sf and vanishes for a uniform pressure
    gI += p.primitiveField()*Sw_.primitiveField();

    if (secondOrder_)
    {
        // Face values at the wet-face centroids, g_f.(x'_f - m_f) added to
        // the linear interpolate (the value at m_f), and the wall pressure
        // integrated over the wall pieces, sum_k [p_P + g_P.(x_k - x_c)] S_k:
        // with the least-squares gradient g exact for a linear pressure, the
        // sum is the exact integral over the geometric fluid polyhedron
        const tmp<volVectorField> tgp(lsqGrad(p, true));
        const volVectorField& gp = tgp();
        const surfaceScalarField& weights = mesh_.weights();
        forAll(own, facei)
        {
            const vector& r = faceValueOffset_[facei];
            if (magSqr(r) == 0) continue;
            const scalar w = weights[facei];
            const vector gf(w*gp[own[facei]] + (1 - w)*gp[nei[facei]]);
            const vector v(theta_[facei]*(gf & r)*Sf[facei]);
            gI[own[facei]] += v;
            gI[nei[facei]] -= v;
        }
        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& fvp = mesh_.boundary()[patchi];
            if (!fvp.coupled()) continue;
            const labelUList& fc = fvp.faceCells();
            const vectorField gN(gp.boundaryField()[patchi].patchNeighbourField());
            const fvsPatchScalarField& tp = theta_.boundaryField()[patchi];
            const fvsPatchVectorField& rp = faceValueOffset_.boundaryField()[patchi];
            const fvsPatchScalarField& wp = weights.boundaryField()[patchi];
            const fvsPatchVectorField& Sfp = Sf.boundaryField()[patchi];
            forAll(fvp, i)
            {
                if (magSqr(rp[i]) == 0) continue;
                const vector gf(wp[i]*gp[fc[i]] + (1 - wp[i])*gN[i]);
                gI[fc[i]] += tp[i]*(gf & rp[i])*Sfp[i];
            }
        }
        for (const label celli : wallCells_)
        {
            gI[celli] += gp[celli] & wallMoment_[celli];
        }
    }

    gI /= mesh_.V();

    g.correctBoundaryConditions();
    return tg;
}


Foam::tmp<Foam::surfaceScalarField> Foam::cutCellGeometry::flux
(
    const volVectorField& U
) const
{
    tmp<surfaceScalarField> tphi(theta_*fvc::flux(U));
    tphi.ref().rename("flux(" + U.name() + ')');

    // Second order: the flux at the wet-face centroid, g_f.(x'_f - m_f)
    // added to the linear interpolate (the value at m_f) as for the other
    // face values. Without it the continuity flux of a cut face is first
    // order, which costs the cut cells and the wall gradient their accuracy
    // and, through the group centroid's distance to its faces, makes linked
    // cells first order.
    if (secondOrder_)
    {
        surfaceScalarField& phi = tphi.ref();
        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        const surfaceScalarField& weights = mesh_.weights();
        const surfaceVectorField& Sf = mesh_.Sf();
        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            if (mesh_.solutionD()[cmpt] < 0) continue;
            const tmp<volVectorField> tg(lsqGrad(U.component(cmpt), true));
            const volVectorField& g = tg();
            forAll(own, facei)
            {
                const vector& r = faceValueOffset_[facei];
                if (magSqr(r) == 0 || theta_[facei] <= 0) continue;
                const scalar w = weights[facei];
                const vector gf(w*g[own[facei]] + (1 - w)*g[nei[facei]]);
                phi[facei] += theta_[facei]*(gf & r)*Sf[facei][cmpt];
            }
            forAll(mesh_.boundary(), patchi)
            {
                const fvPatch& fvp = mesh_.boundary()[patchi];
                if (!fvp.coupled()) continue;
                const labelUList& fc = fvp.faceCells();
                const vectorField gN(g.boundaryField()[patchi].patchNeighbourField());
                const fvsPatchScalarField& tp = theta_.boundaryField()[patchi];
                const fvsPatchVectorField& rp = faceValueOffset_.boundaryField()[patchi];
                const fvsPatchScalarField& wp = weights.boundaryField()[patchi];
                const fvsPatchVectorField& Sfp = Sf.boundaryField()[patchi];
                fvsPatchScalarField& php = phi.boundaryFieldRef()[patchi];
                forAll(fvp, i)
                {
                    if (magSqr(rp[i]) == 0 || tp[i] <= 0) continue;
                    const vector gf(wp[i]*g[fc[i]] + (1 - wp[i])*gN[i]);
                    php[i] += tp[i]*(gf & rp[i])*Sfp[i][cmpt];
                }
            }
        }
    }
    return tphi;
}


Foam::tmp<Foam::surfaceScalarField> Foam::cutCellGeometry::thetaInterpolate
(
    const volScalarField& f
) const
{
    tmp<surfaceScalarField> tf(theta_*fvc::interpolate(f));
    tf.ref().rename("theta*interpolate(" + f.name() + ')');
    return tf;
}


Foam::tmp<Foam::surfaceScalarField> Foam::cutCellGeometry::thetaInterpolate
(
    const tmp<volScalarField>& tf
) const
{
    tmp<surfaceScalarField> tr(thetaInterpolate(tf()));
    tf.clear();
    return tr;
}


Foam::tmp<Foam::surfaceScalarField> Foam::cutCellGeometry::ddtCorr
(
    const volScalarField& rAU,
    const volVectorField& U,
    const surfaceScalarField& phi
) const
{
    // Euler form of fvc::ddtCorr with the theta-weighted flux of U.oldTime()
    static bool warned = false;
    const word schemeName
    (
        mesh_.ddtScheme("ddt(" + alpha_.name() + ',' + U.name() + ')')
    );
    if (schemeName != "Euler" && !warned)
    {
        WarningInFunction
            << "ddtCorr is implemented for the Euler scheme; the Euler form "
            << "is used with the " << schemeName << " scheme" << endl;
        warned = true;
    }

    const dimensionedScalar rDeltaT(1.0/mesh_.time().deltaT());
    const surfaceScalarField phiU0(flux(U.oldTime()));
    surfaceScalarField phi0(phi.oldTime());

    if (moving_ && thetaOldPtr_)
    {
        // The old flux lives on the faces of the previous cut mesh: scale it
        // to the new wet area, and give faces that have just opened the
        // interpolated old velocity (no correction there)
        const surfaceScalarField& thetaOld = thetaOldPtr_();
        forAll(phi0, facei)
        {
            phi0[facei] = (thetaOld[facei] > 0)
              ? phi0[facei]*theta_[facei]/thetaOld[facei]
              : phiU0[facei];
        }
        forAll(phi0.boundaryField(), patchi)
        {
            fvsPatchScalarField& p0 = phi0.boundaryFieldRef()[patchi];
            const fvsPatchScalarField& t0 = thetaOld.boundaryField()[patchi];
            const fvsPatchScalarField& t1 = theta_.boundaryField()[patchi];
            const fvsPatchScalarField& pU = phiU0.boundaryField()[patchi];
            forAll(p0, i)
            {
                p0[i] = (t0[i] > 0) ? p0[i]*t1[i]/t0[i] : pU[i];
            }
        }
    }

    // Same coupling coefficient as ddtScheme::fvcDdtPhiCoeff
    const surfaceScalarField coeff
    (
        1.0
      - min
        (
            mag(phi0 - phiU0)
           /(mag(phi0) + dimensionedScalar(phi0.dimensions(), SMALL)),
            1.0
        )
    );

    // HbyA holds rAU alpha0/dt U.oldTime(), alpha0 the coefficient of
    // U.oldTime() in ddt(alpha, U); mergeH gives a group the volume mean of
    // its cells. Without it the gain of the old flux is rAU/dt ~ 1/(1.5 alpha)
    // at small cut cells, above 1: a difference between the two copies of a
    // coupled face flux then grows from round-off.
    volScalarField alpha0("alpha0", alpha_.oldTime());
    groupMean(alpha0);

    tmp<surfaceScalarField> tcorr
    (
        fvc::interpolate(rAU*alpha0)*coeff*rDeltaT*(phi0 - phiU0)
    );
    tcorr.ref().rename("ddtCorr(" + U.name() + ',' + phi.name() + ')');
    return tcorr;
}


void Foam::cutCellGeometry::syncCoupledFlux(surfaceScalarField& phi) const
{
    const label nInternal = mesh_.nInternalFaces();
    scalarField nbr(mesh_.nBoundaryFaces(), Zero);
    forAll(mesh_.boundary(), patchi)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        if (!isA<processorPolyPatch>(pp) && !isA<cyclicPolyPatch>(pp)) continue;
        SubList<scalar>(nbr, pp.size(), pp.start() - nInternal) =
            phi.boundaryField()[patchi];
    }
    syncTools::swapBoundaryFaceList(mesh_, nbr);
    forAll(mesh_.boundary(), patchi)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        if (!isA<processorPolyPatch>(pp) && !isA<cyclicPolyPatch>(pp)) continue;
        fvsPatchScalarField& php = phi.boundaryFieldRef()[patchi];
        forAll(php, i)
        {
            // The copy across holds the flux out of its own cell: -php[i]
            php[i] = 0.5*(php[i] - nbr[pp.start() - nInternal + i]);
        }
    }
}


void Foam::cutCellGeometry::constrainPressure(fvScalarMatrix& pEqn) const
{
    // Every face of a solid cell is closed, so its row is empty. Give it a
    // diagonal of the typical magnitude of the fluid rows (any non-zero
    // value gives p = 0 there; the magnitude only matters for the linear
    // solver, and the normalisation of its residual: from the diagonal with
    // the coefficients of the coupled faces, the same in every
    // decomposition) and a zero source.
    const scalarField diag0(coupledDiag(pEqn));
    scalarField& diag = pEqn.diag();
    scalarField& source = pEqn.source();

    scalar sumMag = 0;
    label n = 0;
    forAll(diag0, celli)
    {
        if (!solidCells_.test(celli))
        {
            sumMag += mag(diag0[celli]);
            ++n;
        }
    }
    reduce(sumMag, sumOp<scalar>());
    reduce(n, sumOp<label>());
    const scalar scale = (n > 0) ? sumMag/n : 1.0;

    for (const label celli : solidCells_)
    {
        diag[celli] = -scale;
        source[celli] = 0;
    }
}


void Foam::cutCellGeometry::checkRefCell(const label pRefCell) const
{
    bool bad = (pRefCell >= 0 && solidCells_.test(pRefCell));
    reduce(bad, orOp<bool>());
    if (bad)
    {
        FatalErrorInFunction
            << "The pressure reference cell lies inside an immersed body. "
            << "Choose pRefPoint/pRefCell in the fluid." << exit(FatalError);
    }
}


void Foam::cutCellGeometry::forces
(
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& nu,
    List<vector>& Fp,
    List<vector>& Fv,
    const volScalarField* nuMolecularPtr
) const
{
    List<vector> Mp, Mv;
    forces(p, U, nu, Fp, Fv, Mp, Mv, nuMolecularPtr);
}


void Foam::cutCellGeometry::forces
(
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& nu,
    List<vector>& Fp,
    List<vector>& Fv,
    List<vector>& Mp,
    List<vector>& Mv,
    const volScalarField* nuMolecularPtr
) const
{
    Fp.setSize(bodies_.size());
    Fv.setSize(bodies_.size());
    Mp.setSize(bodies_.size());
    Mv.setSize(bodies_.size());
    Fp = Zero;
    Fv = Zero;
    Mp = Zero;
    Mv = Zero;
    const scalarField& V = mesh_.V();
    const vectorField& C = mesh_.cellCentres();
    const vectorField& valueOff = valueOffset_.primitiveField();
    const vectorField& centroidOff = centroidOffset_.primitiveField();

    tmp<volVectorField> tgp;
    if (secondOrder_) tgp = lsqGrad(p, true);

    // zero-shear free slip: relative velocity at the value points
    vectorField zsW;
    vectorField zsUb;
    scalarField zsDun;
    if (zeroShear_) zeroShearRelative(U, zsW, zsUb, zsDun);

    for (const label celli : wallCells_)
    {
        const label bodyi = max(cellBody_[celli], 0);
        if (bodyi >= bodies_.size()) continue;
        // pressure on the wall, p_P Sw, plus with the second-order treatment
        // sum_k g_P.(x_k - x_c) S_k over the wall pieces (as in grad())
        vector fp(p[celli]*Sw_[celli]);

        // moment of fp about the base point x_b of the wall pieces (second
        // order: sum_k (x_k - x_b) ^ S_k (p_P + g_P.(x_k - x_b)); first order:
        // p_P Sw at the foot of the wall normal, no moment) and wall centroid
        // r_w relative to x_b, the lever of the viscous traction
        // (T: sum_k (x_k - x_b) S_k, first order the single piece Sw at x_w)
        vector mp(Zero);
        point xb;
        vector rw;
        tensor T;
        if (secondOrder_)
        {
            T = wallMoment_[celli];
            const vector& g = tgp()[celli];
            xb = C[celli] + valueOff[celli];
            fp += g & T;
            mp = p[celli]*vector(T.yz() - T.zy(), T.zx() - T.xz(), T.xy() - T.yx())
               + (g & wallTorqueMoment_[celli]);
            rw = (T & Sw_[celli])/max(magSqr(Sw_[celli]), VSMALL);
        }
        else
        {
            xb = C[celli] + centroidOff[celli];
            rw = dWall_[celli]*Sw_[celli]/max(Awall_[celli], VSMALL);
            T = rw*Sw_[celli];
        }

        // viscous traction, nu (u_P - u_b)/d Awall
        const label zsk = zeroShear_ ? zsIndex_[celli] : -1;
        vector fv;
        if (zsk >= 0)
        {
            // zero-shear free slip: the wall flux of the Laplacian form,
            // nu (H_w.w_t,w + omega x n + (n.grad(u).n) n) Awall (see
            // zeroShearSource)
            fv = nu[celli]*V[celli]*(zsFriction_[zsk] & zsW[zsk])
               - nu[celli]*zsDun[zsk]*Sw_[celli];
            if (bodies_[bodyi].hasWallVelocity())
            {
                fv -= nu[celli]*(bodies_[bodyi].omega() ^ Sw_[celli]);
            }
        }
        else if (wallVelocity_)
        {
            fv = nu[celli]*noSlipCoeff_[celli]*V[celli]
                *(U[celli] - UbPtr_()[celli]);
        }
        else
        {
            fv = nu[celli]*noSlipCoeff_[celli]*V[celli]*U[celli];
        }

        const vector rb(bodies_[bodyi].relativePosition(xb));
        Fp[bodyi] += fp;
        Fv[bodyi] += fv;
        Mp[bodyi] += (rb ^ fp) + mp;
        Mv[bodyi] += (rb + rw) ^ fv;

        // The traction above is nu du/dn, the wall flux of the Laplacian
        // form of the momentum equation; the stress is nu (grad u + grad u^T)
        // with the molecular viscosity at the wall. On the no-slip wall of a
        // rigid body (u = U + omega ^ r) the transpose part is nu (n ^ omega)
        // per unit area, n out of the body (-Sw/|Sw|): nu omega ^ S_k on
        // piece k. It sums to zero over a closed body, its moment does not
        // (-2 pi nu omega R^2 on a rotating circular cylinder):
        // sum_k (x_k - x_b) ^ (omega ^ S_k) = omega tr(T) - omega & T
        const vector& omega = bodies_[bodyi].omega();
        if (zsk >= 0)
        {
            // Zero-shear free slip: the tangential stress vanishes, so the
            // transpose part nu (grad u)^T.n is minus the Laplacian-form
            // flux above, plus a normal part. Its force sums to zero over
            // a closed wall (div u = 0); its moment is
            // nu sum (u_w x n) Awall = nu sum Sw x u_w, with the fluid
            // velocity on the wall u_w = u_b,w + (I + d H_w)^-1 P w_P and
            // u_b,w = u_b(x_P) - d omega x n (independent of the origin)
            vector uw(zsWallMap_[zsk] & zsW[zsk]);
            if (bodies_[bodyi].hasWallVelocity())
            {
                const vector nf(-Sw_[celli]/max(Awall_[celli], VSMALL));
                uw += zsUb[zsk] - dWall_[celli]*(omega ^ nf);
            }
            // (the viscosity of the flux above, not the molecular one: a
            // free-slip wall has no viscous sublayer)
            Mv[bodyi] += nu[celli]*(Sw_[celli] ^ uw);
        }
        else if (magSqr(omega) > 0)
        {
            const scalar nuMol =
                nuMolecularPtr ? (*nuMolecularPtr)[celli] : nu[celli];
            const vector fr(nuMol*(omega ^ Sw_[celli]));
            Fv[bodyi] += fr;
            Mv[bodyi] += (rb ^ fr) + nuMol*(omega*tr(T) - (omega & T));
        }
    }
    forAll(Fp, i)
    {
        reduce(Fp[i], sumOp<vector>());
        reduce(Fv[i], sumOp<vector>());
        reduce(Mp[i], sumOp<vector>());
        reduce(Mv[i], sumOp<vector>());
    }
}


void Foam::cutCellGeometry::reportForces
(
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& nu,
    const volScalarField* nuMolecularPtr
) const
{
    List<vector> Fp, Fv, Mp, Mv;
    forces(p, U, nu, Fp, Fv, Mp, Mv, nuMolecularPtr);

    forAll(bodies_, bodyi)
    {
        Info<< "IBM force on " << bodies_[bodyi].name()
            << ": pressure " << Fp[bodyi]
            << " viscous " << Fv[bodyi]
            << " total " << Fp[bodyi] + Fv[bodyi] << endl;
        Info<< "IBM moment on " << bodies_[bodyi].name()
            << " about " << bodies_[bodyi].currentCentre()
            << ": pressure " << Mp[bodyi]
            << " viscous " << Mv[bodyi]
            << " total " << Mp[bodyi] + Mv[bodyi] << endl;
    }

    if (writeForces_ && Pstream::master())
    {
        writeForceFile(forcesFilePtr_, "forces.dat", "Forces", "F", Fp, Fv);
        writeForceFile
        (
            momentsFilePtr_, "moments.dat",
            "Moments about the body centres (ibmBody::currentCentre)", "M",
            Mp, Mv
        );
    }
}


void Foam::cutCellGeometry::writeForceFile
(
    autoPtr<OFstream>& filePtr,
    const word& fileBaseName,
    const string& title,
    const word& symbol,
    const List<vector>& Xp,
    const List<vector>& Xv
) const
{
    if (!filePtr)
    {
        const fileName dir
        (
            mesh_.time().globalPath()/"postProcessing"/"ibmForces"
           /mesh_.time().timeName()
        );
        mkDir(dir);
        filePtr.reset(new OFstream(dir/fileBaseName));
        OFstream& os = filePtr();
        os  << "# " << title.c_str() << " per unit density on the immersed bodies" << nl
            << "# time";
        forAll(bodies_, bodyi)
        {
            os  << "  " << bodies_[bodyi].name()
                << ':' << symbol << "p(x y z)  " << bodies_[bodyi].name()
                << ':' << symbol << "v(x y z)";
        }
        os  << nl;
    }
    OFstream& os = filePtr();
    os  << mesh_.time().timeName();
    forAll(bodies_, bodyi)
    {
        os  << ' ' << Xp[bodyi].x() << ' ' << Xp[bodyi].y() << ' ' << Xp[bodyi].z()
            << ' ' << Xv[bodyi].x() << ' ' << Xv[bodyi].y() << ' ' << Xv[bodyi].z();
    }
    os  << nl;
    os.flush();
}


void Foam::cutCellGeometry::writeGeometryFields() const
{
    Info<< "    writing geometry fields to time "
        << mesh_.time().timeName() << endl;
    alpha_.write();
    theta_.write();
    Sw_.write();
    dWall_.write();
    noSlipCoeff_.write();
    bodyIndex_.write();

    if (mergeAlpha_ > 0)
    {
        // Linked groups, comparable between decompositions: role (1 slave,
        // 2 master) and the fluid centroid of the master of every group cell
        volScalarField role
        (
            IOobject("ibmLinkRole", mesh_.time().timeName(), mesh_),
            mesh_,
            dimensionedScalar(dimless, Zero)
        );
        volVectorField master
        (
            IOobject("ibmLinkMaster", mesh_.time().timeName(), mesh_),
            mesh_,
            dimensionedVector(dimLength, Zero)
        );
        const pointField xc(fluidCentroids());
        forAll(linkGroup_, celli)
        {
            if (linkGroup_[celli] < 0) continue;
            role[celli] = (linkFace_[celli] >= 0) ? 1 : 2;
            master[celli] = xc[celli] - linkRel_[celli];
        }
        role.write();
        master.write();
    }

    // Wall-law samples, comparable between decompositions: the distance of
    // the sample from the wall (0: the wall cell samples itself). Only for
    // level-set bodies (the samples need the wall distance).
    bool wallFunctions = false;
    forAll(bodies_, bodyi)
    {
        if (bodies_[bodyi].wallFunction() != ibmBody::NOWALLFUNCTION)
        {
            wallFunctions = true;
        }
        if (!bodies_[bodyi].hasLevelSet())
        {
            wallFunctions = false;
            break;
        }
    }
    if (wallFunctions && wallSampleDistance_ > 0)
    {
        if (!wallSamplesDone_)
        {
            calcWallSamples();
        }
        volScalarField sample
        (
            IOobject("ibmWallSampleDist", mesh_.time().timeName(), mesh_),
            mesh_,
            dimensionedScalar(dimLength, Zero)
        );
        sample.primitiveFieldRef() = sampleDist_;
        sample.write();
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::cutCellGeometry::convectionFlux
(
    const surfaceScalarField& phi
) const
{
    tmp<surfaceScalarField> tf
    (
        new surfaceScalarField
        (
            IOobject
            (
                "ibmConvectionFlux(" + phi.name() + ')',
                phi.instance(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(phi.dimensions(), Zero)
        )
    );
    if (!ibmConvection_)
    {
        return tf;
    }
    surfaceScalarField& f = tf.ref();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const vectorField& off = valueOffset_.primitiveField();
    forAll(own, facei)
    {
        const label o = own[facei], n = nei[facei];
        const label u = (phi[facei] >= 0) ? o : n;
        if
        (
            theta_[facei] > 0 && !solidCells_.test(o) && !solidCells_.test(n)
         && (theta_[facei] < 1 || magSqr(off[u]) > 0)
        )
        {
            f[facei] = phi[facei];
        }
    }
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];
        if (!fvp.coupled()) continue;
        const labelUList& fc = fvp.faceCells();
        const vectorField offN
        (
            valueOffset_.boundaryField()[patchi].patchNeighbourField()
        );
        const fvsPatchScalarField& tp = theta_.boundaryField()[patchi];
        const fvsPatchScalarField& php = phi.boundaryField()[patchi];
        fvsPatchScalarField& fp = f.boundaryFieldRef()[patchi];
        forAll(fvp, i)
        {
            const vector& offU = (php[i] >= 0) ? off[fc[i]] : offN[i];
            if
            (
                tp[i] > 0 && !solidCells_.test(fc[i])
             && (tp[i] < 1 || magSqr(offU) > 0)
            )
            {
                fp[i] = php[i];
            }
        }
    }
    return tf;
}
