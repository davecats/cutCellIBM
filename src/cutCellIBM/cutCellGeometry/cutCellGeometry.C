#include "cutCellGeometry.H"
#include "ibmBodies.H"
#include "fvMesh.H"
#include "Time.H"
#include "fvc.H"
#include "syncTools.H"
#include "volPointInterpolation.H"
#include "pointFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "coupledFvPatchField.H"
#include "boundBox.H"

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
    alphaMin_(getOrDefault<scalar>("alphaMin", 0.02)),
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
    cellBody_(mesh.nCells(), -1),
    nearestBody_(mesh.nCells(), -1),
    solidCells_(mesh.nCells()),
    cutCells_(mesh.nCells()),
    wallCells_(mesh.nCells()),
    fluidVolume_(0),
    forcesFilePtr_(nullptr),
    moving_(false),
    periodicLengths_(getOrDefault<vector>("periodicLengths", Zero)),
    cellLevelSet_(),
    UbPtr_(nullptr),
    massSourcePtr_(nullptr),
    thetaOldPtr_(nullptr),
    alphaOldMin_(1),
    nSwept_(0)
{
    Info<< "Cut-cell immersed boundary: reading " << name() << nl
        << "    alphaMin " << alphaMin_
        << ", nSubSamples " << nSubSamples_
        << ", nBisect " << nBisect_ << endl;

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
        bodies_[bodyi].setPeriodicLengths(periodicLengths_);
        if (bodies_[bodyi].moving()) moving_ = true;
    }
    reduce(moving_, orOp<bool>());

    if (moving_)
    {
        Info<< "    moving bodies: geometry rebuilt every time step"
            << (cmptMax(periodicLengths_) > 0
                ? ", periodic lengths " + Foam::name(periodicLengths_) : "")
            << endl;
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

    if (moving_)
    {
        calcMovingTerms();
    }

    if (writeFields_)
    {
        writeGeometryFields();
    }
}


void Foam::cutCellGeometry::calcGeometry()
{
    // Fluid fractions: union of all bodies. Level-set bodies are combined
    // exactly through the minimum of their level sets; alpha-field bodies
    // and the level-set group are combined through the minimum of their
    // fractions, exact wherever they do not share a cell.
    scalarField alpha(mesh_.nCells(), 1.0);
    scalarField theta(mesh_.nFaces(), 1.0);
    cellBody_ = -1;
    nearestBody_ = -1;

    auto merge = [&](const scalarField& a, const scalarField& t,
                     const labelList& body)
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
            theta[facei] = min(theta[facei], t[facei]);
        }
    };

    label nLevelSet = 0;
    forAll(bodies_, bodyi)
    {
        if (bodies_[bodyi].hasLevelSet()) ++nLevelSet;
    }

    if (nLevelSet > 0)
    {
        scalarField a, t;
        labelList body;
        calcFromLevelSet(a, t, body);
        merge(a, t, body);
        unionLevelSet(mesh_.cellCentres(), cellLevelSet_);
    }

    forAll(bodies_, bodyi)
    {
        if (!bodies_[bodyi].hasLevelSet())
        {
            scalarField a, t;
            calcFromAlphaField
            (
                refCast<const alphaFieldBody>(bodies_[bodyi]), a, t
            );
            merge(a, t, labelList(mesh_.nCells(), bodyi));
        }
    }

    // Absorb slivers (C1): a deactivated cell also has its faces closed, or
    // the discrete divergence would see a flux into a cell without equation
    const cellList& cells = mesh_.cells();
    label nAbsorbed = 0;
    forAll(alpha, celli)
    {
        if (alpha[celli] < alphaMin_)
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
    {
        scalarField thetaUnsynced(theta);
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
                        maxDiff, mag(theta[facei] - thetaUnsynced[facei])
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

    alpha_.primitiveFieldRef() = alpha;
    alpha_.correctBoundaryConditions();

    theta_.primitiveFieldRef() = SubField<scalar>(theta, mesh_.nInternalFaces());
    forAll(theta_.boundaryField(), patchi)
    {
        const fvPatch& fvp = mesh_.boundary()[patchi];
        theta_.boundaryFieldRef()[patchi] =
            SubField<scalar>(theta, fvp.size(), fvp.start());
    }

    calcDerived();
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


void Foam::cutCellGeometry::faceFractions
(
    const scalarField& pointPhi,
    const pointField& edgeCross,
    scalarField& theta
) const
{
    const faceList& faces = mesh_.faces();
    const pointField& points = mesh_.points();
    const labelListList& faceEdges = mesh_.faceEdges();
    const edgeList& edges = mesh_.edges();
    const vectorField& Sf = mesh_.faceAreas();

    theta.setSize(faces.size());

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
        const point& p0 = poly[0];
        for (label i = 1; i + 1 < poly.size(); ++i)
        {
            area += 0.5*((poly[i] - p0) ^ (poly[i+1] - p0));
        }

        const scalar magSf = mag(Sf[facei]);
        theta[facei] =
            (magSf > VSMALL) ? min(max(mag(area)/magSf, 0.0), 1.0) : 0.0;
    }
}


Foam::scalar Foam::cutCellGeometry::sampleCellFraction(const label celli) const
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
    unionLevelSet(pointField(samples), phi);

    label nFluid = 0;
    forAll(phi, i)
    {
        if (phi[i] > 0) ++nFluid;
    }
    return scalar(nFluid)/samples.size();
}


void Foam::cutCellGeometry::calcFromLevelSet
(
    scalarField& alpha,
    scalarField& theta,
    labelList& cellBody
) const
{
    const pointField& points = mesh_.points();
    const edgeList& edges = mesh_.edges();
    const faceList& faces = mesh_.faces();
    const cellList& cells = mesh_.cells();

    // Level set at the mesh points, face centres and cell centres
    scalarField pointPhi, facePhi, cellPhi;
    unionLevelSet(points, pointPhi);
    unionLevelSet(mesh_.faceCentres(), facePhi);
    unionLevelSet(mesh_.cellCentres(), cellPhi);

    // Crossing point on every cut edge, by bisection of the exact level set.
    // All cut edges are bisected together so that the level set is evaluated
    // in batches (matters for STL surfaces).
    DynamicList<label> cutEdges(edges.size()/8);
    forAll(edges, edgei)
    {
        const edge& e = edges[edgei];
        if ((pointPhi[e[0]] > 0) != (pointPhi[e[1]] > 0))
        {
            cutEdges.append(edgei);
        }
    }

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
        unionLevelSet(mid, phiMid);
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

    pointField edgeCross(edges.size(), point::max);
    forAll(cutEdges, i)
    {
        edgeCross[cutEdges[i]] = A[i] + 0.5*(lo[i] + hi[i])*(B[i] - A[i]);
    }

    // Face fractions
    faceFractions(pointPhi, edgeCross, theta);

    // Cell fractions. A cell is a candidate for cutting if its vertices, face
    // centres or centre do not all lie on the same side; otherwise it is
    // wholly fluid or wholly solid (a body smaller than a cell is not
    // resolved, as in the reference implementation).
    alpha.setSize(cells.size());
    cellBody.setSize(cells.size(), -1);

    label nCandidates = 0;
    forAll(cells, celli)
    {
        const cell& c = cells[celli];
        label nPos = 0, nNeg = 0;
        (cellPhi[celli] > 0 ? nPos : nNeg)++;
        for (const label facei : c)
        {
            (facePhi[facei] > 0 ? nPos : nNeg)++;
            const face& f = faces[facei];
            forAll(f, i)
            {
                (pointPhi[f[i]] > 0 ? nPos : nNeg)++;
            }
        }

        if (nPos > 0 && nNeg > 0)
        {
            alpha[celli] = sampleCellFraction(celli);
            ++nCandidates;
        }
        else
        {
            alpha[celli] = (nPos > 0) ? 1.0 : 0.0;
        }
    }

    // Owning body of every cut or solid cell: the body closest to its centre
    List<scalarField> bodyPhi(bodies_.size());
    forAll(bodies_, bodyi)
    {
        if (bodies_[bodyi].hasLevelSet())
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

    reduce(nCandidates, sumOp<label>());
    Info<< "    level-set bodies: " << returnReduce(cutEdges.size(), sumOp<label>())
        << " cut edges, " << nCandidates << " sub-sampled cells" << endl;
}


void Foam::cutCellGeometry::calcFromAlphaField
(
    const alphaFieldBody& body,
    scalarField& alpha,
    scalarField& theta
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
    const scalarField pointPhi(pAlpha.primitiveField() - 0.5);

    const edgeList& edges = mesh_.edges();
    const pointField& points = mesh_.points();
    pointField edgeCross(edges.size(), point::max);
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

    faceFractions(pointPhi, edgeCross, theta);
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

    tmp<surfaceScalarField> tcorr(coeff*rDeltaT*(phi0 - phiU0));
    tcorr.ref().rename("ddtCorr(" + U.name() + ',' + phi.name() + ')');
    return tcorr;
}


void Foam::cutCellGeometry::constrainPressure(fvScalarMatrix& pEqn) const
{
    // Every face of a solid cell is closed, so its row is empty. Give it a
    // diagonal of the typical magnitude of the fluid rows (any non-zero
    // value gives p = 0 there; the magnitude only matters for the linear
    // solver) and a zero source.
    scalarField& diag = pEqn.diag();
    scalarField& source = pEqn.source();

    scalar sumMag = 0;
    label n = 0;
    forAll(diag, celli)
    {
        if (!solidCells_.test(celli))
        {
            sumMag += mag(diag[celli]);
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
    List<vector>& Fv
) const
{
    Fp.setSize(bodies_.size());
    Fv.setSize(bodies_.size());
    Fp = Zero;
    Fv = Zero;
    const scalarField& V = mesh_.V();

    for (const label celli : wallCells_)
    {
        const label bodyi = max(cellBody_[celli], 0);
        if (bodyi >= bodies_.size()) continue;
        // pressure on the wall segment, p_P Sw (Sw points into the body)
        Fp[bodyi] += p[celli]*Sw_[celli];
        // viscous traction, nu (u_P - u_b)/d Awall
        if (moving_)
        {
            Fv[bodyi] += nu[celli]*noSlipCoeff_[celli]*V[celli]
                        *(U[celli] - UbPtr_()[celli]);
        }
        else
        {
            Fv[bodyi] += nu[celli]*noSlipCoeff_[celli]*V[celli]*U[celli];
        }
    }
    forAll(Fp, i)
    {
        reduce(Fp[i], sumOp<vector>());
        reduce(Fv[i], sumOp<vector>());
    }
}


void Foam::cutCellGeometry::reportForces
(
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& nu
) const
{
    List<vector> Fp, Fv;
    forces(p, U, nu, Fp, Fv);

    forAll(bodies_, bodyi)
    {
        Info<< "IBM force on " << bodies_[bodyi].name()
            << ": pressure " << Fp[bodyi]
            << " viscous " << Fv[bodyi]
            << " total " << Fp[bodyi] + Fv[bodyi] << endl;
    }

    if (writeForces_ && Pstream::master())
    {
        if (!forcesFilePtr_)
        {
            const fileName dir
            (
                mesh_.time().globalPath()/"postProcessing"/"ibmForces"
               /mesh_.time().timeName()
            );
            mkDir(dir);
            forcesFilePtr_.reset(new OFstream(dir/"forces.dat"));
            OFstream& os = forcesFilePtr_();
            os  << "# Forces per unit density on the immersed bodies" << nl
                << "# time";
            forAll(bodies_, bodyi)
            {
                os  << "  " << bodies_[bodyi].name()
                    << ":Fp(x y z)  " << bodies_[bodyi].name()
                    << ":Fv(x y z)";
            }
            os  << nl;
        }
        OFstream& os = forcesFilePtr_();
        os  << mesh_.time().timeName();
        forAll(bodies_, bodyi)
        {
            os  << ' ' << Fp[bodyi].x() << ' ' << Fp[bodyi].y() << ' ' << Fp[bodyi].z()
                << ' ' << Fv[bodyi].x() << ' ' << Fv[bodyi].y() << ' ' << Fv[bodyi].z();
        }
        os  << nl;
        os.flush();
    }
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
}
