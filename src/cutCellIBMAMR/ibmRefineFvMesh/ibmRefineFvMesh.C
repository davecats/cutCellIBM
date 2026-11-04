#include "ibmRefineFvMesh.H"
#include "cutCellGeometry.H"
#include "polyTopoChange.H"
#include "surfaceFields.H"
#include "syncTools.H"
#include "wallDist.H"
#include "processorPolyPatch.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(ibmRefineFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, ibmRefineFvMesh, IOobject);
    addToRunTimeSelectionTable(dynamicFvMesh, ibmRefineFvMesh, doInit);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::dictionary Foam::ibmRefineFvMesh::coeffs() const
{
    return IOdictionary
    (
        IOobject
        (
            "dynamicMeshDict",
            time().constant(),
            *this,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    ).optionalSubDict(typeName + "Coeffs");
}


Foam::cutCellGeometry& Foam::ibmRefineFvMesh::ibm() const
{
    const cutCellGeometry* ibmPtr =
        cfindObject<cutCellGeometry>("ibmProperties");
    if (!ibmPtr)
    {
        FatalErrorInFunction
            << typeName << " needs the cut-cell geometry (constant/"
            << "ibmProperties), constructed before the first update"
            << exit(FatalError);
    }
    return const_cast<cutCellGeometry&>(*ibmPtr);
}


Foam::scalarField Foam::ibmRefineFvMesh::cellRadius() const
{
    // Distances in the solved directions only: the thickness of a 2-D
    // mesh does not count
    vector mask(Zero);
    for (direction d = 0; d < vector::nComponents; ++d)
    {
        if (solutionD()[d] == 1) mask[d] = 1;
    }

    const pointField& pts = points();
    const vectorField& C = cellCentres();
    scalarField r(nCells(), Zero);
    forAll(r, celli)
    {
        for (const label pointi : cellPoints()[celli])
        {
            r[celli] = max(r[celli], mag(cmptMultiply(mask, pts[pointi] - C[celli])));
        }
    }
    return r;
}


bool Foam::ibmRefineFvMesh::regionField
(
    const word& name,
    scalarField& value
) const
{
    if (const auto* fp = cfindObject<volScalarField>(name))
    {
        value = fp->primitiveField();
        return true;
    }
    if (const auto* fp = cfindObject<volVectorField>(name))
    {
        value = mag(fp->primitiveField());
        return true;
    }
    if (const auto* fp = cfindObject<volSymmTensorField>(name))
    {
        value = mag(fp->primitiveField());
        return true;
    }
    if (const auto* fp = cfindObject<volTensorField>(name))
    {
        value = mag(fp->primitiveField());
        return true;
    }
    return false;
}


void Foam::ibmRefineFvMesh::dilate(labelList& level, const label nLayers) const
{
    const labelUList& own = faceOwner();
    const labelUList& nei = faceNeighbour();

    for (label layer = 0; layer < nLayers; ++layer)
    {
        labelList neiLevel;
        syncTools::swapBoundaryCellList(*this, level, neiLevel);

        labelList newLevel(level);
        forAll(nei, facei)
        {
            newLevel[own[facei]] = max(newLevel[own[facei]], level[nei[facei]]);
            newLevel[nei[facei]] = max(newLevel[nei[facei]], level[own[facei]]);
        }
        forAll(neiLevel, bFacei)
        {
            const label celli = own[bFacei + nInternalFaces()];
            newLevel[celli] = max(newLevel[celli], neiLevel[bFacei]);
        }
        level.transfer(newLevel);
    }
}


bool Foam::ibmRefineFvMesh::targetLevels
(
    const dictionary& dict,
    const bool solution,
    labelList& refineLevel,
    labelList& keepLevel
) const
{
    const label maxLevel = dict.get<label>("maxLevel");
    const scalar nCutCellLayers =
        dict.getOrDefault<scalar>("nCutCellLayers", 2);
    const label nBufferLayers = dict.getOrDefault<label>("nBufferLayers", 1);

    const labelList& cellLevel = meshCutter_->cellLevel();
    const scalarField r(cellRadius());
    const scalar sqrtD = Foam::sqrt(scalar(nSolutionD()));

    // Distance from the cell centres to the bodies
    scalarField d;
    ibm().wallDistance(cellCentres(), d);
    d = mag(d);

    refineLevel.setSize(nCells());
    keepLevel.setSize(nCells());
    refineLevel = 0;
    keepLevel = 0;

    // The cut cells and the layers around them: within the cell radius of
    // the wall, plus nCutCellLayers finest cell sizes. A cell keeps its
    // level while its parent (twice the radius, its centre one child
    // radius away) would still be refined.
    forAll(r, celli)
    {
        const scalar hFine =
            2*r[celli]/sqrtD/pow(2.0, max(maxLevel - cellLevel[celli], 0));
        const scalar band = nCutCellLayers*hFine;
        if (d[celli] < r[celli] + band)
        {
            refineLevel[celli] = maxLevel;
        }
        if (d[celli] < 3*r[celli] + band)
        {
            keepLevel[celli] = maxLevel;
        }
    }

    // Refinement regions
    labelList fieldRefine(nCells(), Zero);
    labelList fieldKeep(nCells(), Zero);
    bool anyField = false;
    bool complete = true;

    for (const entry& e : dict.subOrEmptyDict("regions"))
    {
        if (!e.isDict()) continue;
        const dictionary& rd = e.dict();
        const word type(rd.get<word>("type"));
        const label level = min(rd.get<label>("level"), maxLevel);

        if (type == "distance")
        {
            const scalar dist = rd.get<scalar>("distance");
            forAll(r, celli)
            {
                if (d[celli] - r[celli] < dist)
                {
                    refineLevel[celli] = max(refineLevel[celli], level);
                }
                if (d[celli] - 3*r[celli] < dist)
                {
                    keepLevel[celli] = max(keepLevel[celli], level);
                }
            }
        }
        else if (type == "field")
        {
            if (!solution) continue;

            const word name(rd.get<word>("field"));
            const scalar lower = rd.get<scalar>("lowerBound");
            const scalar upper = rd.get<scalar>("upperBound");
            scalarField value;
            if (!regionField(name, value))
            {
                Info<< "    region " << e.keyword() << ": no field " << name
                    << " (yet)" << endl;
                complete = false;
                continue;
            }
            anyField = true;
            forAll(value, celli)
            {
                if (value[celli] >= lower && value[celli] <= upper)
                {
                    fieldRefine[celli] = max(fieldRefine[celli], level);
                }
            }
        }
        else
        {
            FatalIOErrorInFunction(rd)
                << "Unknown region type " << type << " of region "
                << e.keyword() << ", expected distance or field"
                << exit(FatalIOError);
        }
    }

    // Solution-based regions with a buffer, and a wider one for keeping
    if (returnReduce(anyField, orOp<bool>()))
    {
        fieldKeep = fieldRefine;
        dilate(fieldRefine, nBufferLayers);
        dilate(fieldKeep, nBufferLayers + 1);
        forAll(refineLevel, celli)
        {
            refineLevel[celli] = max(refineLevel[celli], fieldRefine[celli]);
            keepLevel[celli] = max(keepLevel[celli], fieldKeep[celli]);
        }
    }

    forAll(keepLevel, celli)
    {
        keepLevel[celli] = max(keepLevel[celli], refineLevel[celli]);
    }

    return returnReduce(complete, andOp<bool>());
}


Foam::label Foam::ibmRefineFvMesh::changeTopology(polyTopoChange& meshMod)
{
    const vectorField SfOld(faceAreas());

    // The cut-cell geometry belongs to the old mesh until it is rebuilt
    ibm().outOfDate();

    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);
    updateMesh(map());

    // Faces that are new or whose area vector has changed (split, merged,
    // flipped), since the start of the adaptation
    const labelList& faceMap = map().faceMap();
    const vectorField& Sf = faceAreas();
    bitSet changed(nFaces());
    forAll(faceMap, facei)
    {
        const label oldFacei = faceMap[facei];
        if
        (
            oldFacei < 0
         || changedFaces_.test(oldFacei)
         || mag(Sf[facei] - SfOld[oldFacei]) > 1e-9*mag(Sf[facei])
        )
        {
            changed.set(facei);
        }
    }
    changedFaces_.transfer(changed);

    if (debug)
    {
        checkMesh(true);
        meshCutter_->checkRefinementLevels(-1, labelList());
        checkProcessorPatches();
    }

    return returnReduce(map().nOldCells(), sumOp<label>());
}


void Foam::ibmRefineFvMesh::refine(const labelList& cells)
{
    polyTopoChange meshMod(*this);
    meshCutter_->setRefinement(cells, meshMod);
    const label nOld = changeTopology(meshMod);

    Info<< "    refined from " << nOld << " to "
        << globalData().nTotalCells() << " cells" << endl;
}


void Foam::ibmRefineFvMesh::unrefine(const labelList& splitElems)
{
    polyTopoChange meshMod(*this);
    meshCutter_->setUnrefinement(splitElems, meshMod);
    const label nOld = changeTopology(meshMod);

    Info<< "    unrefined from " << nOld << " to "
        << globalData().nTotalCells() << " cells" << endl;
}


void Foam::ibmRefineFvMesh::correctFluxes(const dictionary& dict)
{
    // The mapped flux is kept where the faces have not changed; on the new,
    // split and merged faces it is the cut-cell flux of the mapped velocity,
    // at every stored time level
    const List<Pair<word>> fluxes
    (
        dict.getOrDefault<List<Pair<word>>>
        (
            "correctFluxes",
            List<Pair<word>>(1, Pair<word>("phi", "U"))
        )
    );

    for (const Pair<word>& names : fluxes)
    {
        surfaceScalarField* phiPtr = getObjectPtr<surfaceScalarField>(names.first());
        const volVectorField* UPtr = cfindObject<volVectorField>(names.second());
        if (!phiPtr || !UPtr) continue;

        surfaceScalarField* phi = phiPtr;
        const volVectorField* U = UPtr;
        while (true)
        {
            const tmp<surfaceScalarField> tphiU(ibm().flux(*U));
            const surfaceScalarField& phiU = tphiU();

            for (label facei = 0; facei < nInternalFaces(); ++facei)
            {
                if (changedFaces_.test(facei)) (*phi)[facei] = phiU[facei];
            }
            forAll(phi->boundaryField(), patchi)
            {
                fvsPatchScalarField& pp = phi->boundaryFieldRef()[patchi];
                const label start = boundaryMesh()[patchi].start();
                forAll(pp, i)
                {
                    if (changedFaces_.test(start + i))
                    {
                        pp[i] = phiU.boundaryField()[patchi][i];
                    }
                }
            }

            if (!phi->nOldTimes() || !U->nOldTimes()) break;
            phi = &phi->oldTime();
            U = &U->oldTime();
        }
    }
}


void Foam::ibmRefineFvMesh::checkProcessorPatches() const
{
    if (!Pstream::parRun()) return;

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    for (const polyPatch& pp : boundaryMesh())
    {
        if (isA<processorPolyPatch>(pp))
        {
            const auto& ppp = refCast<const processorPolyPatch>(pp);
            UOPstream to(ppp.neighbProcNo(), pBufs);
            to << ppp.size() << ppp.nPoints() << ppp.nEdges();
        }
    }
    pBufs.finishedSends();
    for (const polyPatch& pp : boundaryMesh())
    {
        if (isA<processorPolyPatch>(pp))
        {
            const auto& ppp = refCast<const processorPolyPatch>(pp);
            UIPstream from(ppp.neighbProcNo(), pBufs);
            label nf, np, ne;
            from >> nf >> np >> ne;
            if (nf != ppp.size() || np != ppp.nPoints() || ne != ppp.nEdges())
            {
                Pout<< "    patch " << pp.name() << ": faces/points/edges "
                    << ppp.size() << ' ' << ppp.nPoints() << ' '
                    << ppp.nEdges() << ", neighbour " << nf << ' ' << np
                    << ' ' << ne << endl;
            }
        }
    }
}


void Foam::ibmRefineFvMesh::checkCutCells(const label maxLevel) const
{
    const labelList& cellLevel = meshCutter_->cellLevel();
    const bitSet& cut = ibm().cutCells();

    // Cut cells below maxLevel, and cut cells with a face neighbour below
    // maxLevel
    labelList neiLevel;
    syncTools::swapBoundaryCellList(*this, cellLevel, neiLevel);

    label nCoarseCut = 0;
    label nCoarseNeighbour = 0;
    for (const label celli : cut)
    {
        if (cellLevel[celli] < maxLevel) ++nCoarseCut;

        for (const label facei : cells()[celli])
        {
            label other = -1;
            if (facei < nInternalFaces())
            {
                other = cellLevel[faceOwner()[facei] == celli
                    ? faceNeighbour()[facei] : faceOwner()[facei]];
            }
            else if (boundaryMesh()[boundaryMesh().whichPatch(facei)].coupled())
            {
                other = neiLevel[facei - nInternalFaces()];
            }
            if (other >= 0 && other < maxLevel)
            {
                ++nCoarseNeighbour;
                break;
            }
        }
    }
    reduce(nCoarseCut, sumOp<label>());
    reduce(nCoarseNeighbour, sumOp<label>());

    if (nCoarseCut || nCoarseNeighbour)
    {
        WarningInFunction
            << nCoarseCut << " cut cells below maxLevel " << maxLevel << ", "
            << nCoarseNeighbour << " with a face neighbour below it: a body "
            << "has moved out of its refined band (increase nCutCellLayers "
            << "or decrease refineInterval)" << endl;
    }
    else
    {
        Info<< "    all cut cells and their face neighbours at level "
            << maxLevel << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ibmRefineFvMesh::ibmRefineFvMesh(const IOobject& io, const bool doInit)
:
    dynamicFvMesh(io, doInit),
    meshCutter_(hexRef::New(*this)),
    nAdapt_(0),
    changedFaces_()
{
    // The engine splits hexahedra: a cell of the background mesh must have
    // six faces of four points each, and a cell refined before (restart)
    // must still have its eight anchor points
    const labelList& cellLevel = meshCutter_->cellLevel();
    const labelList& pointLevel = meshCutter_->pointLevel();

    labelList nAnchors(nCells(), Zero);
    forAll(pointCells(), pointi)
    {
        for (const label celli : pointCells()[pointi])
        {
            if (pointLevel[pointi] <= cellLevel[celli]) ++nAnchors[celli];
        }
    }

    label nBad = 0;
    forAll(nAnchors, celli)
    {
        const cell& cFaces = cells()[celli];
        bool bad = (nAnchors[celli] != 8 || cFaces.size() < 6);
        for (const label facei : cFaces)
        {
            if (faces()[facei].size() < 4) bad = true;
        }
        if (bad) ++nBad;
    }
    reduce(nBad, sumOp<label>());
    if (nBad)
    {
        FatalErrorInFunction
            << nBad << " cells are not (split) hexahedra: " << typeName
            << " needs a hexahedral background mesh (blockMesh)"
            << exit(FatalError);
    }

    const dictionary dict(coeffs());
    Info<< typeName << ": " << nSolutionD() << "-D, maxLevel "
        << dict.get<label>("maxLevel") << ", "
        << dict.getOrDefault<scalar>("nCutCellLayers", 2)
        << " finest cell layers around the walls, level up to "
        << gMax(cellLevel) << " now" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::ibmRefineFvMesh::update()
{
    const dictionary dict(coeffs());

    const label refineInterval = dict.get<label>("refineInterval");
    const label maxLevel = dict.get<label>("maxLevel");
    const label maxCells = dict.getOrDefault<label>("maxCells", labelMax);
    const bool unrefineCells = dict.getOrDefault<bool>("unrefine", true);

    if (maxLevel < 0 || refineInterval < 0)
    {
        FatalIOErrorInFunction(dict)
            << "maxLevel and refineInterval must not be negative"
            << exit(FatalIOError);
    }

    // The first update of the run, then every refineInterval steps
    const label step = time().timeIndex() - time().startTimeIndex();
    const bool due =
        (nAdapt_ == 0)
     || (refineInterval > 0 && step % refineInterval == 0);

    topoChanging(false);
    if (!due)
    {
        return false;
    }

    Info<< typeName << ": adapting the mesh" << endl;
    bool changed = false;
    changedFaces_.reset();
    changedFaces_.resize(nFaces());

    // Refine until every cell has reached its target (one level per pass)
    for (label pass = 0; pass <= maxLevel; ++pass)
    {
        const bool solution = globalData().nTotalCells() < maxCells;

        labelList refineLevel, keepLevel;
        targetLevels(dict, solution, refineLevel, keepLevel);

        const labelList& cellLevel = meshCutter_->cellLevel();
        DynamicList<label> candidates(nCells()/8);
        forAll(cellLevel, celli)
        {
            if (cellLevel[celli] < refineLevel[celli])
            {
                candidates.append(celli);
            }
        }

        // 2:1 balance: adds the cells the candidates need
        const labelList cellsToRefine
        (
            meshCutter_->consistentRefinement(candidates, true)
        );
        if (returnReduce(cellsToRefine.size(), sumOp<label>()) == 0)
        {
            break;
        }

        refine(cellsToRefine);
        changed = true;
    }

    // Coarsen the cells above the level they keep (one level per update),
    // not while the field of a region is missing (e.g. a function object
    // that has not run yet after a restart)
    labelList refineLevel, keepLevel;
    if
    (
        unrefineCells
     && targetLevels(dict, true, refineLevel, keepLevel)
    )
    {
        const labelList& cellLevel = meshCutter_->cellLevel();
        bitSet keep(nCells());
        forAll(cellLevel, celli)
        {
            if (cellLevel[celli] <= keepLevel[celli]) keep.set(celli);
        }

        const labelList splitElems
        (
            meshCutter_->selectUnrefineElems(keep, bitSet(nPoints(), true))
        );
        if (returnReduce(splitElems.size(), sumOp<label>()))
        {
            unrefine(splitElems);
            changed = true;
        }
    }

    if (++nAdapt_ % 10 == 0)
    {
        // Unrefinement leaves holes in the history
        const_cast<refinementHistoryMultiDim&>
        (
            meshCutter_->history()
        ).compact();
    }

    if (changed)
    {
        // The cut-cell geometry, the fluxes on the changed faces, then the
        // wall distance that depends on the geometry
        ibm().updateMesh();
        correctFluxes(dict);
        HashTable<wallDist*> wallDists(lookupClass<wallDist>());
        forAllIters(wallDists, iter)
        {
            iter.val()->movePoints();
        }
        checkCutCells(maxLevel);
    }

    topoChanging(changed);
    if (changed)
    {
        moving(false);
    }

    return changed;
}


void Foam::ibmRefineFvMesh::updateMesh(const mapPolyMesh& mpm)
{
    // Maps every registered field. The wall distance waits for the cut-cell
    // geometry (patchDistMethods::ibmLevelSet::movePoints).
    dynamicFvMesh::updateMesh(mpm);
    meshCutter_->updateMesh(mpm);

    // The solved and geometric directions are recomputed on demand with a
    // global reduction: do it here on every processor, before code that
    // runs on some processors only (the cut cells) asks for them
    (void)solutionD();
    (void)geometricD();
}


bool Foam::ibmRefineFvMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool valid
) const
{
    // Refinement data with the mesh in the current time directory
    const_cast<hexRef&>(meshCutter_()).setInstance(time().timeName());

    volScalarField level
    (
        IOobject
        (
            "refinementLevel",
            time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        *this,
        dimensionedScalar(dimless, Zero),
        "zeroGradient"
    );
    forAll(level, celli)
    {
        level[celli] = meshCutter_->cellLevel()[celli];
    }

    return
        dynamicFvMesh::writeObject(streamOpt, valid)
     && meshCutter_->write(valid)
     && level.write();
}


// ************************************************************************* //
