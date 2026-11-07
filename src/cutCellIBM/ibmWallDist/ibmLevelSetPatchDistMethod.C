#include "ibmLevelSetPatchDistMethod.H"
#include "meshWavePatchDistMethod.H"
#include "cutCellGeometry.H"
#include "fvMesh.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace patchDistMethods
{
    defineTypeNameAndDebug(ibmLevelSet, 0);
    addToRunTimeSelectionTable(patchDistMethod, ibmLevelSet, dictionary);
}
}


Foam::patchDistMethods::ibmLevelSet::ibmLevelSet
(
    const dictionary& dict,
    const fvMesh& mesh,
    const labelHashSet& patchIDs
)
:
    patchDistMethod(mesh, patchIDs)
{}


bool Foam::patchDistMethods::ibmLevelSet::movePoints()
{
    const cutCellGeometry* ibmPtr =
        mesh_.cfindObject<cutCellGeometry>("ibmProperties");
    return !ibmPtr || ibmPtr->upToDate();
}


bool Foam::patchDistMethods::ibmLevelSet::correct(volScalarField& y)
{
    const cutCellGeometry* ibmPtr =
        mesh_.cfindObject<cutCellGeometry>("ibmProperties");
    if (!ibmPtr)
    {
        FatalErrorInFunction
            << "wallDist method ibmLevelSet needs the cut-cell geometry "
            << "(constant/ibmProperties), constructed before the turbulence "
            << "model" << exit(FatalError);
    }
    const cutCellGeometry& ibm = *ibmPtr;

    // Distance to the wall patches of the mesh, if there are any
    if (patchIDs_.size())
    {
        meshWave(mesh_, patchIDs_).correct(y);
    }
    else
    {
        y = dimensionedScalar(dimLength, GREAT);
    }

    // Distance to the immersed walls at the fluid centroids: the centre of
    // a cut cell may lie in the solid, its fluid centroid never does (for a
    // planar wall its distance is the wall distance alpha V/(2 Awall) of the
    // first-order treatment). Floored in the solid (blanked cells, where the
    // turbulence models must not divide by zero).
    const scalar floor = 1e-9*mag(mesh_.bounds().span());

    scalarField d;
    ibm.wallDistance(ibm.fluidCentroids(), d);

    // Cells of a merged group (small cut cells linked to a master) are one
    // control volume: every cell of the group gets the distance of its
    // centroid
    if (ibm.linking())
    {
        const scalarField& yGroup = ibm.linkWallDist();
        forAll(d, celli)
        {
            if (yGroup[celli] > 0)
            {
                d[celli] = yGroup[celli];
            }
        }
    }

    scalarField& yI = y.primitiveFieldRef();
    forAll(yI, celli)
    {
        yI[celli] = min(yI[celli], max(d[celli], floor));
    }

    forAll(mesh_.boundary(), patchi)
    {
        if (patchIDs_.found(patchi)) continue;
        fvPatchScalarField& yp = y.boundaryFieldRef()[patchi];
        if (yp.empty()) continue;
        scalarField dp;
        ibm.wallDistance(mesh_.boundary()[patchi].Cf(), dp);
        forAll(yp, i)
        {
            yp[i] = min(yp[i], max(dp[i], floor));
        }
    }

    return true;
}


bool Foam::patchDistMethods::ibmLevelSet::correct
(
    volScalarField& y,
    volVectorField& n
)
{
    correct(y);
    n = fvc::grad(y);
    n /= mag(n) + dimensionedScalar(dimless, SMALL);
    n.correctBoundaryConditions();
    return true;
}
