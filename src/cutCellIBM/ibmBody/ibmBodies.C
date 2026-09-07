#include "ibmBodies.H"
#include "fvMesh.H"
#include "Time.H"
#include "searchableSurface.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(cylinderBody, 0);
    addToRunTimeSelectionTable(ibmBody, cylinderBody, dictionary);
    defineTypeNameAndDebug(sphereBody, 0);
    addToRunTimeSelectionTable(ibmBody, sphereBody, dictionary);
    defineTypeNameAndDebug(boxBody, 0);
    addToRunTimeSelectionTable(ibmBody, boxBody, dictionary);
    defineTypeNameAndDebug(planeBody, 0);
    addToRunTimeSelectionTable(ibmBody, planeBody, dictionary);
    defineTypeNameAndDebug(searchableSurfaceBody, 0);
    addToRunTimeSelectionTable(ibmBody, searchableSurfaceBody, dictionary);
    defineTypeNameAndDebug(alphaFieldBody, 0);
    addToRunTimeSelectionTable(ibmBody, alphaFieldBody, dictionary);
}


// * * * * * * * * * * * * * * * * cylinder  * * * * * * * * * * * * * * * * //

Foam::cylinderBody::cylinderBody
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ibmBody(name, mesh, dict),
    centre_(dict.get<point>("centre")),
    axis_(normalised(dict.get<vector>("axis"))),
    radius_(dict.get<scalar>("radius"))
{}


void Foam::cylinderBody::levelSet
(
    const pointField& pts,
    scalarField& phi
) const
{
    forAll(pts, i)
    {
        const vector d(pts[i] - centre_);
        phi[i] = mag(d - (d & axis_)*axis_) - radius_;
    }
}


// * * * * * * * * * * * * * * * * * sphere  * * * * * * * * * * * * * * * * //

Foam::sphereBody::sphereBody
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ibmBody(name, mesh, dict),
    centre_(dict.get<point>("centre")),
    radius_(dict.get<scalar>("radius"))
{}


void Foam::sphereBody::levelSet
(
    const pointField& pts,
    scalarField& phi
) const
{
    forAll(pts, i)
    {
        phi[i] = mag(pts[i] - centre_) - radius_;
    }
}


// * * * * * * * * * * * * * * * * * * box * * * * * * * * * * * * * * * * * //

Foam::boxBody::boxBody
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ibmBody(name, mesh, dict),
    min_(dict.get<point>("min")),
    max_(dict.get<point>("max"))
{}


void Foam::boxBody::levelSet
(
    const pointField& pts,
    scalarField& phi
) const
{
    // Signed distance to an axis-aligned box (negative inside)
    const point c(0.5*(min_ + max_));
    const vector h(0.5*(max_ - min_));
    forAll(pts, i)
    {
        const vector q(cmptMag(pts[i] - c) - h);
        const vector qp(max(q.x(), 0.0), max(q.y(), 0.0), max(q.z(), 0.0));
        phi[i] = mag(qp) + min(cmptMax(q), 0.0);
    }
}


// * * * * * * * * * * * * * * * * * plane * * * * * * * * * * * * * * * * * //

Foam::planeBody::planeBody
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ibmBody(name, mesh, dict),
    point_(dict.get<point>("point")),
    normal_(normalised(dict.get<vector>("normal")))
{}


void Foam::planeBody::levelSet
(
    const pointField& pts,
    scalarField& phi
) const
{
    forAll(pts, i)
    {
        phi[i] = (pts[i] - point_) & normal_;
    }
}


// * * * * * * * * * * * * * * searchableSurface * * * * * * * * * * * * * * //

Foam::searchableSurfaceBody::searchableSurfaceBody
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ibmBody(name, mesh, dict),
    surfPtr_(nullptr)
{
    // Either "surface { type ...; }" or the surface type directly as "type"
    const dictionary& sdict =
        dict.found("surface") ? dict.subDict("surface") : dict;

    surfPtr_ = searchableSurface::New
    (
        sdict.get<word>("type"),
        IOobject
        (
            name,
            mesh.time().constant(),
            "triSurface",
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        sdict
    );

    if (!surfPtr_().hasVolumeType())
    {
        FatalIOErrorInFunction(dict)
            << "The surface of body " << name
            << " does not provide inside/outside information. "
            << "It must be closed." << exit(FatalIOError);
    }
}


void Foam::searchableSurfaceBody::levelSet
(
    const pointField& pts,
    scalarField& phi
) const
{
    List<pointIndexHit> info;
    surfPtr_().findNearest(pts, scalarField(pts.size(), GREAT), info);

    List<volumeType> vt;
    surfPtr_().getVolumeType(pts, vt);

    forAll(pts, i)
    {
        scalar d = info[i].hit() ? mag(info[i].hitPoint() - pts[i]) : GREAT;
        phi[i] = (vt[i] == volumeType::INSIDE) ? -d : d;
    }
}


// * * * * * * * * * * * * * * * * alphaField  * * * * * * * * * * * * * * * //

Foam::alphaFieldBody::alphaFieldBody
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ibmBody(name, mesh, dict),
    fieldName_(dict.get<word>("field")),
    solidFraction_(dict.getOrDefault<bool>("solidFraction", false))
{}


Foam::tmp<Foam::volScalarField> Foam::alphaFieldBody::alpha() const
{
    Info<< "    body " << name_ << ": reading fluid fraction from field "
        << fieldName_ << (solidFraction_ ? " (stored as solid fraction)" : "")
        << endl;

    tmp<volScalarField> talpha
    (
        new volScalarField
        (
            IOobject
            (
                fieldName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        )
    );

    volScalarField& a = talpha.ref();
    if (solidFraction_)
    {
        a = 1.0 - a;
    }
    if (invert_)
    {
        a = 1.0 - a;
    }
    a = max(min(a, 1.0), 0.0);

    return talpha;
}
