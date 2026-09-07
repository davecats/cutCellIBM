#include "ibmBody.H"
#include "fvMesh.H"
#include "searchableSurface.H"

namespace Foam
{
    defineTypeNameAndDebug(ibmBody, 0);
    defineRunTimeSelectionTable(ibmBody, dictionary);
}

const Foam::Enum<Foam::ibmBody::wallType> Foam::ibmBody::wallTypeNames
({
    { wallType::NOSLIP, "noSlip" },
    { wallType::FREESLIP, "freeSlip" },
});


Foam::ibmBody::ibmBody
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    dict_(dict),
    wallType_(wallTypeNames.getOrDefault("wallType", dict, wallType::NOSLIP)),
    invert_(dict.getOrDefault<bool>("invert", false))
{}


Foam::autoPtr<Foam::ibmBody> Foam::ibmBody::New
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
{
    const word type(dict.get<word>("type"));

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(type);

    if (cstrIter.found())
    {
        Info<< "    body " << name << ": type " << type << endl;
        return autoPtr<ibmBody>(cstrIter()(name, mesh, dict));
    }

    // Fall back on the searchableSurface run-time selection table, so that
    // "type triSurfaceMesh;" (STL) or any analytical searchable surface can
    // be given directly.
    if
    (
        searchableSurface::dictConstructorTablePtr_
     && searchableSurface::dictConstructorTablePtr_->found(type)
    )
    {
        Info<< "    body " << name << ": searchableSurface of type " << type
            << endl;
        return autoPtr<ibmBody>
        (
            dictionaryConstructorTablePtr_->cfind("searchableSurface")()
            (name, mesh, dict)
        );
    }

    FatalIOErrorInFunction(dict)
        << "Unknown ibmBody type " << type << nl << nl
        << "Valid types: " << dictionaryConstructorTablePtr_->sortedToc()
        << " or any searchableSurface type "
        << searchableSurface::dictConstructorTablePtr_->sortedToc()
        << exit(FatalIOError);

    return nullptr;
}


void Foam::ibmBody::levelSet(const pointField& pts, scalarField& phi) const
{
    NotImplemented;
}


void Foam::ibmBody::evaluate(const pointField& pts, scalarField& phi) const
{
    phi.setSize(pts.size());
    levelSet(pts, phi);
    if (invert_)
    {
        phi = -phi;
    }
}
