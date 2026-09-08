#include "ibmBody.H"
#include "fvMesh.H"
#include "searchableSurface.H"
#include "quaternion.H"
#include <cmath>

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
    invert_(dict.getOrDefault<bool>("invert", false)),
    motion_(bodyMotion::New(dict)),
    periodicLengths_(Zero),
    t_(0),
    displacement_(Zero),
    velocity_(Zero),
    angle_(0),
    axis_(Zero)
{
    if (moving())
    {
        Info<< "    body " << name << ": moving, motion type "
            << motion_->type() << endl;
    }
}


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
    if (needsTransform())
    {
        pointField bodyPts;
        toBodyFrame(pts, bodyPts);
        levelSet(bodyPts, phi);
    }
    else
    {
        levelSet(pts, phi);
    }
    if (invert_)
    {
        phi = -phi;
    }
}


void Foam::ibmBody::moveTo(const scalar t)
{
    t_ = t;
    displacement_ = motion_->displacement(t);
    velocity_ = motion_->velocity(t);
    const vector& omega = motion_->omega();
    if (mag(omega) > 0)
    {
        axis_ = omega/mag(omega);
        angle_ = motion_->angle(t);
    }
    else
    {
        axis_ = Zero;
        angle_ = 0;
    }
}


bool Foam::ibmBody::needsTransform() const
{
    return
        mag(displacement_) > 0
     || angle_ != 0
     || cmptMax(periodicLengths_) > 0;
}


void Foam::ibmBody::toBodyFrame
(
    const pointField& pts,
    pointField& out
) const
{
    const point c(currentCentre());
    out.setSize(pts.size());

    // Minimum image about the current centre in the periodic directions
    forAll(pts, i)
    {
        vector r(pts[i] - c);
        for (direction d = 0; d < 3; ++d)
        {
            const scalar L = periodicLengths_[d];
            if (L > 0)
            {
                r[d] -= L*std::floor(r[d]/L + 0.5);
            }
        }
        out[i] = c + r;
    }

    // Undo the translation
    if (mag(displacement_) > 0)
    {
        out -= displacement_;
    }

    // Undo the rotation about the initial centre
    if (angle_ != 0)
    {
        const quaternion q(axis_, -angle_);
        const point c0(centre());
        forAll(out, i)
        {
            out[i] = c0 + q.transform(out[i] - c0);
        }
    }
}


void Foam::ibmBody::velocity(const pointField& pts, vectorField& U) const
{
    U.setSize(pts.size());
    const vector& omega = motion_->omega();
    if (mag(omega) > 0)
    {
        const point c(currentCentre());
        forAll(pts, i)
        {
            vector r(pts[i] - c);
            for (direction d = 0; d < 3; ++d)
            {
                const scalar L = periodicLengths_[d];
                if (L > 0)
                {
                    r[d] -= L*std::floor(r[d]/L + 0.5);
                }
            }
            U[i] = velocity_ + (omega ^ r);
        }
    }
    else
    {
        U = velocity_;
    }
}
