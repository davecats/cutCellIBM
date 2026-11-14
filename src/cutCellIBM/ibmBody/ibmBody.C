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

const Foam::Enum<Foam::ibmBody::wallFunctionType>
Foam::ibmBody::wallFunctionNames
({
    { wallFunctionType::NOWALLFUNCTION, "none" },
    { wallFunctionType::SPALARTALLMARAS, "spalartAllmaras" },
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
    wallFunction_
    (
        wallFunctionNames.getOrDefault
        (
            "wallFunction", dict, wallFunctionType::SPALARTALLMARAS
        )
    ),
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
    else if (invariant())
    {
        Info<< "    body " << name << ": invariant motion of type "
            << motion_->type() << ", angular velocity " << motion_->omega()
            << ": static geometry with a wall velocity" << endl;
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
    if (invariant())
    {
        return cmptMax(periodicLengths_) > 0;
    }
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
    if (invariant())
    {
        toBodyFrame(pts, Zero, 0, Zero, out);
    }
    else
    {
        toBodyFrame(pts, displacement_, angle_, axis_, out);
    }
}


void Foam::ibmBody::levelSetChange
(
    const pointField& pts,
    const scalar tau,
    scalarField& dphi
) const
{
    const vector& omega = motion_->omega();
    const vector disp(displacement_ + tau*velocity_);
    const scalar angle = angle_ + tau*mag(omega);
    const vector axis(mag(omega) > 0 ? omega/mag(omega) : vector::zero);

    pointField bodyPts;
    toBodyFrame(pts, Zero, 0, Zero, bodyPts);
    scalarField phi0(pts.size());
    levelSet(bodyPts, phi0);

    toBodyFrame(pts, disp, (mag(omega) > 0 ? angle : 0), axis, bodyPts);
    dphi.setSize(pts.size());
    levelSet(bodyPts, dphi);

    dphi -= phi0;
    if (invert_)
    {
        dphi = -dphi;
    }
}


void Foam::ibmBody::toBodyFrame
(
    const pointField& pts,
    const vector& disp,
    const scalar angle,
    const vector& axis,
    pointField& out
) const
{
    const point c(referencePoint() + disp);
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
    if (mag(disp) > 0)
    {
        out -= disp;
    }

    // Undo the rotation about the initial centre
    if (angle != 0)
    {
        const quaternion q(axis, -angle);
        const point c0(referencePoint());
        forAll(out, i)
        {
            out[i] = c0 + q.transform(out[i] - c0);
        }
    }
}


Foam::vector Foam::ibmBody::relativePosition(const point& x) const
{
    vector r(x - currentCentre());
    for (direction d = 0; d < 3; ++d)
    {
        const scalar L = periodicLengths_[d];
        if (L > 0)
        {
            r[d] -= L*std::floor(r[d]/L + 0.5);
        }
    }
    return r;
}


void Foam::ibmBody::velocity(const pointField& pts, vectorField& U) const
{
    U.setSize(pts.size());
    const vector& omega = motion_->omega();
    if (mag(omega) > 0)
    {
        forAll(pts, i)
        {
            U[i] = velocity_ + (omega ^ relativePosition(pts[i]));
        }
    }
    else
    {
        U = velocity_;
    }
}
