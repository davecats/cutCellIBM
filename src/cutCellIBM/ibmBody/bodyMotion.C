#include "bodyMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

namespace Foam
{
    defineTypeNameAndDebug(bodyMotion, 0);
    defineRunTimeSelectionTable(bodyMotion, dictionary);
    defineTypeNameAndDebug(noMotion, 0);
    addToRunTimeSelectionTable(bodyMotion, noMotion, dictionary);
    defineTypeNameAndDebug(translatingMotion, 0);
    addToRunTimeSelectionTable(bodyMotion, translatingMotion, dictionary);
    defineTypeNameAndDebug(oscillatingMotion, 0);
    addToRunTimeSelectionTable(bodyMotion, oscillatingMotion, dictionary);
    defineTypeNameAndDebug(prescribedMotion, 0);
    addToRunTimeSelectionTable(bodyMotion, prescribedMotion, dictionary);
    defineTypeNameAndDebug(externalMotion, 0);
    addToRunTimeSelectionTable(bodyMotion, externalMotion, dictionary);
}


Foam::bodyMotion::bodyMotion(const dictionary& dict)
:
    dict_(dict),
    omega_(dict.getOrDefault<vector>("angularVelocity", Zero)),
    hasCentreOfRotation_(dict.found("centreOfRotation")),
    centreOfRotation_(dict.getOrDefault<point>("centreOfRotation", Zero))
{}


Foam::autoPtr<Foam::bodyMotion> Foam::bodyMotion::New
(
    const dictionary& bodyDict
)
{
    if (!bodyDict.found("motion"))
    {
        return autoPtr<bodyMotion>(new noMotion(dictionary()));
    }

    const dictionary& dict = bodyDict.subDict("motion");
    const word type(dict.getOrDefault<word>("type", "none"));

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(type);
    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown motion type " << type << nl
            << "Valid types: " << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }
    return autoPtr<bodyMotion>(cstrIter()(dict));
}


void Foam::bodyMotion::setState(const vector&, const vector&)
{
    FatalErrorInFunction
        << "setState is only available for motion type external"
        << exit(FatalError);
}


// none

Foam::noMotion::noMotion(const dictionary& dict) : bodyMotion(dict) {}
bool Foam::noMotion::moving() const { return mag(omega_) > 0; }
Foam::vector Foam::noMotion::displacement(const scalar) const { return Zero; }
Foam::vector Foam::noMotion::velocity(const scalar) const { return Zero; }


// translating

Foam::translatingMotion::translatingMotion(const dictionary& dict)
:
    bodyMotion(dict),
    U_(dict.get<vector>("velocity"))
{}

Foam::vector Foam::translatingMotion::displacement(const scalar t) const
{
    return U_*t;
}

Foam::vector Foam::translatingMotion::velocity(const scalar) const
{
    return U_;
}


// oscillating

Foam::oscillatingMotion::oscillatingMotion(const dictionary& dict)
:
    bodyMotion(dict),
    amplitude_(dict.get<vector>("amplitude")),
    frequency_(dict.get<scalar>("frequency")),
    phase_(dict.getOrDefault<scalar>("phase", 0))
{}

Foam::vector Foam::oscillatingMotion::displacement(const scalar t) const
{
    using constant::mathematical::twoPi;
    return amplitude_*sin(twoPi*frequency_*t + phase_);
}

Foam::vector Foam::oscillatingMotion::velocity(const scalar t) const
{
    using constant::mathematical::twoPi;
    return twoPi*frequency_*amplitude_*cos(twoPi*frequency_*t + phase_);
}


// prescribed

Foam::prescribedMotion::prescribedMotion(const dictionary& dict)
:
    bodyMotion(dict),
    displacement_(Function1<vector>::New("displacement", dict)),
    velocity_(Function1<vector>::New("velocity", dict))
{}

Foam::vector Foam::prescribedMotion::displacement(const scalar t) const
{
    return displacement_->value(t);
}

Foam::vector Foam::prescribedMotion::velocity(const scalar t) const
{
    return velocity_->value(t);
}


// external

Foam::externalMotion::externalMotion(const dictionary& dict)
:
    bodyMotion(dict),
    displacement_(dict.getOrDefault<vector>("displacement", Zero)),
    velocity_(dict.getOrDefault<vector>("velocity", Zero))
{}

Foam::vector Foam::externalMotion::displacement(const scalar) const
{
    return displacement_;
}

Foam::vector Foam::externalMotion::velocity(const scalar) const
{
    return velocity_;
}

void Foam::externalMotion::setState
(
    const vector& displacement,
    const vector& velocity
)
{
    displacement_ = displacement;
    velocity_ = velocity;
}
