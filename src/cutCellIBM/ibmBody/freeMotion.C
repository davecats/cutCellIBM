#include "freeMotion.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(freeMotion, 0);
    addToRunTimeSelectionTable(bodyMotion, freeMotion, dictionary);
}

const Foam::Enum<Foam::freeMotion::couplingType>
Foam::freeMotion::couplingNames
({
    { couplingType::EXPLICIT, "explicit" },
    { couplingType::IMPLICIT, "implicit" },
});


namespace
{
    // A direction mask: every component 0 or 1
    Foam::vector readMask
    (
        const Foam::dictionary& dict,
        const Foam::word& key,
        const Foam::vector& def
    )
    {
        const Foam::vector m(dict.getOrDefault<Foam::vector>(key, def));
        for (Foam::direction i = 0; i < Foam::vector::nComponents; ++i)
        {
            if (m[i] != 0 && m[i] != 1)
            {
                FatalIOErrorInFunction(dict)
                    << key << " " << m << ": the components must be 0 or 1"
                    << Foam::exit(Foam::FatalIOError);
            }
        }
        return m;
    }
}


Foam::freeMotion::freeMotion(const dictionary& dict)
:
    bodyMotion(dict),
    density_(dict.get<scalar>("density")),
    translationMask_(readMask(dict, "translation", vector::one)),
    rotationMask_(readMask(dict, "rotation", Zero)),
    coupling_
    (
        couplingNames.getOrDefault("coupling", dict, couplingType::IMPLICIT)
    ),
    addedMass_(dict.getOrDefault<scalar>("addedMass", 0.5)),
    volume_(dict.getOrDefault<scalar>("volume", -1)),
    inertia_(dict.getOrDefault<scalar>("inertia", -1)),
    gravity_(dict.getOrDefault<vector>("gravity", Zero)),
    externalForce_(dict.getOrDefault<vector>("externalForce", Zero)),
    meanPressureGradient_
    (
        dict.getOrDefault<bool>("meanPressureGradient", true)
    ),
    displacement_(Zero),
    velocity_(dict.getOrDefault<vector>("velocity", Zero)),
    theta_(0),
    displacement0_(Zero),
    velocity0_(velocity_),
    omega0_(Zero),
    theta0_(0),
    residual_(0),
    ma_(Zero),
    Ia_(0),
    hasPrevious_(false),
    Uprev_(Zero),
    Fprev_(Zero),
    Oprev_(0),
    Mprev_(0)
{
    if (invariant_)
    {
        FatalIOErrorInFunction(dict)
            << "a free motion cannot be declared invariant (the body "
            << "translates)" << exit(FatalIOError);
    }
    if (density_ < 0 || addedMass_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "density and addedMass must not be negative"
            << exit(FatalIOError);
    }
    if (cmptSum(rotationMask_) > 1)
    {
        FatalIOErrorInFunction(dict)
            << "rotation " << rotationMask_ << ": at most one free rotation "
            << "axis (rotation about a fixed axis)" << exit(FatalIOError);
    }
    if (coupling_ == EXPLICIT && density_ == 0)
    {
        FatalIOErrorInFunction(dict)
            << "explicit coupling needs a body mass (density > 0)"
            << exit(FatalIOError);
    }
    if (coupling_ == IMPLICIT && density_ == 0 && addedMass_ == 0)
    {
        FatalIOErrorInFunction(dict)
            << "a massless body needs addedMass > 0" << exit(FatalIOError);
    }
    // The initial angular velocity about the free axis only
    omega_ = cmptMultiply(omega_, rotationMask_);
    omega0_ = omega_;
}


Foam::vector Foam::freeMotion::displacement(const scalar) const
{
    return displacement_;
}


Foam::vector Foam::freeMotion::velocity(const scalar) const
{
    return velocity_;
}


Foam::vector Foam::freeMotion::axis() const
{
    return rotationMask_;
}


void Foam::freeMotion::setMassProperties
(
    const scalar volume,
    const scalar inertia
)
{
    if (volume_ <= 0)
    {
        volume_ = volume;
    }
    if (inertia_ < 0 && volume > 0)
    {
        // the body's radius of gyration with the volume in use
        inertia_ = inertia/volume*volume_;
    }
    if (volume_ <= 0)
    {
        FatalIOErrorInFunction(dict_)
            << "the body gives no volume: set volume (and inertia if it "
            << "rotates) in the motion dictionary" << exit(FatalIOError);
    }
    if (rotates() && inertia_ <= 0)
    {
        FatalIOErrorInFunction(dict_)
            << "the body gives no moment of inertia: set inertia in the "
            << "motion dictionary" << exit(FatalIOError);
    }
    if (inertia_ < 0)
    {
        inertia_ = 0;
    }
}


void Foam::freeMotion::readState(const dictionary& dict)
{
    displacement_ = dict.get<vector>("displacement");
    velocity_ = dict.get<vector>("velocity");
    theta_ = dict.get<scalar>("angle");
    omega_ = dict.get<vector>("angularVelocity");
    displacement0_ = displacement_;
    velocity0_ = velocity_;
    theta0_ = theta_;
    omega0_ = omega_;
}


Foam::dictionary Foam::freeMotion::state() const
{
    dictionary d;
    d.add("displacement", displacement_);
    d.add("velocity", velocity_);
    d.add("angle", theta_);
    d.add("angularVelocity", omega_);
    return d;
}


void Foam::freeMotion::advance
(
    const scalar dt,
    const vector& F,
    const vector& M
)
{
    const scalar m = mass();
    const vector Uk(velocity_);
    const vector Ok(omega_);
    const vector e(rotates() ? axis() : vector::zero);
    const scalar I = density_*inertia_;
    const scalar Mk = (M & e);
    const scalar wk = (Ok & e);

    if (implicit())
    {
        // Bounds of the effective added mass: the given coefficient, and a
        // cap against a secant from two nearly equal iterates
        const scalar maMin = addedMass_*volume_;
        const scalar IaMin = addedMass_*inertia_;
        // Secant estimate of -dt dF/dU. The estimate in use decays by the
        // factor DECAY per update and is raised at once by a larger secant:
        // the loads of an iterate carry the noise of an unconverged flow and
        // jump where cells change their status (the geometry moves with the
        // velocity), and an estimate below the true one is what makes the
        // iteration overshoot (one above it only slows it down). Secants
        // from velocity differences below SECANT of the velocity scale are
        // noise and not taken; without them the estimate keeps decaying, so
        // that a large one cannot freeze the body. Capped at CAP (m + m_min).
        static const scalar DECAY = 0.8;
        static const scalar SECANT = 1e-4;
        static const scalar CAP = 1e3;
        const scalar k2 = max(inertia_, 0)/volume_;
        const scalar Uscale = max(mag(Uk), Foam::sqrt(k2)*mag(Ok));
        forAll(translationMask_, i)
        {
            if (translationMask_[i] == 0) continue;
            scalar ma = max(DECAY*ma_[i], maMin);
            const scalar dU = Uk[i] - Uprev_[i];
            if (hasPrevious_ && mag(dU) > SECANT*max(Uscale, VSMALL))
            {
                const scalar est = -dt*(F[i] - Fprev_[i])/dU;
                ma = max(ma, min(est, CAP*(m + maMin)));
            }
            ma_[i] = ma;
        }
        if (rotates())
        {
            scalar Ia = max(DECAY*Ia_, IaMin);
            const scalar dw = wk - Oprev_;
            if (hasPrevious_ && mag(dw) > SECANT*max(mag(wk), VSMALL))
            {
                const scalar est = -dt*(Mk - Mprev_)/dw;
                Ia = max(Ia, min(est, CAP*(I + IaMin)));
            }
            Ia_ = Ia;
        }
        Uprev_ = Uk;
        Fprev_ = F;
        Oprev_ = wk;
        Mprev_ = Mk;
        hasPrevious_ = true;
    }

    // Translation, in the free directions: (m + m_a) U^{k+1} = m U^n +
    // m_a U^k + dt F (explicit: m_a = 0); the others keep their value
    forAll(translationMask_, i)
    {
        if (translationMask_[i] > 0)
        {
            const scalar ma = implicit() ? ma_[i] : 0;
            velocity_[i] = (m*velocity0_[i] + ma*Uk[i] + dt*F[i])/(m + ma);
        }
    }
    displacement_ = displacement0_ + dt*velocity_;

    // Rotation about the free axis, the same scheme with I and I_a
    if (rotates())
    {
        const scalar Ia = implicit() ? Ia_ : 0;
        const scalar w = (I*(omega0_ & e) + Ia*wk + dt*Mk)/(I + Ia);
        omega_ = w*e;
        theta_ = theta0_ + dt*w;
    }

    // Change of the iterate against the velocity scale of the body
    const scalar k = (volume_ > 0 ? Foam::sqrt(max(inertia_, 0)/volume_) : 0);
    const scalar dU = max(mag(velocity_ - Uk), k*mag(omega_ - Ok));
    const scalar Uref = max(max(mag(velocity_), k*mag(omega_)), VSMALL);
    residual_ = dU/Uref;
}


void Foam::freeMotion::predict
(
    const scalar dt,
    const vector& F,
    const vector& M
)
{
    displacement0_ = displacement_;
    velocity0_ = velocity_;
    omega0_ = omega_;
    theta0_ = theta_;
    // the secant pairs belong to one time step; the added mass is kept
    hasPrevious_ = false;
    advance(dt, F, M);
}


void Foam::freeMotion::correct
(
    const scalar dt,
    const vector& F,
    const vector& M
)
{
    advance(dt, F, M);
}
