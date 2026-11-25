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
    inertia_(Zero),
    hasInertia_(dict.found("inertia") || dict.found("inertiaTensor")),
    gravity_(dict.getOrDefault<vector>("gravity", Zero)),
    externalForce_(dict.getOrDefault<vector>("externalForce", Zero)),
    meanPressureGradient_
    (
        dict.getOrDefault<bool>("meanPressureGradient", true)
    ),
    displacement_(Zero),
    velocity_(dict.getOrDefault<vector>("velocity", Zero)),
    orientation_(quaternion::I),
    displacement0_(Zero),
    velocity0_(velocity_),
    omega0_(Zero),
    orientation0_(quaternion::I),
    residual_(0),
    ma_(Zero),
    Ia_(Zero),
    hasPrevious_(false),
    Uprev_(Zero),
    Fprev_(Zero),
    Oprev_(Zero),
    Mprev_(Zero)
{
    if (dict.found("inertia") && dict.found("inertiaTensor"))
    {
        FatalIOErrorInFunction(dict)
            << "give either inertia (isotropic) or inertiaTensor"
            << exit(FatalIOError);
    }
    if (dict.found("inertia"))
    {
        inertia_ = dict.get<scalar>("inertia")*symmTensor::I;
    }
    else if (dict.found("inertiaTensor"))
    {
        inertia_ = dict.get<symmTensor>("inertiaTensor");
    }
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
    // The initial angular velocity about the free axes only
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


Foam::tensor Foam::freeMotion::globalInertia(const quaternion& q) const
{
    const tensor R(q.R());
    return density_*(R & tensor(inertia_) & R.T());
}


void Foam::freeMotion::setMassProperties
(
    const scalar volume,
    const symmTensor& inertia
)
{
    if (volume_ <= 0)
    {
        volume_ = volume;
    }
    if (!hasInertia_ && volume > 0)
    {
        // the body's radii of gyration with the volume in use
        inertia_ = inertia/volume*volume_;
    }
    if (volume_ <= 0)
    {
        FatalIOErrorInFunction(dict_)
            << "the body gives no volume: set volume (and inertia if it "
            << "rotates) in the motion dictionary" << exit(FatalIOError);
    }
    if (rotates())
    {
        // Positive about every free axis, at any orientation: positive
        // definite
        const vector ev(eigenValues(inertia_));
        if (cmptMin(ev) <= 0)
        {
            FatalIOErrorInFunction(dict_)
                << "the body gives no moment of inertia, or one that is not "
                << "positive definite (" << inertia_ << "): set inertia or "
                << "inertiaTensor in the motion dictionary"
                << exit(FatalIOError);
        }
    }
}


void Foam::freeMotion::readState(const dictionary& dict)
{
    displacement_ = dict.get<vector>("displacement");
    velocity_ = dict.get<vector>("velocity");
    omega_ = dict.get<vector>("angularVelocity");
    if (dict.found("orientation"))
    {
        orientation_ = dict.get<quaternion>("orientation");
    }
    else
    {
        // A state of the single-axis implementation
        const scalar angle = dict.get<scalar>("angle");
        if (cmptSum(rotationMask_) != 1 && angle != 0)
        {
            FatalIOErrorInFunction(dict)
                << "a state with an angle needs exactly one free rotation "
                << "axis" << exit(FatalIOError);
        }
        orientation_ =
            (angle != 0) ? quaternion(rotationMask_, angle) : quaternion::I;
    }
    displacement0_ = displacement_;
    velocity0_ = velocity_;
    orientation0_ = orientation_;
    omega0_ = omega_;
}


Foam::dictionary Foam::freeMotion::state() const
{
    dictionary d;
    d.add("displacement", displacement_);
    d.add("velocity", velocity_);
    d.add("orientation", orientation_);
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

    // Inertia of the body in the global frame at the old orientation and at
    // the current iterate of the new one, and the unscaled one (fluid
    // density) for the bound of the added inertia
    const tensor I0(globalInertia(orientation0_));
    const tensor I1(globalInertia(orientation_));
    const tensor R1(orientation_.R());
    const tensor J1(R1 & tensor(inertia_) & R1.T());

    if (implicit())
    {
        // Bounds of the effective added mass: the given coefficient, and a
        // cap against a secant from two nearly equal iterates
        const scalar maMin = addedMass_*volume_;
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
        const scalar k2 = tr(inertia_)/3/volume_;
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
        // The same per free rotational direction (global axes), with the
        // diagonal of the inertia
        forAll(rotationMask_, i)
        {
            if (rotationMask_[i] == 0) continue;
            const scalar IaMin = addedMass_*J1(i, i);
            scalar Ia = max(DECAY*Ia_[i], IaMin);
            const scalar dw = Ok[i] - Oprev_[i];
            if (hasPrevious_ && mag(dw) > SECANT*max(mag(Ok), VSMALL))
            {
                const scalar est = -dt*(M[i] - Mprev_[i])/dw;
                Ia = max(Ia, min(est, CAP*(I1(i, i) + IaMin)));
            }
            Ia_[i] = Ia;
        }
        Uprev_ = Uk;
        Fprev_ = F;
        Oprev_ = Ok;
        Mprev_ = M;
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

    // Rotation about the free axes, the same scheme for the angular
    // momentum: (I^{n+1} + I_a) omega^{k+1} = I^n omega^n + I_a omega^k + dt M
    // with the components about the fixed axes zero (the system restricted
    // to the free ones). I^{n+1} = I(q^{n+1}) depends on the new orientation
    // q^{n+1} = q(omega^{n+1} dt) q^n: the two are iterated to consistency
    // (one pass for an isotropic body), so that I^n omega^n of the next step
    // is the angular momentum and a torque-free body conserves it (the
    // gyroscopic term is the change of I with the orientation)
    if (rotates())
    {
        const vector Ia(implicit() ? Ia_ : vector::zero);
        vector b((I0 & omega0_) + dt*M);
        for (direction i = 0; i < 3; ++i)
        {
            b[i] = (rotationMask_[i] > 0) ? b[i] + Ia[i]*Ok[i] : 0;
        }

        static const label NITER = 100;
        static const scalar TOL = 1e-13;
        quaternion q(orientation_);
        vector w(Ok);
        scalar change = GREAT;
        // Every pass solves with the inertia at the orientation of the
        // previous one; converged when the angular velocity repeats
        for (label iter = 0; iter < NITER && change > 0; ++iter)
        {
            tensor A(globalInertia(q));
            label nFree = 0;
            direction iFree = 0;
            for (direction i = 0; i < 3; ++i)
            {
                A(i, i) += Ia[i];
                if (rotationMask_[i] > 0)
                {
                    ++nFree;
                    iFree = i;
                    continue;
                }
                for (direction j = 0; j < 3; ++j)
                {
                    A(i, j) = (i == j ? 1 : 0);
                    A(j, i) = (i == j ? 1 : 0);
                }
            }
            vector wNew(Zero);
            if (nFree == 1)
            {
                // One free axis: a scalar equation
                wNew[iFree] = b[iFree]/A(iFree, iFree);
            }
            else
            {
                wNew = inv(A) & b;
            }
            const scalar dw = mag(wNew - w);
            w = wNew;
            const scalar mw = mag(w);
            q =
                (mw > 0)
              ? normalised(quaternion(w/mw, dt*mw)*orientation0_)
              : orientation0_;
            // (the first pass solved with the orientation of the iterate,
            // not the one of this omega)
            change = (iter > 0 && dw <= TOL*max(mw, VSMALL)) ? 0 : max(dw, SMALL);
        }
        if (change > 0)
        {
            WarningInFunction
                << "angular velocity and orientation not consistent after "
                << NITER << " iterations (change " << change << ")" << endl;
        }
        omega_ = w;
        orientation_ = q;
    }

    // Change of the iterate against the velocity scale of the body
    const scalar k =
        (volume_ > 0 ? Foam::sqrt(max(tr(inertia_)/3, 0)/volume_) : 0);
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
    orientation0_ = orientation_;
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
