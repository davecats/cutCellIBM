/*---------------------------------------------------------------------------*\
    cutCellIBM: cut-cell immersed boundary method for OpenFOAM

Class
    Foam::freeMotion

Description
    Rigid body moving freely under the force and moment of the flow on it
    (Newton-Euler), integrated by the solver from the IBM wall forces
    (cutCellGeometry::moveFreeBodies). All quantities are per unit fluid
    density, as the forces of cutCellGeometry: the body's mass is
    density*volume with density = rho_s/rho.

    Forces on the body (per unit fluid density):
        F = F_p + F_v                         (IBM wall pressure and viscous)
          + volume*G                          (mean pressure gradient of
                                               ibmMeanVelocityForce, if any:
                                               the periodic pressure of the
                                               solver excludes it)
          + density*volume*g                  (gravity, optional)
          + externalForce                     (optional, constant)
        M = M_p + M_v                         (about the current centre)

    Translation and rotation are restricted to the directions of
    "translation" and "rotation" (components 0 or 1 in the global axes);
    at most one rotational direction may be free, so that the rotation is
    about a fixed axis and its angle is a scalar. The body is moved with
    d^{n+1} = d^n + dt U^{n+1} (implicit Euler for the position), which is
    the displacement the wall flux U.Sw of cutCellGeometry assumes.

    Coupling (the velocity at the new time):
    - explicit: from the force of the previous time step,
          m (U^{n+1} - U^n)/dt = F^n;
      stable only if the body is heavy enough against its added mass
      (rho_s/rho well above the added-mass coefficient, 0.5 for a sphere in
      unbounded fluid, more when confined).
    - implicit: the body velocity is corrected in every PIMPLE outer
      corrector with the force of the flow just solved, with an effective
      added mass m_a per free direction (fluid density 1):
          (m + m_a) U^{n+1,k+1} = m U^n + m_a U^{n+1,k} + dt F^{n+1,k}.
      At convergence (U^{n+1,k+1} = U^{n+1,k}) this is the implicit Euler
      step m (U^{n+1} - U^n)/dt = F^{n+1}, whatever m_a. It is a Newton
      step for m (U - U^n) - dt F(U) = 0 with m_a = -dt dF/dU, which holds
      the added mass of the fluid and the drag (dt times the drag
      coefficient; dominant in creeping flow). m_a is the secant estimate
      -dt dF/dU from the last two iterates of the step, bounded below by
      addedMass*volume and kept from step to step: it is raised at once by a
      larger secant and decays by a factor 0.8 per update otherwise (noise
      of the unconverged flow and jumps of the loads where cells change
      status must not pull it below the true value, which makes the
      iteration overshoot, nor hold it far above, which freezes the body). The error contracts by
      |m_a - M|/(m + m_a) per corrector (M the true -dt dF/dU), < 1 for any
      m_a >= M, so the iteration converges for any rho_s, also rho_s/rho < 1
      and down to zero (a bubble). The rotation uses the same scheme with
      the inertia (bound addedMass*inertia). The corrector of the first
      outer iteration is the predictor, from the loads at the old time
      level.

    Mass properties: "volume" and "inertia" (per unit density: the volume
    and int r_perp^2 dV about the free rotation axis through the centre) of
    the part of the body inside the domain; for a sphere the defaults are
    4/3 pi R^3 and 2/5 R^2 volume. With a symmetry plane through the body
    give the values of the part in the domain (half a sphere: volume
    2/3 pi R^3; the inertia default follows the volume).

    Usage (a body of constant/ibmProperties):
    \verbatim
    motion
    {
        type            free;
        density         2;          // rho_s/rho
        translation     (1 0 0);    // free directions (default (1 1 1))
        rotation        (0 0 1);    // free axis (default (0 0 0))
        velocity        (1.5 0 0);  // initial (default 0)
        angularVelocity (0 0 1);    // initial (default 0)
        coupling        implicit;   // or explicit (default implicit)
        addedMass       0.5;        // lower bound of m_a/volume (implicit;
                                    // default 0.5)
        volume          0.0335;     // optional (default: from the body)
        inertia         0.00134;    // optional (default: from the body)
        gravity         (0 0 0);    // optional
        externalForce   (0 0 0);    // optional, per unit fluid density
        meanPressureGradient yes;   // optional (default yes)
    }
    \endverbatim
    The state (displacement, velocity, angle, angular velocity) is written
    to <time>/uniform/ibmFreeBodies and read on restart.

\*---------------------------------------------------------------------------*/

#ifndef freeMotion_H
#define freeMotion_H

#include "bodyMotion.H"
#include "Enum.H"

namespace Foam
{

class freeMotion
:
    public bodyMotion
{
public:

    enum couplingType { EXPLICIT, IMPLICIT };
    static const Enum<couplingType> couplingNames;

private:

    scalar density_;
    vector translationMask_;
    vector rotationMask_;
    couplingType coupling_;
    scalar addedMass_;
    scalar volume_;
    scalar inertia_;
    vector gravity_;
    vector externalForce_;
    bool meanPressureGradient_;

    // State at the new time level (current iterate) and at the old one
    vector displacement_;
    vector velocity_;
    scalar theta_;           // rotation angle about axis()
    vector displacement0_;
    vector velocity0_;
    vector omega0_;
    scalar theta0_;

    //- Change of the velocities in the last correction, relative to the
    //  largest velocity scale (for the report)
    scalar residual_;

    // Implicit coupling: effective added mass per translational direction
    // and added inertia (secant estimates), and the previous iterate
    vector ma_;
    scalar Ia_;
    bool hasPrevious_;
    vector Uprev_;
    vector Fprev_;
    scalar Oprev_;
    scalar Mprev_;

public:

    TypeName("free");

    explicit freeMotion(const dictionary& dict);

    virtual vector displacement(const scalar t) const;
    virtual vector velocity(const scalar t) const;

    virtual bool rotates() const { return cmptMax(rotationMask_) > 0; }
    virtual vector axis() const;
    virtual scalar angle(const scalar t) const { return theta_; }

    // Access

        scalar density() const { return density_; }
        scalar addedMass() const { return addedMass_; }
        couplingType coupling() const { return coupling_; }
        bool implicit() const { return coupling_ == IMPLICIT; }
        bool meanPressureGradient() const { return meanPressureGradient_; }
        scalar volume() const { return volume_; }
        scalar inertia() const { return inertia_; }
        scalar mass() const { return density_*volume_; }
        scalar residual() const { return residual_; }

        //- Effective added mass per direction and added inertia in use
        //  (implicit coupling)
        const vector& effectiveAddedMass() const { return ma_; }
        scalar effectiveAddedInertia() const { return Ia_; }

        //- Body force on the body per unit fluid density (gravity and the
        //  constant external force)
        vector bodyForce() const
        {
            return density_*volume_*gravity_ + externalForce_;
        }

    // Set-up

        //- Mass properties from the body if not given (volume <= 0: none)
        void setMassProperties(const scalar volume, const scalar inertia);

        //- State from a restart dictionary
        void readState(const dictionary& dict);

        //- State as a dictionary
        dictionary state() const;

    // Time integration

        //- Start a time step of size dt: keep the old state and predict the
        //  new one from force F and moment M at the old time level
        void predict(const scalar dt, const vector& F, const vector& M);

        //- Correct the new state from force F and moment M of the current
        //  iterate (implicit coupling)
        void correct(const scalar dt, const vector& F, const vector& M);

private:

        //- Velocity update (implicit: with the secant added mass)
        void advance(const scalar dt, const vector& F, const vector& M);
};

} // End namespace Foam

#endif
