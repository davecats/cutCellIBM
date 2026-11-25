/*---------------------------------------------------------------------------*\
    cutCellIBM: cut-cell immersed boundary method for OpenFOAM

Class
    Foam::bodyMotion

Description
    Rigid motion of an immersed body: a translation d(t) with velocity
    u(t) = d'(t), and optionally a rotation with constant angular velocity
    omega about the body's own centre. The velocity must be the exact
    derivative of the displacement: the first supplies the wall condition,
    the second (through the wall flux) the mass source, and if they disagree
    the two halves of the boundary condition pull against each other.

    Selected with the "motion" sub-dictionary of a body:
    \verbatim
    motion { type none; }                                   // default
    motion { type translating;  velocity (0.37 0 0); }
    motion { type oscillating;  amplitude (1 0 0); frequency 0.05; phase 0; }
    motion { type prescribed;   displacement table (...); velocity table (...); }
    motion { type external; }   // state set from code, see setState()
    motion { type free; density 2; ... }    // driven by the flow, see freeMotion
    \endverbatim
    Any type accepts "angularVelocity (0 0 w)" (rad/s, about the body centre)
    and "centreOfRotation (x y z)", the initial position of the point the
    body rotates about and is translated with (default: the body centre).
    The level set is periodic about this point.

    "invariant yes" (any type, default no) declares that the motion leaves
    the body unchanged as a set: a cylinder or a body of revolution turning
    about its axis, a sphere turning about its centre, a cylinder sliding
    along its axis, a plane sliding in itself. The geometry is then built
    once, with the body as defined, and the motion only supplies the wall
    velocity u_b = u + omega ^ (x - c) of the no-slip wall; there is no wall
    flux, and a steady solver can run it. cutCellGeometry checks the
    invariance at start-up and at every update (fatal otherwise).

    The "external" type is the hook for bodies driven by the solver itself,
    e.g. a point particle integrated with the flow: the driver calls
    setState(displacement, velocity) before cutCellGeometry::update().

\*---------------------------------------------------------------------------*/

#ifndef bodyMotion_H
#define bodyMotion_H

#include "dictionary.H"
#include "point.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "Function1.H"
#include "quaternion.H"

namespace Foam
{

class bodyMotion
{
protected:

    const dictionary dict_;
    vector omega_;

    //- The motion leaves the body unchanged as a set
    bool invariant_;

    //- Initial centre of rotation, if given
    bool hasCentreOfRotation_;
    point centreOfRotation_;

public:

    TypeName("bodyMotion");

    declareRunTimeSelectionTable
    (
        autoPtr,
        bodyMotion,
        dictionary,
        (const dictionary& dict),
        (dict)
    );

    explicit bodyMotion(const dictionary& dict);

    //- Select from the "motion" sub-dictionary of a body (none if absent)
    static autoPtr<bodyMotion> New(const dictionary& bodyDict);

    virtual ~bodyMotion() = default;

    //- True if the body moves at all
    virtual bool moving() const { return true; }

    //- True if the motion is declared to leave the body unchanged
    bool invariant() const { return invariant_; }

    //- Displacement of the body centre from its initial position
    virtual vector displacement(const scalar t) const = 0;

    //- Velocity of the body centre
    virtual vector velocity(const scalar t) const = 0;

    //- Angular velocity about the body centre
    const vector& omega() const { return omega_; }

    //- True if the centre of rotation is given (else the body centre)
    bool hasCentreOfRotation() const { return hasCentreOfRotation_; }

    //- Initial centre of rotation
    const point& centreOfRotation() const { return centreOfRotation_; }

    //- True if the body turns (the geometry is rotated by angle() about
    //  axis())
    virtual bool rotates() const { return mag(omega_) > 0; }

    //- Unit axis the rotation angle is measured about (omega/|omega|)
    virtual vector axis() const { return omega_/max(mag(omega_), VSMALL); }

    //- Rotation angle about axis() at time t
    virtual scalar angle(const scalar t) const { return mag(omega_)*t; }

    //- Orientation at time t: the rotation from the body at rest to the
    //  body at t, about the reference point (identity if the body does not
    //  turn). Default: angle(t) about axis(); a free body integrates it.
    virtual quaternion orientation(const scalar t) const
    {
        return rotates() ? quaternion(axis(), angle(t)) : quaternion::I;
    }

    //- Set the state from outside (external motion only)
    virtual void setState(const vector& displacement, const vector& velocity);
};


class noMotion : public bodyMotion
{
public:
    TypeName("none");
    explicit noMotion(const dictionary& dict);
    virtual bool moving() const;
    virtual vector displacement(const scalar t) const;
    virtual vector velocity(const scalar t) const;
};


class translatingMotion : public bodyMotion
{
    vector U_;
public:
    TypeName("translating");
    explicit translatingMotion(const dictionary& dict);
    virtual vector displacement(const scalar t) const;
    virtual vector velocity(const scalar t) const;
};


class oscillatingMotion : public bodyMotion
{
    vector amplitude_;
    scalar frequency_;
    scalar phase_;
public:
    TypeName("oscillating");
    explicit oscillatingMotion(const dictionary& dict);
    virtual vector displacement(const scalar t) const;
    virtual vector velocity(const scalar t) const;
};


class prescribedMotion : public bodyMotion
{
    autoPtr<Function1<vector>> displacement_;
    autoPtr<Function1<vector>> velocity_;
public:
    TypeName("prescribed");
    explicit prescribedMotion(const dictionary& dict);
    virtual vector displacement(const scalar t) const;
    virtual vector velocity(const scalar t) const;
};


class externalMotion : public bodyMotion
{
    vector displacement_;
    vector velocity_;
public:
    TypeName("external");
    explicit externalMotion(const dictionary& dict);
    virtual vector displacement(const scalar t) const;
    virtual vector velocity(const scalar t) const;
    virtual void setState(const vector& displacement, const vector& velocity);
};

} // End namespace Foam

#endif
