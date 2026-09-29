/*---------------------------------------------------------------------------*\
    cutCellIBM: cut-cell immersed boundary method for OpenFOAM

Class
    Foam::ibmBody

Description
    An immersed body. Two kinds exist:

    - level-set bodies: the body is the region where a scalar level set
      phi(x) is <= 0 (phi > 0 in the fluid). Several level-set bodies are
      combined by taking the minimum of their level sets, so their union is
      exact. The level set need not be a distance function; only its sign and
      the location of its zero level are used (bisection along mesh edges).

    - alpha-field bodies (see alphaFieldBody): the fluid volume fraction is
      read from a field and the face fractions are reconstructed from it.

    Each body carries a wall condition, noSlip (default) or freeSlip, which
    selects whether the viscous traction of the wall segment is added to the
    momentum diagonal (see cutCellGeometry).

    A body may move rigidly (see bodyMotion). The level set of a moving
    body is evaluated by transforming the query points into the body frame:
    minimum-image wrap in the periodic directions about the current centre,
    then the inverse translation and rotation. Static bodies skip all of it.

    Dictionary entries common to all bodies:
    \verbatim
    <name>
    {
        type        cylinder | sphere | box | plane | <searchableSurface type>
                    | alphaField;
        wallType    noSlip;      // or freeSlip (optional)
        invert      false;       // true: fluid *inside* the surface (optional)
        motion      { type translating; velocity (0.5 0 0); }   // optional
        ...type-specific entries
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef ibmBody_H
#define ibmBody_H

#include "dictionary.H"
#include "pointField.H"
#include "scalarField.H"
#include "runTimeSelectionTables.H"
#include "Enum.H"
#include "autoPtr.H"
#include "bodyMotion.H"

namespace Foam
{

class fvMesh;

class ibmBody
{
public:

    enum wallType
    {
        NOSLIP,
        FREESLIP
    };

    static const Enum<wallType> wallTypeNames;


protected:

    const word name_;
    const fvMesh& mesh_;
    const dictionary dict_;
    wallType wallType_;
    bool invert_;

    //- Rigid motion (none by default)
    autoPtr<bodyMotion> motion_;

    //- Periodic lengths of the domain (0: not periodic), set by the geometry
    vector periodicLengths_;

    // Current motion state, set by moveTo()
    scalar t_;
    vector displacement_;
    vector velocity_;
    scalar angle_;
    vector axis_;


public:

    TypeName("ibmBody");

    declareRunTimeSelectionTable
    (
        autoPtr,
        ibmBody,
        dictionary,
        (const word& name, const fvMesh& mesh, const dictionary& dict),
        (name, mesh, dict)
    );

    ibmBody(const word& name, const fvMesh& mesh, const dictionary& dict);

    static autoPtr<ibmBody> New
    (
        const word& name,
        const fvMesh& mesh,
        const dictionary& dict
    );

    virtual ~ibmBody() = default;


    // Member functions

        const word& name() const { return name_; }

        wallType wall() const { return wallType_; }

        const dictionary& dict() const { return dict_; }

        //- True if the body is described by a level set
        virtual bool hasLevelSet() const { return true; }

        //- Level set (before inversion): positive in the fluid, negative in
        //  the body. Implemented by level-set bodies.
        virtual void levelSet(const pointField& pts, scalarField& phi) const;

        //- Level set including the invert switch and the body motion:
        //  positive in the fluid
        void evaluate(const pointField& pts, scalarField& phi) const;

        //- Reference point of the body at rest (centre of rotation, and
        //  centre of the periodic image); as given by the body definition
        virtual point centre() const = 0;


    // Motion

        //- True if the body moves
        virtual bool moving() const { return motion_->moving(); }

        bodyMotion& motion() { return motion_(); }
        const bodyMotion& motion() const { return motion_(); }

        //- Set the periodic lengths of the domain (0: not periodic)
        void setPeriodicLengths(const vector& L) { periodicLengths_ = L; }

        //- Update the motion state to time t
        virtual void moveTo(const scalar t);

        //- Initial reference point of the motion: the centre of rotation
        //  if the motion gives one, otherwise the body centre
        point referencePoint() const
        {
            return motion_->hasCentreOfRotation()
                 ? motion_->centreOfRotation() : centre();
        }

        //- Current centre (of rotation and of the periodic image), initial
        //  reference point plus displacement
        point currentCentre() const { return referencePoint() + displacement_; }

        const vector& displacement() const { return displacement_; }

        //- Translational velocity of the body
        const vector& velocityOfCentre() const { return velocity_; }

        //- Angular velocity about the current centre
        const vector& omega() const { return motion_->omega(); }

        //- True if the body rotates: a level-set body with a non-zero
        //  angular velocity (the wall flux then needs the rotational part,
        //  see cutCellGeometry::calcMovingTerms)
        bool rotating() const
        {
            return hasLevelSet() && moving() && mag(omega()) > 0;
        }

        //- Position relative to the current centre, x - c, with the minimum
        //  image in the periodic directions
        vector relativePosition(const point& x) const;

        //- Body velocity at the given points (rigid: U + omega x r)
        virtual void velocity(const pointField& pts, vectorField& U) const;

        //- Query points transformed into the body frame at rest
        void toBodyFrame(const pointField& pts, pointField& out) const;

        //- True if the points must be transformed before evaluation
        bool needsTransform() const;
};

} // End namespace Foam

#endif
