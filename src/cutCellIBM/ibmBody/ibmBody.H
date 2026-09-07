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

    Dictionary entries common to all bodies:
    \verbatim
    <name>
    {
        type        cylinder | sphere | box | plane | <searchableSurface type>
                    | alphaField;
        wallType    noSlip;      // or freeSlip (optional)
        invert      false;       // true: fluid *inside* the surface (optional)
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

        //- Level set including the invert switch: positive in the fluid
        void evaluate(const pointField& pts, scalarField& phi) const;
};

} // End namespace Foam

#endif
