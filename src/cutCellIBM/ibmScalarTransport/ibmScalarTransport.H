/*---------------------------------------------------------------------------*\
    cutCellIBM: cut-cell immersed boundary method for OpenFOAM

Class
    Foam::ibmScalarTransport

Description
    One passive scalar transported by the cut-cell flux, with the immersed
    wall seen either as a fixed-value wall or as a zero-flux wall, per body.
    Several scalars are held by ibmScalarTransportList, which reads
    constant/scalarTransportProperties.

    The equation is the cut-cell form of scalar transport,

        ddt(alpha, T) + div(phi, T) - laplacian(theta D, T)
      + Sp(D wallCoeff, T) + Sp(D blankCoeff, T)
     == D (wallCoeff + blankCoeff) T_wall + q alpha + fvOptions(T),

    where wallCoeff = Awall/(dWall V) on the wall cells of bodies seen as
    fixed-value walls (zero for zero-gradient walls, which need no term at
    all: the wall segment simply carries no flux), blankCoeff drives the
    solid interior to T_wall, q is a prescribed source per unit fluid volume
    in a circular region, and fvOptions(T) collects the sources and sinks of
    constant/fvOptions that list this field (semi-implicit, coded, ...);
    these act per cell volume, not per fluid volume, so they should be kept
    away from the body or made alpha-aware. Cells containing the wall are
    fluid cells: the scalar is transported into their fluid part through
    the wet fraction of their faces, and never across the wall segment.

    For a moving body no mass-source term appears: Reynolds transport over
    the fluid part of a cell gives exactly the form above (the relative flux
    through the wall is zero), and a uniform T is preserved because
    (alpha^{n+1} - alpha^n) V/dt = S and sum_f theta_f phi_f = -S cancel.

    A budget is reported at every step: fluid content sum(alpha V T), the
    injected amount (q and fvOptions), the amount that left through
    fixed-value walls, and the swept-volume defect a moving body carries,
    sum((alpha_prev - alpha^n) V T^n): the old fluid volume of the momentum
    and scalar equations is integrated from the wall flux and is only
    first-order consistent with the previous geometry, and a cell absorbed
    into the body hands over no content. The four must balance to
    round-off; for static bodies the defect is zero.

    Dictionary of one scalar (the field is read from the start time):
    \verbatim
    DT      0.1;                                        // diffusivity
    source  { centre (5 8.5 0); radius 0.5; rate 1; }   // q [T/s], optional
    walls                                               // per body, optional
    {
        cylinder { type fixedValue; value 0; }
        ".*"     { type zeroGradient; }
    }
    \endverbatim
    A body not matched in "walls" takes "scalarWall fixedValue|zeroGradient;
    scalarWallValue 0;" from its own entry in constant/ibmProperties, and
    zeroGradient if neither is given.

\*---------------------------------------------------------------------------*/

#ifndef ibmScalarTransport_H
#define ibmScalarTransport_H

#include "volFields.H"
#include "surfaceFields.H"
#include "fvOptions.H"
#include "cutCellGeometry.H"

namespace Foam
{

class ibmScalarTransport
{
    // Private data

        const cutCellGeometry& ibm_;
        const fvMesh& mesh_;
        const surfaceScalarField& phi_;
        fv::options& fvOptions_;

        volScalarField T_;
        dimensionedScalar DT_;

        //- Source per unit fluid volume, q
        volScalarField q_;

        //- Per body: fixed-value wall (true) or zero-flux wall, and value
        boolList bodyFixed_;
        scalarList bodyValue_;

        //- Awall/(dWall V) on fixed-value wall cells, 0 elsewhere [1/m^2]
        volScalarField wallCoeff_;

        //- Wall value of the body owning each cell (0 in the fluid)
        volScalarField Twall_;

        //- True for a steady-state time scheme (no budget in time)
        bool steady_;

        // Budget
        scalar content0_;
        scalar injected_;
        scalar wallFlux_;
        scalar sweptDefect_;
        scalarField alphaPrev_;


    // Private member functions

        //- Mask the wall coefficient and the wall value per body
        void updateWall();

        scalar content() const;


public:

    TypeName("ibmScalarTransport");

    //- Construct the scalar "name" from its dictionary
    ibmScalarTransport
    (
        const cutCellGeometry& ibm,
        const surfaceScalarField& phi,
        fv::options& fvOptions,
        const word& name,
        const dictionary& dict
    );

    virtual ~ibmScalarTransport() = default;


    // Member Functions

        const volScalarField& T() const { return T_; }

        const word& name() const { return T_.name(); }

        //- Solve one step (or one outer iteration) and report the budget
        void solve();
};

} // End namespace Foam

#endif
