/*---------------------------------------------------------------------------*\
    cutCellIBM: cut-cell immersed boundary method for OpenFOAM

Class
    Foam::fv::ibmMeanVelocityForce

Description
    Constant flow rate forcing for the cut-cell immersed boundary, after
    meanVelocityForce: a uniform body force along flowDir is adjusted so
    that the bulk velocity of the fluid sits on Ubar.

    Differences from meanVelocityForce:
    - the bulk velocity is averaged over the fluid volume alpha V, and the
      force acts on the fluid volume of each cell only, so that solid and
      cut cells receive no spurious momentum (alpha is the field ibmAlpha
      registered by cutCellGeometry);
    - the force increment is computed with the mobility of a uniform
      velocity change, 1/(A - H1) (the SIMPLEC coefficient), instead of 1/A.
      For a uniform force the neighbour couplings of the momentum matrix
      cancel, so this is the exact bulk response; 1/A underestimates it by
      the diffusion number nu dt/h^2 with implicit viscosity, which makes
      the standard controller overshoot;
    - the increments of the predictor and of every corrector call are
      accumulated within a time step instead of overwritten;
    - a zero target is allowed ("Ubar (0 0 0); flowDir (1 0 0);"), the
      natural choice for a body moving through an otherwise quiescent fluid.
    - with "includeMovingBodies yes" Ubar is the mean velocity of the fluid
      and of the moving bodies together (the total flux through a section
      of a suspension, J for a body in a pipe): the fluid target is
      (|Ubar| (Vf + sum V_b) - sum V_b U_b.flowDir)/Vf, V_b the discrete
      solid volume of a moving body in the domain and U_b the velocity of
      its centre (selectionMode all; bodies with a geometry that moves);
    - the current force per unit mass, flowDir*gradP, is registered as the
      uniformDimensionedVectorField ibmMeanPressureGradient: the periodic
      pressure of the solver leaves out this mean gradient, and freely
      moving bodies (freeMotion) add it on their volume.

    Usage, in constant/fvOptions:
    \verbatim
    momentumSource
    {
        type            ibmMeanVelocityForce;
        selectionMode   all;
        fields          (U);
        Ubar            (0.5 0 0);
        flowDir         (1 0 0);    // required when Ubar is zero
        relaxation      1.0;
        includeMovingBodies no;     // optional (default no)
    }
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef ibmMeanVelocityForce_H
#define ibmMeanVelocityForce_H

#include "cellSetOption.H"
#include "fvOptionList.H"

namespace Foam
{
namespace fv
{

class ibmMeanVelocityForce
:
    public cellSetOption
{
    // Private data

        vector Ubar_;
        vector flowDir_;
        scalar relaxation_;

        //- Accumulated force, and increment of the current step
        scalar gradP0_;
        scalar dGradP_;

        //- Mobility of a uniform velocity change, 1/(A - H1)
        autoPtr<volScalarField> rAPtr_;

        //- Fluid volume of the selected cells
        scalar Vfluid_;

        //- True for a steady-state time scheme (no alpha/deltaT floor)
        bool steady_;

        //- Ubar is the mean over the fluid and the moving bodies
        bool includeMovingBodies_;


    // Private member functions

        const volScalarField& alpha() const;

        scalar magUbarAve(const volVectorField& U) const;

        void writeProps(const scalar gradP) const;

        //- Target of the fluid-volume average along flowDir
        scalar targetUbar() const;

        //- Set the registered ibmMeanPressureGradient to the current force
        void publishGradient() const;


public:

    TypeName("ibmMeanVelocityForce");

    ibmMeanVelocityForce
    (
        const word& sourceName,
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh
    );

    virtual ~ibmMeanVelocityForce() = default;


    // Member Functions

        //- Current force per unit density
        scalar gradP() const { return gradP0_ + dGradP_; }

        virtual void correct(volVectorField& U);

        //- Add the gradient increments of the current outer corrector to H
        //  (the explicit part of the momentum equation, per unit volume) as
        //  the source alpha flowDir dGradP. correct() applies them to U as
        //  the response s = alpha dGradP/(A - H1); they enter the momentum
        //  source only at the next assembly (constrain). With them in H the
        //  pressure corrector keeps the shift exactly, HbyA(U + s) =
        //  HbyA(U) + s, and the increments of all correctors add up to the
        //  force the flow takes. Without them only the fraction H1/A of a
        //  shift survives a pressure corrector: nearly all of it in viscous
        //  flow (nu dt/h^2 >> 1), almost nothing in transient flow, where
        //  the summed increments overshoot (~ the number of correctors) and
        //  the outer iteration diverges (V6G Re 128, 2026-11-25)
        void addCorrectorSource(volVectorField& H) const;

        //- addCorrectorSource of every active ibmMeanVelocityForce of a list
        //  (the pressure correctors of ibmPimpleFoam and ibmSimpleFoam)
        static void addCorrectorSources
        (
            optionList& options,
            volVectorField& H
        );

        virtual void addSup(fvMatrix<vector>& eqn, const label fieldi);

        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        virtual void constrain(fvMatrix<vector>& eqn, const label fieldi);

        virtual bool read(const dictionary& dict);
};

} // End namespace fv
} // End namespace Foam

#endif
