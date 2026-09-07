/*---------------------------------------------------------------------------*\
    cutCellIBM: cut-cell immersed boundary method for OpenFOAM

Class
    Foam::cutCellGeometry

Description
    Cut-cell immersed-boundary geometry and operators, after the method of
    .asset/fvm/SIMPLEC_faceArea_noSlip.ipynb (see doc/freeslip_ibm_part2.pdf).

    Every face carries its fluid area fraction theta_f in [0,1], every cell
    its fluid volume fraction alpha_P in [0,1]. The fluid part of a cut cell is
    a polyhedron bounded by the wet parts of its faces and by a wall segment,
    and geometric closure

        sum_f theta_f Sf + Sw = 0

    defines the wall area vector Sw (pointing from the fluid into the body).
    Nothing else is needed: with Sw defined this way the wall pressure force is
    absorbed into a theta-weighted Gauss gradient, no-penetration follows from
    theta-weighted continuity, and the no-slip wall traction is a single
    diagonal term  nu A_wall/d_wall  with  d_wall = alpha V/(2 A_wall).

    All fields are computed once at construction from constant/ibmProperties:
    \verbatim
    alphaMin      0.02;   // cells with less fluid are absorbed into the body
    nSubSamples   32;     // per direction, for the cell fluid fraction
    nBisect       40;     // bisection steps locating the wall on cut edges
    writeFields   true;   // write ibmAlpha, ibmTheta, ... at start-up
    writeForces   true;   // postProcessing/ibmForces/<time>/forces.dat

    bodies
    {
        cylinder { type cylinder; centre (5 5 0); axis (0 0 1); radius 1; }
        hull     { type triSurfaceMesh; file "hull.stl"; wallType noSlip; }
        given    { type alphaField; field alpha.body; }
    }
    \endverbatim

    Registered fields (on the mesh registry, available to fvOptions and
    function objects): ibmAlpha, ibmTheta, ibmSw, ibmWallDistance,
    ibmNoSlipCoeff, ibmBlankCoeff, ibmBodyIndex.

    The operators are used by the solvers as
    \verbatim
        fvm::ddt(ibm.alpha(), U)
      + fvm::div(phi, U)                       // phi = ibm.flux(U)
      - fvm::laplacian(ibm.theta()*nuf, U)
      + fvm::Sp(nu*ibm.noSlipCoeff(), U)       // wall traction
      + fvm::Sp(nu*ibm.blankCoeff(), U)        // u = 0 in solid cells
     == -ibm.grad(p)
    \endverbatim
    and, in the pressure equation, the laplacian coefficient is
    ibm.theta()*interpolate(rAtU) and the solid rows are made regular with
    ibm.constrainPressure(pEqn).

\*---------------------------------------------------------------------------*/

#ifndef cutCellGeometry_H
#define cutCellGeometry_H

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"
#include "bitSet.H"
#include "OFstream.H"
#include "ibmBody.H"

namespace Foam
{

class alphaFieldBody;

class cutCellGeometry
:
    public IOdictionary
{
    // Private data

        const fvMesh& mesh_;

        PtrList<ibmBody> bodies_;

        scalar alphaMin_;
        label nSubSamples_;
        label nBisect_;
        bool writeFields_;
        bool writeForces_;

        //- Fluid volume fraction
        volScalarField alpha_;

        //- Fluid area fraction of every face
        surfaceScalarField theta_;

        //- Wall area vector of every cell, Sw = -sum_f theta_f Sf
        volVectorField Sw_;

        //- Wall area |Sw|
        volScalarField Awall_;

        //- Wall-normal distance of the cell average, alpha V/(2 Awall)
        volScalarField dWall_;

        //- Awall/(dWall V) on no-slip wall cells, 0 elsewhere [1/m^2]
        volScalarField noSlipCoeff_;

        //- sum_f |Sf| deltaCoeff_f / V in solid cells, 0 elsewhere [1/m^2]
        volScalarField blankCoeff_;

        //- Index of the body owning each cut/solid/wall cell (-1: fluid)
        volScalarField bodyIndex_;

        labelList cellBody_;

        bitSet solidCells_;
        bitSet cutCells_;
        bitSet wallCells_;

        scalar fluidVolume_;

        mutable autoPtr<OFstream> forcesFilePtr_;


    // Private member functions

        //- Minimum over the level sets of all level-set bodies
        void unionLevelSet(const pointField& pts, scalarField& phi) const;

        //- Fluid fraction of every face from the level set at the mesh points
        //  and the crossing point on every cut edge
        void faceFractions
        (
            const scalarField& pointPhi,
            const pointField& edgeCross,
            scalarField& theta
        ) const;

        //- Fluid fractions from the level-set bodies
        void calcFromLevelSet
        (
            scalarField& alpha,
            scalarField& theta,
            labelList& cellBody
        ) const;

        //- Fluid fractions from an alpha-field body
        void calcFromAlphaField
        (
            const alphaFieldBody& body,
            scalarField& alpha,
            scalarField& theta
        ) const;

        //- Fluid fraction of one cell by sub-sampling the level set
        scalar sampleCellFraction(const label celli) const;

        //- Everything that follows from alpha and theta
        void calcDerived();

        void report() const;

        //- No copy construct
        cutCellGeometry(const cutCellGeometry&) = delete;
        void operator=(const cutCellGeometry&) = delete;


public:

    TypeName("cutCellGeometry");

    // Constructors

        //- Construct from mesh, reading constant/ibmProperties
        explicit cutCellGeometry(const fvMesh& mesh);


    //- Destructor
    virtual ~cutCellGeometry() = default;


    // Member Functions

        // Access

            const fvMesh& mesh() const { return mesh_; }
            const PtrList<ibmBody>& bodies() const { return bodies_; }
            const volScalarField& alpha() const { return alpha_; }
            const surfaceScalarField& theta() const { return theta_; }
            const volVectorField& Sw() const { return Sw_; }
            const volScalarField& Awall() const { return Awall_; }
            const volScalarField& dWall() const { return dWall_; }
            const volScalarField& noSlipCoeff() const { return noSlipCoeff_; }
            const volScalarField& blankCoeff() const { return blankCoeff_; }
            const labelList& cellBody() const { return cellBody_; }
            const bitSet& solidCells() const { return solidCells_; }
            const bitSet& cutCells() const { return cutCells_; }
            const bitSet& wallCells() const { return wallCells_; }
            scalar alphaMin() const { return alphaMin_; }

            //- Total fluid volume, sum(alpha V), global
            scalar fluidVolume() const { return fluidVolume_; }


        // Operators

            //- Cut-cell Gauss gradient, sum_f theta_f (p_f - p_P) Sf / V.
            //  Exactly zero for a uniform field; includes the wall force.
            tmp<volVectorField> grad(const volScalarField& p) const;

            //- Volume flux through the wet part of every face
            tmp<surfaceScalarField> flux(const volVectorField& U) const;

            //- theta_f * interpolate(f), the diffusivity of a cut laplacian
            tmp<surfaceScalarField> thetaInterpolate
            (
                const volScalarField& f
            ) const;

            tmp<surfaceScalarField> thetaInterpolate
            (
                const tmp<volScalarField>& tf
            ) const;

            //- Transient Rhie-Chow correction (Euler), theta-consistent
            //  replacement of fvc::ddtCorr(U, phi)
            tmp<surfaceScalarField> ddtCorr
            (
                const volVectorField& U,
                const surfaceScalarField& phi
            ) const;

            //- Give the (empty) rows of solid cells a regular diagonal
            void constrainPressure(fvScalarMatrix& pEqn) const;

            //- Fatal error if the pressure reference cell is not in the fluid
            void checkRefCell(const label pRefCell) const;

            //- Zero a field in the solid cells
            template<class Type>
            void zeroSolid
            (
                GeometricField<Type, fvPatchField, volMesh>& f
            ) const;

            //- Force on every body (per unit density): pressure and viscous
            void forces
            (
                const volScalarField& p,
                const volVectorField& U,
                const volScalarField& nu,
                List<vector>& Fp,
                List<vector>& Fv
            ) const;

            //- Print the forces and append them to the forces file
            void reportForces
            (
                const volScalarField& p,
                const volVectorField& U,
                const volScalarField& nu
            ) const;


        // Write

            //- Write the geometry fields to the current time directory
            void writeGeometryFields() const;

            //- Required by regIOobject
            virtual bool writeData(Ostream&) const { return true; }
};

} // End namespace Foam

#ifdef NoRepository
    #include "cutCellGeometryTemplates.C"
#endif

#endif
