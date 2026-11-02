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
    ibmNoSlipCoeff, ibmBlankCoeff, ibmBodyIndex; for moving bodies also
    ibmUb (body velocity at the cell centres) and ibmMassSource.

    Moving bodies (doc/moving_ibm.pdf): the mesh never moves, the fractions
    do. update(), called at the start of every time step, moves the bodies
    to the current time, recomputes the geometry and forms
        S = U . Sw + S_rot            (wall flux, volume-integrated;
                                       S_rot the flux of omega x (x - c)
                                       through the wall, from the open
                                       wet faces),
        alpha^n V = alpha^{n+1} V - dt S   (swept volume, stored as
                                            alpha.oldTime()),
        u_b at the cell centres.
    Continuity then reads sum_f theta_f phi_f = -S, i.e. the pressure
    equation gains massSource() = S/V on its right-hand side; the no-slip
    wall traction and the blanking get the source nu c u_b through
    movingWallSource(); fvm::ddt(alpha, U) picks up alpha.oldTime()
    automatically. For static bodies update() returns at once and none of
    the moving-body fields exist.

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

        //- Awall/(dWall V) on every wall cell, 0 elsewhere [1/m^2]: the
        //  coefficient of a Dirichlet wall condition, before any masking
        volScalarField wallCoeff_;

        //- wallCoeff_ on no-slip wall cells, 0 elsewhere [1/m^2]
        volScalarField noSlipCoeff_;

        //- sum_f |Sf| deltaCoeff_f / V in solid cells, 0 elsewhere [1/m^2]
        volScalarField blankCoeff_;

        //- Face fractions and wet-face centroids of the geometric cut, over
        //  all mesh faces, before sliver absorption and coupled-face sync
        scalarField thetaGeo_;
        pointField wetCentreGeo_;

        //- Index of the body owning each cut/solid/wall cell (-1: fluid)
        volScalarField bodyIndex_;


    // Second-order wall treatment (secondOrder true), see doc/secondOrderStudy.md

        bool secondOrder_;

        //- Fluid centroid minus cell centre
        volVectorField centroidOffset_;

        //- Wet-face centroid minus face centre
        surfaceVectorField wetFaceOffset_;

        //- Centroid-to-centroid vector across each face (owner to neighbour)
        surfaceVectorField dCentroid_;

        //- 1/(n.d) on faces between two fluid cells, mesh value elsewhere
        surfaceScalarField cutDeltaCoeffs_;

        //- Union level set at the fluid centroids (level-set bodies)
        scalarField centroidLevelSet_;

        //- Offset of the wet-face centroid from the point the linear face
        //  interpolate belongs to, x'_f - m_f, m_f = w_f x_c,O + (1-w_f) x_c,N;
        //  zero on faces that need no correction
        surfaceVectorField faceValueOffset_;

        //- Level set of the geometric cut at the mesh points (positive in
        //  the fluid) and crossing point of every cut edge
        scalarField pointPhiGeo_;
        pointField edgeCrossGeo_;

        //- Wall pieces of every wall cell k (the wall polygon(s) of the
        //  geometric cut and the geometric wet parts of the faces closed by
        //  the post-processing): sum_k (x_k - x_c) S_k and sum_k S_k
        tensorField wallMoment_;
        vectorField wallPieceSum_;

        //- Work arrays filled by faceFractions/calcFromLevelSet
        mutable vectorField wetOffsetAll_;

        labelList cellBody_;

        //- Nearest level-set body of every cell (-1 if none)
        mutable labelList nearestBody_;

        bitSet solidCells_;
        bitSet cutCells_;
        bitSet wallCells_;

        scalar fluidVolume_;

        //- Mesh sizes the geometry was computed for (see upToDate())
        label nCellsGeom_;
        label nFacesGeom_;
        label nPointsGeom_;

        mutable autoPtr<OFstream> forcesFilePtr_;


    // Moving bodies only

        //- True if any body moves
        bool moving_;

        //- Periodic lengths of the domain (0: not periodic)
        vector periodicLengths_;

        //- Union level set at the cell centres (pressure reference choice)
        scalarField cellLevelSet_;

        //- Body velocity at the cell centres
        autoPtr<volVectorField> UbPtr_;

        //- Mass source per unit volume, S/V = (U . Sw + S_rot)/V  [1/s]
        autoPtr<volScalarField> massSourcePtr_;

        //- Face fractions of the previous geometry (for ddtCorr)
        autoPtr<surfaceScalarField> thetaOldPtr_;

        //- Diagnostics of the last update: smallest old fluid fraction and
        //  number of cells swept through by more than they hold
        scalar alphaOldMin_;
        label nSwept_;


    // Wall law of the turbulent no-slip walls with a wall function (the
    // analytic law of the wall of the SA model), demand-driven

        //- Cell sizes in wall units, h+ = u_tau h/nu, between which the
        //  weight of the wall law rises from 0 (wall resolved by the mesh)
        //  to 1 (wall modelled)
        scalar hPlusResolved_;
        scalar hPlusModelled_;

        //- Distance from the wall of the velocity sample of the wall law,
        //  in wall-normal extents h of the wall cell (default 1.5; 0: the
        //  wall cell itself)
        scalar wallSampleDistance_;

        //- Wall-normal extent h of every wall cell (0 elsewhere)
        mutable scalarField wallCellSize_;

        //- Cell whose value samples the wall law of every wall cell (-1:
        //  the wall cell itself), and the distance of its fluid centroid
        //  from the wall
        mutable labelList sampleCell_;
        mutable scalarField sampleDist_;

        //- Whether the above belong to the current geometry
        mutable bool wallSamplesDone_;


    // Cell linking (merging) of small cut cells, mergeAlpha > 0 only

        //- Cells with 0 < alpha < mergeAlpha_ are linked to a fluid face
        //  neighbour (0: off)
        scalar mergeAlpha_;

        //- Link conductance relative to the larger diagonal of the pair
        scalar mergePenalty_;

        //- Face linking a cell to its neighbour on the way to its master
        //  (-1: not linked)
        labelList linkFace_;

        //- Master of every linked cell (end of the chain), -1 otherwise
        labelList linkMaster_;

        //- Distance to the wall of the centroid of the merged fluid volume
        //  (master and slaves), on every cell of a group (0 elsewhere)
        scalarField linkWallDist_;


    // Private member functions

        //- Minimum over the level sets of all level-set bodies
        void unionLevelSet(const pointField& pts, scalarField& phi) const;

        //- Fluid fraction and centroid of the wet part of every face from
        //  the level set at the mesh points and the crossing point on every
        //  cut edge
        void faceFractions
        (
            const scalarField& pointPhi,
            const pointField& edgeCross,
            scalarField& theta,
            pointField& wetCentre
        ) const;

        //- Fluid fractions from the level-set bodies
        void calcFromLevelSet
        (
            scalarField& alpha,
            scalarField& theta,
            pointField& wetCentre,
            labelList& cellBody,
            vectorField& centroidOffset,
            scalarField& pointPhi,
            pointField& edgeCross
        ) const;

        //- Fluid fractions from an alpha-field body
        void calcFromAlphaField
        (
            const alphaFieldBody& body,
            scalarField& alpha,
            scalarField& theta,
            pointField& wetCentre,
            scalarField& pointPhi,
            pointField& edgeCross
        ) const;

        //- Fluid fraction of one cell by sub-sampling the level set, and the
        //  centroid of its fluid part
        scalar sampleCellFraction(const label celli, point& centroid) const;

        //- Centroid vectors and cut delta coefficients of every face, and
        //  the offsets of the face-value correction of the gradient
        void calcCentroidGeometry();

        //- Wall pieces of every wall cell (wallMoment_, wallPieceSum_)
        void calcWallPieces();

        //- Masters and link faces of the small cut cells (mergeAlpha > 0)
        void calcLinks();

        //- Wall-normal extent and wall-law sample cell of every wall cell
        //  (on first use after a geometry update)
        void calcWallSamples() const;

        //- The open internal face from a cut cell to the fluid neighbour
        //  with more fluid, preferring the direction away from the wall and
        //  large wet areas; -1 if there is none
        label outwardFace(const label celli) const;

        //- u_tau of a wall cell from the SA law of the wall with the
        //  tangential slip (along the cell's wall) of the sample cell at
        //  the distance y: the cell itself at dWall, or sampleCell_ at
        //  sampleDist_; 0 if the cell has no wall function
        scalar wallFunctionUTau
        (
            const volVectorField& U,
            const volScalarField& nu,
            const label celli,
            const label samplei,
            const scalar y
        ) const;

        //- Everything that follows from alpha and theta
        void calcDerived();

        //- The whole geometry: fractions of all bodies, slivers, sync,
        //  derived quantities. Called at construction, at every update
        //  of a moving body and after a change of the mesh topology.
        void calcGeometry();

        //- Mass source, swept old volume and body velocity (moving only)
        void calcMovingTerms();

        //- Distance from x0 along n to the boundary of the bounding box of
        //  cell celli (solved directions only)
        scalar cellExitDistance
        (
            const label celli,
            const point& x0,
            const vector& n
        ) const;

        //- Flux of the rotational body velocity through the wall of every
        //  cell (chord and closed faces), from its open wet faces
        void rotationalWallFlux(scalarField& Srot) const;

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
            const volScalarField& wallCoeff() const { return wallCoeff_; }
            const volScalarField& noSlipCoeff() const { return noSlipCoeff_; }
            const volScalarField& blankCoeff() const { return blankCoeff_; }
            const labelList& cellBody() const { return cellBody_; }
            const bitSet& solidCells() const { return solidCells_; }
            const bitSet& cutCells() const { return cutCells_; }
            const bitSet& wallCells() const { return wallCells_; }
            scalar alphaMin() const { return alphaMin_; }

            //- Face fractions of the geometric cut (all mesh faces), before
            //  sliver absorption and coupled-face synchronisation
            const scalarField& thetaGeo() const { return thetaGeo_; }

            //- Centroids of the geometric wet faces (all mesh faces; the
            //  face centre where the face is wholly wet or dry)
            const pointField& wetCentreGeo() const { return wetCentreGeo_; }

            //- Total fluid volume, sum(alpha V), global
            scalar fluidVolume() const { return fluidVolume_; }

            //- Signed distance to the union of the bodies at the given
            //  points (positive in the fluid), from the level sets, which
            //  are distance functions for every level-set body type. Fatal
            //  error if a body has no level set (alpha-field body).
            void wallDistance(const pointField& pts, scalarField& d) const;

            //- Fluid centroids of all cells (cell centre in full and solid
            //  cells), whether or not the second-order treatment is on
            tmp<pointField> fluidCentroids() const;


        // Second-order wall treatment

            bool secondOrder() const { return secondOrder_; }

            const volVectorField& centroidOffset() const { return centroidOffset_; }
            const surfaceVectorField& wetFaceOffset() const { return wetFaceOffset_; }
            const surfaceVectorField& dCentroid() const { return dCentroid_; }
            const surfaceScalarField& cutDeltaCoeffs() const { return cutDeltaCoeffs_; }
            const surfaceVectorField& faceValueOffset() const { return faceValueOffset_; }

            //- sum_k (x_k - x_c) S_k over the wall pieces of every cell
            const tensorField& wallMoment() const { return wallMoment_; }

            //- sum_k S_k over the wall pieces (equals Sw by closure)
            const vectorField& wallPieceSum() const { return wallPieceSum_; }

            //- Least-squares gradient about the fluid centroids, over the
            //  open faces to fluid neighbours (exact for a linear field)
            tmp<volVectorField> lsqGrad(const volScalarField& f) const;

            //- Explicit remainder of the centroid-consistent face gradient:
            //  n.grad f = (f_N - f_P)/(n.d) + [g_f.n - (g_f.d)/(n.d)],
            //  the bracket, with g_f the interpolated least-squares gradient
            //- Gradient averaged over the fluid part of every cell, with a
            //  Dirichlet value on the immersed wall:
            //  (sum_f theta_f vf_f Sf + vf_w Sw)/(alpha V), zero in the solid.
            //  For quantities of the turbulence models (velocity gradient
            //  with the wall velocity, nuTilda with zero), where fvc::grad
            //  would treat the immersed wall as zero-gradient.
            template<class Type>
            tmp
            <
                GeometricField
                <
                    typename outerProduct<vector, Type>::type,
                    fvPatchField,
                    volMesh
                >
            >
            fluidGrad
            (
                const GeometricField<Type, fvPatchField, volMesh>& vf,
                const Field<Type>& wallValue
            ) const;

            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
            snGradCorrection
            (
                const GeometricField<Type, fvPatchField, volMesh>& vf
            ) const;


        // Moving bodies

            //- True if any body moves
            bool moving() const { return moving_; }

            //- Move the bodies to the current time and rebuild the
            //  geometry and the moving-wall terms. No-op for static bodies.
            void update();

            //- Whether the geometry belongs to the current mesh: false
            //  between a change of the mesh topology (adaptive refinement)
            //  and updateMesh()
            bool upToDate() const;

            //- Rebuild the geometry after a change of the mesh topology,
            //  with the bodies where they are (level-set bodies only)
            void updateMesh();

            //- Body velocity at the cell centres
            const volVectorField& Ub() const;

            //- Mass source per unit volume S/V, to be added to the
            //  right-hand side of the pressure equation
            const volScalarField& massSource() const;

            //- Momentum source of the moving no-slip wall and of the
            //  blanking, nu (noSlipCoeff + blankCoeff) u_b
            tmp<volVectorField> movingWallSource
            (
                const volScalarField& nu
            ) const;

            //- Same, with the wall viscosity of a turbulent wall treatment
            //  for the traction and nu for the blanking
            tmp<volVectorField> movingWallSource
            (
                const volScalarField& nuWall,
                const volScalarField& nu
            ) const;

            //- Pressure reference cell: the user's cell if it is fluid,
            //  otherwise the fluid cell farthest from the bodies
            label pressureRefCell(const label userRefCell) const;

            //- Smallest old fluid fraction of the last update
            scalar alphaOldMin() const { return alphaOldMin_; }

            //- Cells swept through by more than they hold in the last update
            label nSwept() const { return nSwept_; }


        // Operators

            //- Cut-cell Gauss gradient, sum_f theta_f (p_f - p_P) Sf / V.
            //  Exactly zero for a uniform field; includes the wall force.
            //  With secondOrder: face values moved from the centroid
            //  midpoint to the wet-face centroid and the wall pressure
            //  integrated over the wall pieces, exact for a linear field
            //  (the result is then V_geo grad p / V).
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

            //- Viscosity of the no-slip wall traction nu_w Awall/dWall (u_P - u_b)
            //  for turbulent flow, the immersed counterpart of a nut wall
            //  function: nu everywhere, and on the wall cells of no-slip bodies
            //  with a wall function (wallFunction spalartAllmaras, default)
            //  nu_w = u_tau^2 dWall/|u_t|, |u_t| the cell's tangential slip and
            //  u_tau that of the wall law (wallLaw), capped at
            //  nu + kappa u_tau max(h, dWall)
            tmp<volScalarField> wallViscosity
            (
                const volVectorField& U,
                const volScalarField& nu
            ) const;

            //- Wall law of the wall cells of no-slip bodies with a wall
            //  function: the weight w of the wall law against the solution,
            //  0 for h+ <= hPlusResolved, 1 for h+ >= hPlusModelled and
            //  smooth in log h+ in between (h the extent of the cell along
            //  the wall normal), and u_tau^2 = (1 - w) u_tau,cell^2
            //  + w u_tau,sample^2 from the SA law of the wall with the wall
            //  cell's own value at dWall and with the sample at about
            //  wallSampleDistance h from the wall (the wall cell itself if
            //  wallSampleDistance h <= dWall, or the sample is in the cell's
            //  linked group, in the solid or on another processor). Both 0
            //  in all other cells
            void wallLaw
            (
                const volVectorField& U,
                const volScalarField& nu,
                scalarField& uTau,
                scalarField& weight
            ) const;

            //- Link the small cut cells to their masters in an assembled
            //  equation: a conductance c on the link face, so that the sum of
            //  the two rows is the equation of the merged cell (c cancels)
            //  and the two values coincide as c grows (cell merging without
            //  changing the matrix sparsity). Apply after relax(), before
            //  solve(). No-op unless mergeAlpha > 0.
            template<class Type>
            void linkCells(fvMatrix<Type>& eqn) const;

            //- Replace the values of every merged group (master and slaves)
            //  by their fluid-volume-weighted average: for cell averages of
            //  the merged control volume, e.g. a gradient over the fluid,
            //  whose contributions of the faces inside the group cancel.
            //  No-op unless mergeAlpha > 0.
            template<class Type>
            void groupAverage(Field<Type>& f) const;

            //- Replace the values of every merged group by their
            //  cell-volume-weighted average: for per-cell-volume quantities
            //  such as ibm.grad(p) (a force over V).
            //  No-op unless mergeAlpha > 0.
            template<class Type>
            void groupMean(GeometricField<Type, fvPatchField, volMesh>& f)
            const;

            //- Diagonal coefficient A (or H1) of the merged control volumes
            //  from the unlinked matrix: (sum_i V_i A_i + sum_l (upper_l
            //  + lower_l))/sum_i V_i over the cells i and the internal link
            //  faces l of the group, on every cell of the group. With mergeH,
            //  the pressure-velocity coupling of the merged cell; rAU from
            //  the linked matrix would hold the penalty and decouple the
            //  pressure along a strip of linked cells.
            //  No-op unless mergeAlpha > 0.
            template<class Type>
            void mergeA(const fvMatrix<Type>& eqn, volScalarField& A) const;

            //- H of the merged control volumes from the unlinked matrix,
            //  (sum_i V_i H_i + sum_l (upper_l psi_nei + lower_l psi_own))
            //  /sum_i V_i: the coupling through the internal link faces
            //  removed. No-op unless mergeAlpha > 0.
            template<class Type>
            void mergeH
            (
                const fvMatrix<Type>& eqn,
                GeometricField<Type, fvPatchField, volMesh>& H
            ) const;

            //- Master of every linked cell, -1 otherwise
            const labelList& linkMaster() const { return linkMaster_; }

            //- True if small cut cells are linked (mergeAlpha > 0)
            bool linking() const { return mergeAlpha_ > 0; }

            //- Wall distance of the merged centroid on the cells of a group
            //  (master and slaves), 0 elsewhere
            const scalarField& linkWallDist() const { return linkWallDist_; }

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
