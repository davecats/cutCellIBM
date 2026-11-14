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
    mergeAlpha    0.5;    // cells with less fluid are linked to a master
                          // (merged control volume; 0: no linking)
    alphaMin      0;      // cells with less fluid are absorbed into the
                          // body (default 0 with linking: only cells
                          // without fluid; 0.02 without)
    nSubSamples   32;     // per direction, for the cell fluid fraction
    nBisect       40;     // bisection steps locating the wall on cut edges
    writeFields   true;   // write ibmAlpha, ibmTheta, ... at start-up
    writeForces   true;   // postProcessing/ibmForces/<time>/forces.dat and
                          // moments.dat (torque about each body's centre)
    zeroShear     false;  // free-slip walls with zero tangential stress
                          // (curvature friction, zeroShearSource); false:
                          // du/dn = 0, not stress-free on a curved wall

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

    Invariant bodies (motion with "invariant yes", see bodyMotion): a body
    whose motion leaves it unchanged as a set (a cylinder turning about or
    sliding along its axis, a sphere turning about its centre) is not
    moving. Its geometry is built once, there is no wall flux, and only
    ibmUb exists: update() moves the bodies to the current time and
    re-evaluates u_b, which enters the no-slip traction, the blanking, the
    forces and the wall functions exactly as for a moving body. The
    invariance is checked at construction, by moving the body over a
    virtual time in which its wall points travel one cell size, and at
    every update for the actual motion state: the level set at the edge
    crossings of the cut must not change by more than
    invarianceTolerance (default 1e-6) times the edge length. Alpha-field
    bodies cannot be invariant (no level set).

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
#include "labelPair.H"
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

        //- Convective face values of cut and linked cells reconstructed from
        //  the value point of the upwind cell (ibmConvection, default secondOrder)
        bool ibmConvection_;

        //- With ibmConvection, the tangential velocity on the faces of wall
        //  cells with a wall function from the SA law of the wall instead of
        //  the least-squares gradient, with the weight of the wall law
        //  (ibmConvectionWallLaw, default false: no clear gain, see
        //  doc/secondOrderStudy.md)
        bool ibmConvectionWallLaw_;

        //- Zero tangential stress on the free-slip walls (zeroShear, default
        //  false; doc/freeSlipPlan.md). Without it a free-slip wall has no
        //  wall flux in the Laplacian form, i.e. du/dn = 0, which on a curved
        //  wall is not stress-free: P du/dn = H.(u - u_b)_t + omega x n is
        //  (H the curvature tensor of the wall, positive on a convex body,
        //  n into the fluid, omega the angular velocity of the body)
        bool zeroShear_;

        //- With zeroShear and secondOrder: the tangential velocity at the
        //  wall-piece centroid, u_P + delta.grad u (delta the tangential
        //  offset of the wall centroid from the value point, lsqGrad),
        //  instead of u_P (zeroShearTangential, default true)
        bool zeroShearTangential_;

        //- Fluid centroid minus cell centre
        volVectorField centroidOffset_;

        //- Point the cell value belongs to minus cell centre: the fluid
        //  centroid, and on every cell of a linked group (master and
        //  slaves, which share one value) the centroid of the group
        volVectorField valueOffset_;

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

        //- sum_k (x_k - x_c) ((x_k - x_c) ^ S_k): moment of the linear part
        //  of the wall pressure about the value point (torque)
        tensorField wallTorqueMoment_;

        //- Zero-shear free slip (zeroShear_): the free-slip wall cells of
        //  bodies with a smooth level set, their index in zsCells_ (-1
        //  elsewhere), and per such cell the friction tensor
        //  (Awall/V) H_w (I + d H_w)^-1 [1/m^2], its implicit part (largest
        //  |eigenvalue|) and the map (I + d H_w)^-1 P from the relative
        //  velocity at the value point to its tangential wall value
        //  (d = dWall, H_w the curvature tensor at the wall)
        labelList zsCells_;
        labelList zsIndex_;
        symmTensorField zsFriction_;
        scalarField zsImplicit_;
        symmTensorField zsWallMap_;
        vectorField zsOffset_;
        bool zsReported_ = false;

        //- Work arrays filled by faceFractions/calcFromLevelSet
        mutable vectorField wetOffsetAll_;

        labelList cellBody_;

        //- Nearest level-set body of every cell (-1 if none)
        mutable labelList nearestBody_;

        bitSet solidCells_;
        bitSet cutCells_;
        bitSet wallCells_;

        scalar fluidVolume_;

        //- False from a change of the mesh topology until updateMesh()
        bool upToDate_;

        mutable autoPtr<OFstream> forcesFilePtr_;
        mutable autoPtr<OFstream> momentsFilePtr_;


    // Moving bodies only

        //- True if any body moves (the geometry changes)
        bool moving_;

        //- True if any wall has a velocity (moving or invariant bodies):
        //  ibmUb exists
        bool wallVelocity_;

        //- Largest change of the level set of an invariant body at the
        //  edge crossings, relative to the edge length
        scalar invarianceTolerance_;

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

        //- Solid cells of the previous geometry (cells just uncovered)
        bitSet solidOld_;

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
        //  the wall cell itself or a sample on another processor), and the
        //  distance of its fluid centroid from the wall
        mutable labelList sampleCell_;
        mutable scalarField sampleDist_;

        //- Samples on other processors: the slot of every wall cell (-1:
        //  none), and per processor the local cells whose values go there
        //  and the slots that receive values from there
        mutable labelList sampleSlot_;
        mutable label nSampleSlots_;
        mutable labelListList sampleSendCells_;
        mutable labelListList sampleRecvSlots_;

        //- Whether the above belong to the current geometry
        mutable bool wallSamplesDone_;


    // Least-squares gradient (lsqGrad): its geometric part, demand-driven,
    // rebuilt on the first call after the geometry or the mesh changes

        //- Whether the data below belong to the current geometry
        mutable bool lsqDone_;

        //- Open internal faces of the fit, in face order, and their
        //  weighted centroid vectors w d, w = 1/|d|^2
        mutable labelList lsqFaces_;
        mutable vectorField lsqWd_;

        //- The same on the coupled patches (patch-local face indices)
        mutable labelListList lsqPatchFaces_;
        mutable List<vectorField> lsqPatchWd_;

        //- Inverse of the regularised normal matrix of every cell with a fit
        mutable symmTensorField lsqInvA_;
        mutable bitSet lsqOk_;

        //- Fluid cells without a fit (slots), which take the mean gradient
        //  of their fluid face neighbours with one: the cell and the number
        //  of those neighbours of every slot, (slot, neighbour cell) over
        //  the internal faces in face order, and per coupled patch
        //  (slot, patch face)
        mutable labelList lsqFallbackCells_;
        mutable labelList lsqFallbackCount_;
        mutable labelPairList lsqFallback_;
        mutable List<labelPairList> lsqPatchFallback_;

        //- Cells whose gradient the second-order corrections read (cells
        //  with a value offset, wall cells, linked groups, the cells of the
        //  faces with a correction and the fallback cells with their
        //  neighbours), and the entries of lsqFaces_/lsqPatchFaces_ that
        //  touch them
        mutable labelList lsqNearCells_;
        mutable labelList lsqNearFaces_;
        mutable labelListList lsqNearPatchFaces_;

        //- Faces on which the centroid-consistent face-gradient correction
        //  (snGradCorrection) does not vanish: |n - d/(n.d)| > snGradTol;
        //  coupled faces if either side needs it
        mutable labelList snGradFaces_;
        mutable labelListList snGradPatchFaces_;


    // Cell linking (merging) of small cut cells, mergeAlpha > 0 only

        //- Cells with 0 < alpha < mergeAlpha_ are linked to a fluid face
        //  neighbour (0: off; default 0.5)
        scalar mergeAlpha_;

        //- Link conductance relative to the larger diagonal of the pair
        scalar mergePenalty_;

        //- Face linking a cell to its neighbour on the way to its master
        //  (-1: not linked): an internal face or a coupled boundary face
        //  (processor, cyclic), the master may be on another processor
        labelList linkFace_;

        //- Group of every cell of a merged group (master and slaves): the
        //  global index of the master cell (globalIndex of the cells), the
        //  same in every decomposition up to the numbering; -1 outside
        labelList linkGroup_;

        //- Group of the cell across every coupled boundary face (-1 for
        //  cells outside a group and for the other boundary faces)
        labelList nbrLinkGroup_;

        //- Coupled boundary faces that carry a link (from either side)
        boolList boundaryLink_;

        //- Fluid centroid of every group cell relative to that of its master
        //  (in the frame of the cell: cyclics are translations)
        vectorField linkRel_;

        //- Distance to the wall of the centroid of the merged fluid volume
        //  (master and slaves), on every cell of a group (0 elsewhere)
        scalarField linkWallDist_;

        //- Local slot of the group of every cell (-1 outside), the group
        //  of every slot, and the exchange of the group sums: per
        //  processor, the slots whose partial sums go to the processor of
        //  the master (sendSlots_) and the slots of local masters that
        //  receive partial sums from it (recvSlots_)
        labelList groupSlot_;
        labelList slotGroup_;
        labelListList sendSlots_;
        labelListList recvSlots_;

        //- Cells of the split groups, the groups with a link over a
        //  coupled face (processor, cyclic), and the number of these groups
        //  on all processors
        labelList splitGroupCells_;
        label nSplitGroups_;


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

        //- Fluid centroid of one cell of an alpha-field body: sub-samples of
        //  the trilinear interpolant of the point level set (its edge
        //  crossings are those of the face polygons). Axis-aligned hexahedra
        //  only; false (centroid unchanged) for other cells.
        bool sampleInterpolatedCentroid
        (
            const label celli,
            const scalarField& pointPhi,
            point& centroid
        ) const;

        //- Warn if the laplacian/snGrad schemes do not match secondOrder
        //  (cutCellCorrected with the switch on, another scheme without)
        void checkSchemes() const;

        //- Value points (valueOffset_): fluid centroids, group centroids on
        //  the cells of linked groups
        void calcValueOffset();

        //- True for a face between two cells of the same linked group
        bool insideGroup(const label facei) const;

        //- Centroid vectors and cut delta coefficients of every face, and
        //  the offsets of the face-value correction of the gradient, from
        //  the value points (from the fluid centroids inside a group)
        void calcCentroidGeometry();

        //- Wall pieces of every wall cell (wallMoment_, wallPieceSum_)
        void calcWallPieces();

        //- Curvature friction of the zero-shear free-slip walls (zsCells_,
        //  zsFriction_, zsImplicit_, zsWallMap_, zsOffset_)
        void calcZeroShear();

        //- Relative velocity u - u_b and body velocity u_b of zsCells_ at
        //  the value point, moved by zsOffset_ along the wall (second order
        //  with zeroShearTangential: u_P + delta.grad u); u_b zero for
        //  static bodies
        void zeroShearRelative
        (
            const volVectorField& U,
            vectorField& wP,
            vectorField& ubP
        ) const;

        //- Masters and link faces of the small cut cells (mergeAlpha > 0)
        void calcLinks();

        //- Wall-normal extent and wall-law sample cell of every wall cell
        //  (on first use after a geometry update)
        void calcWallSamples() const;

        //- The geometric part of lsqGrad, the cells and faces the
        //  second-order corrections read (lsqDone_)
        void calcLsq() const;

        //- lsqGrad on every cell (near false), or on lsqNearCells_ only,
        //  zero elsewhere (near true; the same values there)
        tmp<volVectorField> lsqGrad
        (
            const volScalarField& f,
            const bool near
        ) const;

        //- The open internal or coupled face from a cut cell to the fluid
        //  neighbour with more fluid, preferring the direction away from
        //  the wall and large wet areas, ties decided by the face centre;
        //  -1 if there is none. nbrAlpha: alpha across the boundary faces.
        label outwardFace
        (
            const label celli,
            const scalarField& nbrAlpha
        ) const;

        //- Sum over every group: on entry the contribution of every cell, on
        //  exit the sum over its group on every cell of a group (also
        //  across processors); the other cells keep their values
        template<class Type>
        void groupSum(UList<Type>& f) const;

        //- u_tau of a wall cell from the SA law of the wall with the
        //  tangential slip (along the cell's wall) of the sample velocity
        //  Us at the distance y: the cell itself at dWall, or the sample
        //  at sampleDist_; 0 if the cell has no wall function
        scalar wallFunctionUTau
        (
            const vector& Us,
            const volScalarField& nu,
            const label celli,
            const scalar y
        ) const;

        //- Wall-law profile of the velocity in the wall cells, for the
        //  convective face values of U: u_tau and weight of the wall law
        //  (wallLaw), the wall normal (out of the fluid) and the body
        //  velocity, with their values across the coupled faces; unset
        //  without a wall-function cell
        struct wallProfileData
        {
            autoPtr<volScalarField> uTau, w;
            autoPtr<volVectorField> n, Ub;
            bool valid() const { return uTau.valid(); }
        };

        //- wallProfileData of a velocity field; nothing for other fields
        template<class Type>
        void wallProfile
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const volScalarField&,
            wallProfileData&
        ) const
        {}

        void wallProfile
        (
            const volVectorField& U,
            const volScalarField& nu,
            wallProfileData& d
        ) const;

        //- Face value e (extrapolated with the least-squares gradient) of a
        //  wall cell with value v at the distance y from the wall, r from
        //  the value point to the wet-face centroid: (1 - w) e + w e_law,
        //  e_law with the tangential slip v_t u+(y_f+)/u+(y+) of the SA law,
        //  y_f = y - r.n, and the normal component of e. Velocities only;
        //  other fields keep e
        template<class Type>
        Type wallProfileValue
        (
            const Type& e,
            const Type&,
            const vector&,
            const vector&,
            const vector&,
            const scalar,
            const scalar,
            const scalar,
            const scalar
        ) const
        {
            return e;
        }

        vector wallProfileValue
        (
            const vector& e,
            const vector& v,
            const vector& r,
            const vector& n,
            const vector& Ub,
            const scalar y,
            const scalar uTau,
            const scalar w,
            const scalar nu
        ) const;

        //- Everything that follows from alpha and theta
        void calcDerived();

        //- The whole geometry: fractions of all bodies, slivers, sync,
        //  derived quantities. Called at construction, at every update
        //  of a moving body and after a change of the mesh topology.
        void calcGeometry();

        //- Mass source, swept old volume and body velocity (moving only)
        void calcMovingTerms();

        //- Body velocity u_b at the wall points of the cells a body with a
        //  wall velocity owns, zero elsewhere (ibmUb)
        void calcWallVelocity();

        //- Check that the invariant bodies keep their shape under their
        //  motion state advanced over a time tau (fatal otherwise)
        void checkInvariance(const scalar tau) const;

        //- Virtual time over which the fastest wall point of an invariant
        //  body travels one edge length of the cut (start-up check)
        scalar invarianceTime() const;

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
            const volVectorField& valueOffset() const { return valueOffset_; }
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

            bool ibmConvection() const { return ibmConvection_; }

            //- Zero tangential stress on the free-slip walls (zeroShear)
            bool zeroShear() const { return zeroShear_; }

            //- Wall flux of the zero-shear free-slip walls in the Laplacian
            //  form, nu (Awall/V) (H_w.(u_w - u_b,w)_t + omega x n), with the
            //  tangential wall value from the value point by the Robin factor
            //  (I + d H_w)^-1 (exact for a rigidly rotating flow): implicit
            //  with the largest |eigenvalue| of the friction tensor, the rest
            //  explicit. Add to the momentum matrix when zeroShear() is true.
            //  Bodies without a smooth level set (triSurfaceMesh, alpha
            //  field) keep du/dn = 0.
            tmp<fvVectorMatrix> zeroShearSource
            (
                const volScalarField& nu,
                const volVectorField& U
            ) const;

            //- Explicit correction of a convection term on the faces of cut
            //  and linked cells: phi (v_f - vSchemeFace) with v_f the value
            //  of the upwind cell at its value point extrapolated to the
            //  wet-face centroid with the least-squares gradient,
            //  v_U + g_U.(x'_f - x_v,U); zero elsewhere and without
            //  ibmConvection. Add fvc::surfaceIntegrate of it to the equation.
            //  With the viscosity nu (velocity of a turbulent flow) and
            //  ibmConvectionWallLaw, the tangential velocity on the faces of
            //  wall cells with a wall function follows the SA law of the
            //  wall with the weight of the wall law (wallProfileValue).
            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
            convectionCorrection
            (
                const surfaceScalarField& phi,
                const GeometricField<Type, fvPatchField, volMesh>& vf,
                const GeometricField<Type, fvsPatchField, surfaceMesh>& vSchemeFace,
                const volScalarField* nu = nullptr
            ) const;

            //- convectionCorrection against the case's scheme div(phi,vf),
            //  integrated over the cells: add it to the matrix of
            //  fvm::div(phi, vf) when ibmConvection() is set
            template<class Type>
            tmp<GeometricField<Type, fvPatchField, volMesh>>
            convectionSource
            (
                const surfaceScalarField& phi,
                const GeometricField<Type, fvPatchField, volMesh>& vf,
                const volScalarField* nu = nullptr
            ) const;

            //- convectionSource of the velocity of a turbulent flow, with
            //  the wall law on the faces of wall cells with a wall function
            template<class Type>
            tmp<GeometricField<Type, fvPatchField, volMesh>>
            convectionSource
            (
                const surfaceScalarField& phi,
                const GeometricField<Type, fvPatchField, volMesh>& vf,
                const volScalarField& nu
            ) const
            {
                return convectionSource(phi, vf, &nu);
            }

            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
            snGradCorrection
            (
                const GeometricField<Type, fvPatchField, volMesh>& vf
            ) const;


        // Moving bodies

            //- True if any body moves (the geometry changes)
            bool moving() const { return moving_; }

            //- True if any wall has a velocity: moving or invariant bodies.
            //  Ub() and movingWallSource() exist exactly then.
            bool wallVelocity() const { return wallVelocity_; }

            //- Move the bodies to the current time and rebuild the
            //  geometry and the moving-wall terms; for invariant bodies
            //  only the wall velocity. No-op for static bodies.
            void update();

            //- Whether the geometry belongs to the current mesh: false
            //  between a change of the mesh topology (adaptive refinement)
            //  and updateMesh()
            bool upToDate() const { return upToDate_; }

            //- Mark the geometry as belonging to a previous mesh (called by
            //  the adaptive mesh before every change of the topology)
            void outOfDate() { upToDate_ = false; }

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
            //  replacement of interpolate(rAU)*fvc::ddtCorr(U, phi): the
            //  coefficient of U.oldTime() in ddt(alpha, U) is the old fluid
            //  fraction (the group mean in a merged control volume), as rho
            //  in the compressible ddtCorr(rho, U, phi)
            tmp<surfaceScalarField> ddtCorr
            (
                const volScalarField& rAU,
                const volVectorField& U,
                const surfaceScalarField& phi
            ) const;

            //- One value of the flux on each coupled face: the two copies
            //  (processor, cyclic) set to their antisymmetric mean. Nothing in
            //  the pressure correction removes a difference between them, and
            //  the ddt correction carries it to the next time step.
            void syncCoupledFlux(surfaceScalarField& phi) const;

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

            //- Diagonal of a matrix with the coefficients of the coupled
            //  faces (processor, cyclic), which the matrix holds in
            //  internalCoeffs and an internal face in the diagonal: the same
            //  in every decomposition
            template<class Type>
            tmp<scalarField> coupledDiag(const fvMatrix<Type>& eqn) const;

            //- Values of the cells the moving bodies have just uncovered
            //  (solid before the last update, fluid now): the mean over
            //  their open faces of the fluid neighbours (also across coupled
            //  faces), layer by layer. For the pressure, which is defined
            //  up to a constant: the solid value would make a jump of the
            //  size of the pressure level (set by the reference cell) seen
            //  by grad(p) before the first pressure solve of the step.
            template<class Type>
            void fillUncovered
            (
                GeometricField<Type, fvPatchField, volMesh>& f
            ) const;

            //- Link the small cut cells to their masters in an assembled
            //  equation: a conductance c on the link face, so that the sum of
            //  the two rows is the equation of the merged cell (c cancels)
            //  and the two values coincide as c grows (cell merging without
            //  changing the matrix sparsity). Apply after relax(), before
            //  solve(). groupUniform (the momentum predictor, whose solution
            //  the pressure corrector replaces by values uniform in the
            //  groups): the residual of every group cell at the current psi,
            //  less its volume share of the group residual, goes to the
            //  source, so that a converged solution is uniform in the
            //  groups and depends neither on the penalty nor on the
            //  convergence of the linear solver (steady SIMPLE: the
            //  decomposition). No-op unless mergeAlpha > 0.
            template<class Type>
            void linkCells
            (
                fvMatrix<Type>& eqn,
                const bool groupUniform = false
            ) const;

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
            //  + lower_l))/sum_i V_i over the cells i and the link faces l
            //  of the group (on a coupled link face each side adds the
            //  coupling of its own row, -boundaryCoeffs), on every cell of
            //  the group. With mergeH,
            //  the pressure-velocity coupling of the merged cell; rAU from
            //  the linked matrix would hold the penalty and decouple the
            //  pressure along a strip of linked cells.
            //  No-op unless mergeAlpha > 0.
            template<class Type>
            void mergeA(const fvMatrix<Type>& eqn, volScalarField& A) const;

            //- H of the merged control volumes from the unlinked matrix,
            //  (sum_i V_i H_i + sum_l (upper_l psi_nei + lower_l psi_own))
            //  /sum_i V_i: the coupling through the link faces removed. No-op unless mergeAlpha > 0.
            template<class Type>
            void mergeH
            (
                const fvMatrix<Type>& eqn,
                GeometricField<Type, fvPatchField, volMesh>& H
            ) const;

            //- Group (global index of the master) of every cell of a merged
            //  group, -1 otherwise
            const labelList& linkGroup() const { return linkGroup_; }

            //- Face linking every cell to its neighbour on the way to its
            //  master (an internal or a coupled boundary face), -1 if the
            //  cell is not linked
            const labelList& linkFace() const { return linkFace_; }

            //- Group of the cell across every coupled boundary face, -1
            //  for cells outside a group and on the other boundary faces
            const labelList& nbrLinkGroup() const { return nbrLinkGroup_; }

            //- Coupled boundary faces that carry a link (from either side)
            const boolList& boundaryLink() const { return boundaryLink_; }

            //- Local slot of the group of every cell, -1 outside a group
            const labelList& groupSlot() const { return groupSlot_; }

            //- Number of local slots (groups with a cell on this processor)
            label nGroupSlots() const { return slotGroup_.size(); }

            //- Cells of the split groups: the groups with a link over a
            //  coupled face (processor, cyclic), which the smoothers and
            //  preconditioners treat explicitly
            const labelList& splitGroupCells() const
            {
                return splitGroupCells_;
            }

            //- Number of split groups on all processors
            label nSplitGroups() const { return nSplitGroups_; }

            //- Complete the sums over the groups held in the local slots:
            //  on entry the sum over the local cells of every group, on
            //  exit the sum over all its cells on every processor
            template<class Type>
            void groupSlotSum(List<Type>& s) const;

            //- Link conductance relative to the larger diagonal of the pair
            scalar mergePenalty() const { return mergePenalty_; }

            //- True if small cut cells are linked (mergeAlpha > 0)
            bool linking() const { return mergeAlpha_ > 0; }

            //- Wall distance of the merged centroid on the cells of a group
            //  (master and slaves), 0 elsewhere
            const scalarField& linkWallDist() const { return linkWallDist_; }

            //- Force on every body (per unit density): pressure and viscous.
            //  nu: viscosity of the wall traction (with a wall law the wall
            //  viscosity); nuMolecularPtr: molecular viscosity of the
            //  transpose part of the stress on rotating bodies (default nu)
            void forces
            (
                const volScalarField& p,
                const volVectorField& U,
                const volScalarField& nu,
                List<vector>& Fp,
                List<vector>& Fv,
                const volScalarField* nuMolecularPtr = nullptr
            ) const;

            //- Force and moment on every body (per unit density), pressure
            //  and viscous; moments about the current centre of each body
            //  (ibmBody::currentCentre: the centre of rotation if the motion
            //  gives one, otherwise the body centre), minimum image in the
            //  periodic directions
            void forces
            (
                const volScalarField& p,
                const volVectorField& U,
                const volScalarField& nu,
                List<vector>& Fp,
                List<vector>& Fv,
                List<vector>& Mp,
                List<vector>& Mv,
                const volScalarField* nuMolecularPtr = nullptr
            ) const;

            //- Print the forces and moments and append them to the forces
            //  and moments files
            void reportForces
            (
                const volScalarField& p,
                const volVectorField& U,
                const volScalarField& nu,
                const volScalarField* nuMolecularPtr = nullptr
            ) const;

            //- Append one line of per-body pressure and viscous vectors to a
            //  file in postProcessing/ibmForces/<start time> (opened on the
            //  first call, with a header)
            void writeForceFile
            (
                autoPtr<OFstream>& filePtr,
                const word& fileBaseName,
                const string& title,
                const word& symbol,
                const List<vector>& Xp,
                const List<vector>& Xv
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
