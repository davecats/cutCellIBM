/*---------------------------------------------------------------------------*\
    cutCellIBM: cut-cell immersed boundary method for OpenFOAM

Class
    Foam::ibmRefineFvMesh

Description
    Adaptive refinement of a hexahedral background mesh for cut-cell
    immersed boundaries, in 2-D (one cell thick, empty front and back) and
    3-D. Every adaptation brings each cell to a target level:

      - the cut cells and nCutCellLayers layers of cells around them (the
        cells closer to a body than their own half diagonal plus
        nCutCellLayers finest cell sizes) are at maxLevel, the finest level;
      - optional refinement regions raise the level elsewhere, either by
        the distance from the bodies or by the value of a solution field
        (a vector or tensor field by its magnitude); a region's level is at
        most maxLevel;
      - everything else may be coarsened back towards level 0.

    The refinement engine (hexRef3D in 3-D, hexRef4 in 2-D, with the 2:1
    balance across faces) is taken from multiDimAMR, see ../hexRef/README.

    The first update of a run (the first time step or iteration) refines
    until every cell has reached its target, which builds the mesh around
    the bodies from a uniform background mesh. Later updates, every
    refineInterval time steps (iterations of a steady solver), refine and
    coarsen to the targets of the current solution; refineInterval 0 keeps
    the mesh of the first update. After every change of the topology the
    cut-cell geometry and the wall distance are rebuilt, and the solver
    recomputes the flux (mesh.topoChanging()).

    Load balancing (parallel runs, enableBalancing true): after an
    adaptation that changed the mesh, if the number of cells per processor
    exceeds its mean by more than allowedImbalance, the mesh and all fields
    are redistributed (fvMeshDistribute). The decomposition (balanceMethod,
    default scotch, coefficients from system/decomposeParDict if needed)
    keeps every cell that descends from one cell of the background mesh on
    one processor, so that the refinement history stays complete and the
    cells can still be coarsened. The mesh is redistributed only if the
    new decomposition is balanced (imbalance at most allowedImbalance) or
    at least halves the imbalance. The redistribution comes before the
    rebuild of the cut-cell geometry, which is done once.

    constant/dynamicMeshDict (re-read at every update):
    \verbatim
    dynamicFvMesh   ibmRefineFvMesh;

    ibmRefineFvMeshCoeffs
    {
        maxLevel        4;      // finest level: the cut cells
        nCutCellLayers  2;      // finest cells on either side of the walls
        refineInterval  10;     // 0: only the first update
        maxCells        1000000;// no solution-based refinement beyond
        unrefine        true;   // coarsen where the targets have dropped
        nBufferLayers   1;      // around the solution-based regions

        enableBalancing  true;  // optional (default false)
        allowedImbalance 0.2;   // max/mean cells per processor - 1
        balanceMethod    scotch;// optional (default scotch)

        regions
        {
            near
            {
                type        distance;   // |distance to the bodies| < distance
                distance    0.1;
                level       2;
            }
            wake
            {
                type        field;      // lowerBound <= value <= upperBound
                field       vorticity;  // e.g. from a function object
                lowerBound  10;
                upperBound  1e15;
                level       3;
            }
        }
    }
    \endverbatim

    Bodies must be level-set bodies (analytic shapes or surfaces): their
    signed distance places the refinement, and alpha-field bodies cannot be
    rebuilt on a changed mesh.

SourceFiles
    ibmRefineFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef ibmRefineFvMesh_H
#define ibmRefineFvMesh_H

#include "dynamicFvMesh.H"
#include "hexRef.H"
#include "bitSet.H"

namespace Foam
{

class cutCellGeometry;
class polyTopoChange;

class ibmRefineFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Refinement engine (cell and point levels, history)
        autoPtr<hexRef> meshCutter_;

        //- Number of adaptations done
        label nAdapt_;

        //- Faces created or changed by the current adaptation
        bitSet changedFaces_;

        //- True while fvMeshDistribute changes the mesh: the refinement
        //  levels and history are distributed afterwards, as a whole
        bool distributing_;


    // Private member functions

        //- The coefficients, re-read from constant/dynamicMeshDict
        dictionary coeffs() const;

        //- The cut-cell geometry of the case
        cutCellGeometry& ibm() const;

        //- Distance from the cell centre to the farthest cell point in the
        //  solved directions (the radius of the cell)
        scalarField cellRadius() const;

        //- Field value of the region "field" per cell, false if the field
        //  does not exist (yet)
        bool regionField(const word& name, scalarField& value) const;

        //- Raise every cell to the highest level among its face neighbours,
        //  nLayers times
        void dilate(labelList& level, const label nLayers) const;

        //- Target level of every cell (refine below it) and the level it
        //  keeps (no coarsening below it, a band wider than the target
        //  against cells being coarsened and refined in turn). False if
        //  the field of a region does not exist (yet).
        bool targetLevels
        (
            const dictionary& dict,
            const bool solution,
            labelList& refineLevel,
            labelList& keepLevel
        ) const;

        //- Apply a topology change, map the fields, mark the changed faces;
        //  returns the old number of cells
        label changeTopology(polyTopoChange& meshMod);

        //- Split the cells, map the fields
        void refine(const labelList& cells);

        //- Merge the cells around the split points/edges, map the fields
        void unrefine(const labelList& splitElems);

        //- The fluxes and their velocities (correctFluxes, default phi U)
        static List<Pair<word>> fluxNames(const dictionary& dict);

        //- |phi + phi_nbr|/(|phi| + |phi_nbr|) on the boundary faces, zero
        //  on the uncoupled ones: the two sides of a coupled face carry
        //  opposite fluxes
        tmp<scalarField> coupledFluxMismatch
        (
            const surfaceScalarField& phi
        ) const;

        //- Fluxes on the changed faces from the mapped velocity
        void correctFluxes(const dictionary& dict);

        //- Wall distance of the turbulence models, after the cut-cell
        //  geometry
        void updateWallDistance();

        //- Number of cells of every processor, and the imbalance:
        //  max/mean - 1
        static scalar imbalance(const labelList& nCellsProc);
        labelList nCellsProc() const;

        //- Evaluate the coupled patches of the volume fields of one type
        //  (fvMeshDistribute leaves the new processor patches zero)
        template<class Type>
        void correctCoupledPatches();

        //- Redistribute the mesh and the fields if the imbalance exceeds
        //  allowedImbalance and a better decomposition exists. Returns
        //  true if the mesh has been redistributed.
        bool balance(const dictionary& dict);

        //- Debug: sums over the cells that a redistribution keeps
        void printFieldSums(const dictionary& dict) const;

        //- Debug: compare the processor patches with their neighbours
        void checkProcessorPatches() const;

        //- Warn about cut cells below maxLevel
        void checkCutCells(const label maxLevel) const;

        //- No copy construct
        ibmRefineFvMesh(const ibmRefineFvMesh&) = delete;

        //- No copy assignment
        void operator=(const ibmRefineFvMesh&) = delete;


public:

    TypeName("ibmRefineFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit ibmRefineFvMesh(const IOobject& io, const bool doInit=true);


    virtual ~ibmRefineFvMesh() = default;


    // Member Functions

        const hexRef& meshCutter() const { return meshCutter_(); }

        //- Adapt the mesh. Returns true if the topology changed.
        virtual bool update();

        //- Map the fields and the refinement levels
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Write the mesh, the refinement data and the cell level
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool valid
        ) const;
};

} // End namespace Foam

#endif
