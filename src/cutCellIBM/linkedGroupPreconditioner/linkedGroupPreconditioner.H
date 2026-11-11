/*---------------------------------------------------------------------------*\
    cutCellIBM: cut-cell immersed boundary method for OpenFOAM

Class
    Foam::linkedGroupPreconditioner

Description
    Preconditioner for matrices with linked cells (cutCellGeometry::linkCells)
    whose links cross coupled faces (processor, cyclic).

    A link is a conductance c = mergePenalty x diagonal between two cells.
    The smoothers and preconditioners of OpenFOAM treat the coupled faces
    explicitly (block Jacobi over the processors): a Jacobi step leaves the
    mean of two cells coupled by c almost unchanged (factor c/(d + c)), a
    very small eigenvalue of the preconditioned matrix for every group with
    a link over a coupled face (split group). GAMG has such a mode on every
    level, as its agglomeration stays within a processor: the two clusters
    that hold the cells of the link are again coupled by c.

    This preconditioner wraps a base preconditioner (GAMG, DIC, DILU, ...)
    and splits the residual r of the cells of the split groups into its mean
    over the group and the rest, with P the average over every split group
    (the identity elsewhere):

        w = Q M'^-1 Q^T r + (1 - P) D^-1 (1 - P) r

    - Fields uniform in the split groups do not feel the conductances, so
      for them the base preconditioner M' is built on a copy of the matrix
      in which the links over coupled faces are softened to softPenalty x
      diagonal (those inside a processor stay: incomplete Cholesky resolves
      them): no stiff coupling over a coupled face is left on any level.
      Q^T hands the residual of the merged equation of a group (the sum
      over its cells) to its cells in proportion to their coefficients to
      the cells outside the group, and Q returns the mean of the result
      with the same weights: a uniform correction then solves every row
      of the group, whatever the softened conductance.
    - The rest, the differences inside the split groups, is dominated by
      the conductances and local: one Jacobi step with the diagonal D of
      the matrix.

    The two subspaces are nearly orthogonal in the energy of the matrix
    (cosine about sqrt(d/c)), so the preconditioned matrix is as well
    conditioned as that of a decomposition without split groups. The
    preconditioner is symmetric if the base is. It changes the convergence
    of the linear solver only, not the equations; without split groups
    (serial runs without cyclic links, decompositions that keep the groups
    together) it is the base preconditioner. The averages use the exchange
    of the group sums of cutCellGeometry (two small messages per
    application between the processors that share a split group).

    Usage, in system/fvSolution:
    \verbatim
    p
    {
        solver          PCG;
        preconditioner
        {
            preconditioner  linkedGroup;
            base            { preconditioner GAMG; smoother DICGaussSeidel; }
            softPenalty     1;      // default
        }
        ...
    }
    \endverbatim
    base: the name of the base preconditioner or its dictionary. For
    asymmetric matrices (PBiCGStab) the same with e.g. base DILU.

\*---------------------------------------------------------------------------*/

#ifndef linkedGroupPreconditioner_H
#define linkedGroupPreconditioner_H

#include "lduMatrix.H"

namespace Foam
{

class cutCellGeometry;

class linkedGroupPreconditioner
:
    public lduMatrix::preconditioner
{
    // Private data

        //- Controls handed to the base preconditioner
        dictionary baseDict_;

        //- Link conductance over the coupled faces, relative to the
        //  diagonal, in the matrix of the base preconditioner
        scalar softPenalty_;

        //- Copy of the matrix and of its interface coefficients with the
        //  softened links, and the solver that presents them to the base
        //  preconditioner (only with split groups)
        autoPtr<lduMatrix> softMatrixPtr_;
        FieldField<Field, scalar> softBouCoeffs_;
        FieldField<Field, scalar> softIntCoeffs_;
        autoPtr<lduMatrix::solver> softSolverPtr_;

        //- Base preconditioner
        autoPtr<lduMatrix::preconditioner> basePtr_;

        //- Cut-cell geometry (nullptr: no split groups in this matrix)
        const cutCellGeometry* ibmPtr_;

        //- Reciprocal of the number of cells of every split group, per
        //  local group slot
        List<solveScalar> rN_;

        //- Reciprocal of the diagonal of the cells of the split groups
        List<solveScalar> rD_;

        //- Weight of every cell of a split group in its merged equation:
        //  its share of the coefficients of the group to the cells outside
        //  (sum 1 over the group)
        List<solveScalar> weight_;


    // Private member functions

        //- Soften the links over the coupled faces in copies of the matrix
        //  and its interface coefficients; returns the number of faces
        label soften(const lduMatrix::solver& sol);

        //- The weights of the cells of the split groups
        void calcWeights(const lduMatrix::solver& sol);

        //- r with the sum over every split group handed to its cells by
        //  their weights (rMean, all cells), and r less its mean over the
        //  group divided by the diagonal (rest, cells of the split groups)
        void split
        (
            const solveScalarField& r,
            solveScalarField& rMean,
            List<solveScalar>& rest
        ) const;

        //- w, the base preconditioner applied to rMean: its weighted mean
        //  over every split group plus the rest less its mean
        void combine
        (
            solveScalarField& w,
            const List<solveScalar>& rest
        ) const;


public:

    TypeName("linkedGroup");

    linkedGroupPreconditioner
    (
        const lduMatrix::solver& sol,
        const dictionary& controls
    );

    virtual ~linkedGroupPreconditioner() = default;


    // Member functions

        virtual void precondition
        (
            solveScalarField& wA,
            const solveScalarField& rA,
            const direction cmpt = 0
        ) const;

        virtual void preconditionT
        (
            solveScalarField& wT,
            const solveScalarField& rT,
            const direction cmpt = 0
        ) const;

        #if OPENFOAM >= 2306
        virtual void setFinished(const solverPerformance& perf) const
        {
            basePtr_->setFinished(perf);
        }
        #endif
};

} // End namespace Foam

#endif
