/*---------------------------------------------------------------------------*\
    cutCellIBM: cut-cell immersed boundary method for OpenFOAM
\*---------------------------------------------------------------------------*/

#include "linkedGroupPreconditioner.H"
#include "cutCellGeometry.H"
#include "fvMesh.H"
#include "Vector2D.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linkedGroupPreconditioner, 0);

    lduMatrix::preconditioner::
        addsymMatrixConstructorToTable<linkedGroupPreconditioner>
        addlinkedGroupPreconditionerSymMatrixConstructorToTable_;

    lduMatrix::preconditioner::
        addasymMatrixConstructorToTable<linkedGroupPreconditioner>
        addlinkedGroupPreconditionerAsymMatrixConstructorToTable_;

    //- Presents the matrix with the softened links to the base
    //  preconditioner, which takes its matrix from a solver
    class linkedGroupBaseSolver
    :
        public lduMatrix::solver
    {
    public:

        TypeName("linkedGroupBase");

        linkedGroupBaseSolver
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& controls
        )
        :
            lduMatrix::solver
            (
                fieldName,
                matrix,
                interfaceBouCoeffs,
                interfaceIntCoeffs,
                interfaces,
                controls
            )
        {}

        virtual solverPerformance solve
        (
            scalarField&,
            const scalarField&,
            const direction
        ) const
        {
            NotImplemented;
            return solverPerformance();
        }
    };

    defineTypeNameAndDebug(linkedGroupBaseSolver, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::linkedGroupPreconditioner::soften
(
    const lduMatrix::solver& sol
)
{
    const cutCellGeometry& ibm = *ibmPtr_;
    const fvMesh& mesh = ibm.mesh();
    const boolList& boundaryLink = ibm.boundaryLink();
    const lduInterfaceFieldPtrsList& interfaces = sol.interfaces();

    // linkCells added the conductance c = mergePenalty x diagonal to the
    // coefficient of the cell across (boundary coefficient) and to the
    // diagonal (internal coefficient): take the same off both, down to
    // softPenalty x diagonal (the row sums stay)
    const scalar keep = softPenalty_/ibm.mergePenalty();
    if (keep >= 1)
    {
        return 0;
    }

    softMatrixPtr_.reset(new lduMatrix(sol.matrix()));
    auto copyCoeffs = []
    (
        const FieldField<Field, scalar>& from,
        FieldField<Field, scalar>& to
    )
    {
        to.setSize(from.size());
        forAll(from, patchi)
        {
            if (from.set(patchi))
            {
                to.set(patchi, new scalarField(from[patchi]));
            }
        }
    };
    copyCoeffs(sol.interfaceBouCoeffs(), softBouCoeffs_);
    copyCoeffs(sol.interfaceIntCoeffs(), softIntCoeffs_);
    scalarField& diag = softMatrixPtr_->diag();

    label nSoft = 0;
    forAll(mesh.boundary(), patchi)
    {
        if (!interfaces.set(patchi)) continue;
        const fvPatch& fvp = mesh.boundary()[patchi];
        const labelUList& fc = fvp.faceCells();
        scalarField& bou = softBouCoeffs_[patchi];
        scalarField& ic = softIntCoeffs_[patchi];
        const label bStart = fvp.start() - mesh.nInternalFaces();
        forAll(fc, i)
        {
            if (!boundaryLink[bStart + i]) continue;
            const scalar less = (1 - keep)*bou[i];
            bou[i] -= less;
            ic[i] -= less;
            diag[fc[i]] -= less;
            ++nSoft;
        }
    }

    softSolverPtr_.reset
    (
        new linkedGroupBaseSolver
        (
            sol.fieldName(),
            softMatrixPtr_(),
            softBouCoeffs_,
            softIntCoeffs_,
            interfaces,
            dictionary::null
        )
    );
    return nSoft;
}


void Foam::linkedGroupPreconditioner::split
(
    const solveScalarField& r,
    solveScalarField& rMean,
    List<solveScalar>& rest
) const
{
    const labelList& cells = ibmPtr_->splitGroupCells();
    const labelList& slot = ibmPtr_->groupSlot();

    List<solveScalar> sum(rN_.size(), 0);
    for (const label celli : cells)
    {
        sum[slot[celli]] += r[celli];
    }
    ibmPtr_->groupSlotSum(sum);

    rMean = r;
    rest.setSize(cells.size());
    forAll(cells, k)
    {
        const label celli = cells[k];
        const label s = slot[celli];
        rMean[celli] = weight_[k]*sum[s];
        rest[k] = rD_[k]*(r[celli] - rN_[s]*sum[s]);
    }
}


void Foam::linkedGroupPreconditioner::combine
(
    solveScalarField& w,
    const List<solveScalar>& rest
) const
{
    const labelList& cells = ibmPtr_->splitGroupCells();
    const labelList& slot = ibmPtr_->groupSlot();

    // The weighted mean of w and the sum of the rest over every group
    List<Vector2D<solveScalar>> sum
    (
        rN_.size(),
        Vector2D<solveScalar>(0, 0)
    );
    forAll(cells, k)
    {
        Vector2D<solveScalar>& s = sum[slot[cells[k]]];
        s.x() += weight_[k]*w[cells[k]];
        s.y() += rest[k];
    }
    ibmPtr_->groupSlotSum(sum);
    forAll(cells, k)
    {
        const label celli = cells[k];
        const Vector2D<solveScalar>& s = sum[slot[celli]];
        w[celli] = s.x() + rest[k] - rN_[slot[celli]]*s.y();
    }
}


void Foam::linkedGroupPreconditioner::calcWeights
(
    const lduMatrix::solver& sol
)
{
    const cutCellGeometry& ibm = *ibmPtr_;
    const fvMesh& mesh = ibm.mesh();
    const labelList& cells = ibm.splitGroupCells();
    const labelList& slot = ibm.groupSlot();
    const labelList& group = ibm.linkGroup();
    const labelList& nbrGroup = ibm.nbrLinkGroup();
    const lduMatrix& matrix = sol.matrix();
    const labelUList& own = matrix.lduAddr().lowerAddr();
    const labelUList& nei = matrix.lduAddr().upperAddr();
    const scalarField& diag = matrix.diag();
    const scalarField& upper = matrix.upper();
    const scalarField& lower = matrix.lower();
    const FieldField<Field, scalar>& bouCoeffs = sol.interfaceBouCoeffs();
    const lduInterfaceFieldPtrsList& interfaces = sol.interfaces();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    const label nInternal = mesh.nInternalFaces();

    // For every cell of a split group: e, the sum of the coefficients of
    // its row that multiply the cells of its group (the conductances
    // cancel: what is left ties the cell to the cells outside, with the
    // sign of the diagonal), and the size of its coefficients to the cells
    // outside. Per group: (sum of e, sum of the sizes)
    List<solveScalar> e(cells.size());
    List<Vector2D<solveScalar>> sum
    (
        ibm.nGroupSlots(),
        Vector2D<solveScalar>(0, 0)
    );
    forAll(cells, k)
    {
        const label celli = cells[k];
        scalar inside = diag[celli];
        scalar outside = 0;
        for (const label facei : mesh.cells()[celli])
        {
            scalar coeff = 0;
            bool same = false;
            if (facei < nInternal)
            {
                // Coefficient of the other cell in the row of this cell
                const bool owner = (own[facei] == celli);
                const label other = owner ? nei[facei] : own[facei];
                coeff = owner ? upper[facei] : lower[facei];
                same = (group[other] == group[celli]);
            }
            else
            {
                const label patchi = pbm.whichPatch(facei);
                if (!interfaces.set(patchi)) continue;
                coeff = -bouCoeffs[patchi][facei - pbm[patchi].start()];
                same = (nbrGroup[facei - nInternal] == group[celli]);
            }
            if (same)
            {
                inside += coeff;
            }
            else
            {
                outside += mag(coeff);
            }
        }
        e[k] = max((diag[celli] < 0 ? -1 : 1)*inside, scalar(0));
        sum[slot[celli]] += Vector2D<solveScalar>(e[k], outside);
    }
    ibm.groupSlotSum(sum);

    // A group without coefficients to the outside (or with e lost in the
    // round-off of the conductances) takes the plain mean
    weight_.setSize(cells.size());
    forAll(cells, k)
    {
        const label s = slot[cells[k]];
        weight_[k] =
            (sum[s].x() > 1e-6*sum[s].y() && sum[s].y() > 0)
          ? e[k]/sum[s].x()
          : rN_[s];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::linkedGroupPreconditioner::linkedGroupPreconditioner
(
    const lduMatrix::solver& sol,
    const dictionary& controls
)
:
    lduMatrix::preconditioner(sol),
    baseDict_(),
    softPenalty_(controls.getOrDefault<scalar>("softPenalty", 1)),
    softMatrixPtr_(nullptr),
    softBouCoeffs_(),
    softIntCoeffs_(),
    softSolverPtr_(nullptr),
    basePtr_(nullptr),
    ibmPtr_(nullptr),
    rN_(),
    rD_()
{
    // The base preconditioner, from its name or its dictionary
    {
        const entry& e = controls.lookupEntry("base", keyType::LITERAL);
        if (e.isDict())
        {
            baseDict_.add("preconditioner", e.dict());
        }
        else
        {
            baseDict_.add("preconditioner", word(e.stream()));
        }
    }

    // The geometry, if this is a matrix over the cells of its mesh with
    // split groups
    const lduMatrix& matrix = sol.matrix();
    const fvMesh* meshPtr = dynamic_cast<const fvMesh*>(&matrix.mesh());
    if (meshPtr && matrix.diag().size() == meshPtr->nCells())
    {
        const cutCellGeometry* ibmPtr =
            meshPtr->cfindObject<cutCellGeometry>("ibmProperties");
        if (ibmPtr && ibmPtr->nSplitGroups() > 0)
        {
            ibmPtr_ = ibmPtr;
        }
    }

    if (!ibmPtr_)
    {
        basePtr_ = lduMatrix::preconditioner::New(sol, baseDict_);
        return;
    }

    const label nSoft = returnReduce(soften(sol), sumOp<label>());
    basePtr_ = lduMatrix::preconditioner::New
    (
        softSolverPtr_ ? softSolverPtr_() : sol,
        baseDict_
    );

    // Number of cells of every split group and the diagonal of its cells
    const labelList& cells = ibmPtr_->splitGroupCells();
    const labelList& slot = ibmPtr_->groupSlot();
    const scalarField& diag = matrix.diag();
    List<solveScalar> n(ibmPtr_->nGroupSlots(), 0);
    rD_.setSize(cells.size());
    forAll(cells, k)
    {
        n[slot[cells[k]]] += 1;
        rD_[k] = 1/diag[cells[k]];
    }
    ibmPtr_->groupSlotSum(n);
    rN_.setSize(n.size());
    rN_ = 0;
    for (const label celli : cells)
    {
        rN_[slot[celli]] = 1/n[slot[celli]];
    }
    calcWeights(sol);

    if (debug)
    {
        Info<< "linkedGroup: " << sol.fieldName() << ", "
            << ibmPtr_->nSplitGroups() << " split groups, " << nSoft
            << " coupled link faces softened to " << softPenalty_
            << " x diagonal" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::linkedGroupPreconditioner::precondition
(
    solveScalarField& wA,
    const solveScalarField& rA,
    const direction cmpt
) const
{
    if (!ibmPtr_)
    {
        basePtr_->precondition(wA, rA, cmpt);
        return;
    }

    solveScalarField rMean;
    List<solveScalar> rest;
    split(rA, rMean, rest);
    basePtr_->precondition(wA, rMean, cmpt);
    combine(wA, rest);
}


void Foam::linkedGroupPreconditioner::preconditionT
(
    solveScalarField& wT,
    const solveScalarField& rT,
    const direction cmpt
) const
{
    if (!ibmPtr_)
    {
        basePtr_->preconditionT(wT, rT, cmpt);
        return;
    }

    solveScalarField rMean;
    List<solveScalar> rest;
    split(rT, rMean, rest);
    basePtr_->preconditionT(wT, rMean, cmpt);
    combine(wT, rest);
}


// ************************************************************************* //
