/*---------------------------------------------------------------------------*\
    cutCellIBM: cut-cell immersed boundary method for OpenFOAM

Class
    Foam::ibmScalarTransportList

Description
    The passive scalars of constant/scalarTransportProperties, each an
    ibmScalarTransport solved in turn by solve(). The list is empty (and
    solve() does nothing) when the file does not exist.

    \verbatim
    scalars
    {
        T { DT 0.1; source {...}; walls {...}; }
        c { DT 0.01; walls {...}; }
    }
    \endverbatim
    The single-scalar form "field T; DT 0.1; source {...};" at the top level
    of the file is still accepted.

\*---------------------------------------------------------------------------*/

#ifndef ibmScalarTransportList_H
#define ibmScalarTransportList_H

#include "IOdictionary.H"
#include "PtrList.H"
#include "ibmScalarTransport.H"

namespace Foam
{

class ibmScalarTransportList
:
    public PtrList<ibmScalarTransport>
{
    //- constant/scalarTransportProperties (empty if absent)
    IOdictionary dict_;

public:

    //- Construct from the geometry, the flux and the fvOptions
    ibmScalarTransportList
    (
        const cutCellGeometry& ibm,
        const surfaceScalarField& phi,
        fv::options& fvOptions
    );

    //- Solve every scalar
    void solve();
};

} // End namespace Foam

#endif
