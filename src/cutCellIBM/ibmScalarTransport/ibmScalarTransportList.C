#include "ibmScalarTransportList.H"
#include "Time.H"

Foam::ibmScalarTransportList::ibmScalarTransportList
(
    const cutCellGeometry& ibm,
    const surfaceScalarField& phi,
    fv::options& fvOptions
)
:
    PtrList<ibmScalarTransport>(),
    dict_
    (
        IOobject
        (
            "scalarTransportProperties",
            ibm.mesh().time().constant(),
            ibm.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    )
{
    if (const dictionary* sp = dict_.findDict("scalars"))
    {
        for (const entry& e : *sp)
        {
            if (!e.isDict())
            {
                FatalIOErrorInFunction(*sp)
                    << "Entry " << e.keyword()
                    << " of scalars is not a dictionary"
                    << exit(FatalIOError);
            }
            append
            (
                new ibmScalarTransport
                (
                    ibm, phi, fvOptions, e.keyword(), e.dict()
                )
            );
        }
    }
    else if (dict_.found("field") || dict_.found("DT"))
    {
        // single-scalar form
        append
        (
            new ibmScalarTransport
            (
                ibm, phi, fvOptions, dict_.getOrDefault<word>("field", "T"), dict_
            )
        );
    }

    if (size())
    {
        Info<< "Passive scalars: " << size() << nl << endl;
    }
}


void Foam::ibmScalarTransportList::solve()
{
    for (ibmScalarTransport& s : *this)
    {
        s.solve();
    }
}
