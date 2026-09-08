Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Immersed boundary geometry (constant/ibmProperties)
cutCellGeometry ibm(mesh);

// A static solid carries neither velocity nor pressure. With moving bodies
// the initial field is kept as given: the solid interior is driven to the
// body velocity by the blanking source from the first step on, and a
// free stream u = u_b stays exact only if the solid starts at u_b too.
if (!ibm.moving())
{
    ibm.zeroSolid(U);
    ibm.zeroSolid(p);
}

Info<< "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    ibm.flux(U)
);


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
ibm.checkRefCell(pRefCell);
mesh.setFluxRequired(p.name());

// With moving bodies the reference cell may be swallowed by a body; the
// user's choice is kept while it is fluid, else the farthest fluid cell
const label pRefCellUser = pRefCell;
(void)pRefCellUser;


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

if
(
    mesh.lookupObject<IOdictionary>("turbulenceProperties")
        .get<word>("simulationType") != "laminar"
)
{
    WarningInFunction
        << "Only simulationType laminar is supported by the cut-cell immersed"
        << " boundary: the turbulence transport equations do not see the body"
        << endl;
}

#include "createFvOptions.H"

autoPtr<ibmScalarTransport> scalarPtr(ibmScalarTransport::New(ibm, phi, fvOptions));
