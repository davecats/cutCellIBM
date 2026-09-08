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

if (ibm.moving())
{
    FatalErrorInFunction
        << "ibmSimpleFoam is a steady solver: moving bodies need ibmPimpleFoam"
        << exit(FatalError);
}

// The solid carries neither velocity nor pressure
ibm.zeroSolid(U);
ibm.zeroSolid(p);

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
setRefCell(p, simple.dict(), pRefCell, pRefValue);
ibm.checkRefCell(pRefCell);
mesh.setFluxRequired(p.name());


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
