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

// Turbulent flow: the wall traction takes the wall viscosity of the immersed
// wall treatment (cutCellGeometry::wallViscosity) and the transpose part of
// the viscous stress is added; laminar flow is unchanged
const bool ibmTurbulent =
    mesh.lookupObject<IOdictionary>("turbulenceProperties")
        .get<word>("simulationType") != "laminar";

if (ibmTurbulent && turbulence->type().substr(0, 3) != "ibm")
{
    WarningInFunction
        << "The turbulence model " << turbulence->type() << " does not see the"
        << " immersed walls (its equations, wall distance and wall conditions"
        << " ignore them). Use an ibm model, e.g. ibmSpalartAllmaras, with"
        << " libs (\"libcutCellIBMTurbulence.so\");" << endl;
}

#include "createFvOptions.H"

// Passive scalars of constant/scalarTransportProperties (empty if absent)
ibmScalarTransportList scalars(ibm, phi, fvOptions);
