/*---------------------------------------------------------------------------*\
Application
    ibmSetGeometry

Description
    Compute and write the cut-cell immersed-boundary geometry fields
    (ibmAlpha, ibmTheta, ibmSw, ibmWallDistance, ibmNoSlipCoeff,
    ibmBodyIndex) from constant/ibmProperties, without running a solver.
    Useful to inspect the body representation, or to produce an alpha field
    that can be fed back through a body of type alphaField.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cutCellGeometry.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Compute and write the cut-cell immersed boundary geometry fields"
    );

    argList::addBoolOption
    (
        "writeFaceCentres",
        "write the face centres (surfaceVectorField Cf) for post-processing"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    cutCellGeometry ibm(mesh);
    if (!ibm.getOrDefault<bool>("writeFields", true))
    {
        ibm.writeGeometryFields();
    }

    if (args.found("writeFaceCentres"))
    {
        surfaceVectorField Cf
        (
            IOobject("Cf", runTime.timeName(), mesh),
            mesh.Cf()
        );
        Cf.write();
    }

    // Uniform pressure and velocity checks of the operators
    {
        volScalarField one
        (
            IOobject("one", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar(dimless, 1.0),
            "zeroGradient"
        );
        volVectorField g(ibm.grad(one));
        Info<< "Cut-cell gradient of a uniform field, max magnitude: "
            << gMax(mag(g.primitiveField())) << " (zero by closure)" << endl;
    }

    Info<< nl << "End\n" << endl;
    return 0;
}
