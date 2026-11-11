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
#include "wallDist.H"

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

    argList::addBoolOption
    (
        "linkZone",
        "write the faceZone ibmLinkFaces of the faces linking small cut cells "
        "to their masters (serial), for decomposePar to keep every linked "
        "group on one processor: constraints { links { type "
        "preserveFaceZones; zones (ibmLinkFaces); } }"
    );

    argList::addBoolOption
    (
        "writeWallDist",
        "write the wall distance of the turbulence models (yWall), with the "
        "wallDist method of system/fvSchemes"
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

    if (args.found("writeWallDist"))
    {
        volScalarField y
        (
            IOobject("yWall", runTime.timeName(), mesh),
            wallDist::New(mesh).y()
        );
        Info<< "Wall distance: min " << gMin(y.primitiveField())
            << ", max " << gMax(y.primitiveField()) << endl;
        y.write();
    }

    if (args.found("linkZone"))
    {
        // A linked group split over processors is coupled through the
        // processor faces by the stiff link conductance, which the
        // processor-local smoothers and preconditioners treat explicitly:
        // the pressure solve then needs many more iterations
        // (doc/crossProcessorLinksHandout.md). decomposePar keeps the two
        // cells of every face of this zone together (preserveFaceZones).
        if (Pstream::parRun())
        {
            FatalErrorInFunction
                << "-linkZone is for the undecomposed case (before decomposePar)"
                << exit(FatalError);
        }
        const labelList& linkFace = ibm.linkFace();
        labelHashSet faceSet;
        forAll(linkFace, celli)
        {
            const label facei = linkFace[celli];
            if (facei >= 0 && facei < mesh.nInternalFaces())
            {
                faceSet.insert(facei);
            }
        }
        const labelList addr(faceSet.sortedToc());
        const boolList flipMap(addr.size(), false);

        const word zoneName("ibmLinkFaces");
        faceZoneMesh& zones = mesh.faceZones();
        const label zonei = zones.findZoneID(zoneName);
        if (zonei < 0)
        {
            const label n = zones.size();
            zones.setSize(n + 1);
            zones.set(n, new faceZone(zoneName, addr, flipMap, n, zones));
        }
        else
        {
            zones[zonei].resetAddressing(addr, flipMap);
        }
        zones.instance() = mesh.facesInstance();
        zones.writeOpt(IOobject::AUTO_WRITE);
        zones.write();
        Info<< "faceZone " << zoneName << ": " << addr.size()
            << " link faces written to " << zones.objectPath() << endl;
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
