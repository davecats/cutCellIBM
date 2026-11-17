/*---------------------------------------------------------------------------*\
Application
    ibmTruncationTC

Description
    Local truncation error of the ibmPimpleFoam operators (laminar, no
    linking) for the exact Taylor-Couette solution of
    tutorials/order/taylorCouette: Ri 0.5, Ro 1.5, centre (1.6 1.6), inner
    cylinder rotating at Omega = 1, outer at rest (hard-wired), nu from
    transportProperties. U and p exact at the value points, U = Ub in solid
    cells. Writes to the first time step, per unit FLUID volume:
      Rvisc  = (-lap(theta nu, U) + nu noSlip U + nu blank U - movingWallSource) / alpha
               (exact viscous term is zero: u_theta = A r + B/r is harmonic)
      Rconv  = div(phi, U)/alpha - (u.grad u)_exact = div(phi,U)/alpha + grad p_exact
      Rpres  = grad(p)/alpha - grad p_exact
      Rcont  = (div(phi) + S)/alpha  with phi = ibm.flux(U_exact)
      RcontFace, RcontWall: Rcont split into the error of the face fluxes
               (against 3-point Gauss over the wet segment) and that of the
               wall source (against the one closing the exact face fluxes)
    and Xvalue (value points), ibmAlpha.
    With linking (mergeAlpha > 0), -groupValue master|centroid puts one exact
    value into every cell of a group, at the fluid centroid of the master or
    at the centroid of the merged fluid volume (default own: every cell at its
    own value point), and Rcont, RcontFace, Rvisc are also written as group
    averages over the merged fluid volume (RcontGroup, RcontFaceGroup,
    RviscGroup: the residual of the merged equation per unit fluid volume;
    the faces inside a group cancel). Analysis:
    validation/taylorCouette_truncation.py (doc/secondOrderStudy.md).

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "cutCellGeometry.H"

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "groupValue", "own|master|centroid",
        "value point of linked groups (default own)"
    );
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary transport
    (
        IOobject("transportProperties", runTime.constant(), mesh,
                 IOobject::MUST_READ, IOobject::NO_WRITE)
    );
    const dimensionedScalar nu("nu", dimViscosity, transport);
    const scalar Ri = 0.5, Ro = 1.5, Om = 1.0;
    const point c0(1.6, 1.6, 0);
    const scalar A = -Om*Ri*Ri/(Ro*Ro - Ri*Ri), B = Om*Ri*Ri*Ro*Ro/(Ro*Ro - Ri*Ri);

    cutCellGeometry ibm(mesh);
    ++runTime;
    ibm.update();

    const vectorField xc
    (
        mesh.C().primitiveField()
      + (ibm.secondOrder() ? ibm.valueOffset() : ibm.centroidOffset()).primitiveField()
    );
    const scalarField& alpha = ibm.alpha().primitiveField();

    // Value point of linked groups (serial: the group is the master's index)
    const word groupValue(args.getOrDefault<word>("groupValue", "own"));
    vectorField xv(xc);
    const labelList& grp = ibm.linkGroup();
    if (ibm.linking() && groupValue != "own")
    {
        const vectorField xf
        (
            mesh.C().primitiveField() + ibm.centroidOffset().primitiveField()
        );
        const scalarField& V = mesh.V();
        vectorField sumX(mesh.nCells(), Zero);
        scalarField sumW(mesh.nCells(), Zero);
        forAll(grp, i)
        {
            if (grp[i] < 0) continue;
            sumX[grp[i]] += alpha[i]*V[i]*xf[i];
            sumW[grp[i]] += alpha[i]*V[i];
        }
        forAll(grp, i)
        {
            if (grp[i] < 0) continue;
            xv[i] =
                (groupValue == "master")
              ? xf[grp[i]]
              : sumX[grp[i]]/sumW[grp[i]];
        }
        Info<< "group value points: " << groupValue << endl;
    }

    volVectorField U
    (
        IOobject("U", runTime.timeName(), mesh),
        mesh, dimensionedVector(dimVelocity, Zero), "zeroGradient"
    );
    volScalarField p
    (
        IOobject("p", runTime.timeName(), mesh),
        mesh, dimensionedScalar(dimPressure/dimDensity, Zero), "zeroGradient"
    );
    vectorField gpEx(mesh.nCells(), Zero);
    forAll(xv, i)
    {
        const scalar X = xv[i].x() - c0.x(), Y = xv[i].y() - c0.y();
        const scalar r = Foam::sqrt(X*X + Y*Y);
        const scalar ut = A*r + B/r;
        p[i] = A*A*r*r/2 + 2*A*B*Foam::log(r) - B*B/(2*r*r);
        const scalar dpdr = ut*ut/r;
        gpEx[i] = vector(dpdr*X/r, dpdr*Y/r, 0);
        if (alpha[i] > 0) U[i] = vector(-ut*Y/r, ut*X/r, 0);
        else U[i] = ibm.wallVelocity() ? ibm.Ub()[i] : vector::zero;
    }
    U.correctBoundaryConditions();
    p.correctBoundaryConditions();

    const surfaceScalarField phi("phi", ibm.flux(U));
    const volScalarField nuv
    (
        IOobject("nu", runTime.timeName(), mesh), mesh, nu, "zeroGradient"
    );

    fvVectorMatrix visc
    (
      - fvm::laplacian(ibm.thetaInterpolate(nuv), U)
      + fvm::Sp(nuv*ibm.noSlipCoeff(), U)
      + fvm::Sp(nuv*ibm.blankCoeff(), U)
    );
    if (ibm.wallVelocity()) visc -= ibm.movingWallSource(nuv);
    tmp<fvVectorMatrix> conv(ibm.convection(phi, U));

    volVectorField rc(conv() & U);
    const volVectorField rv(visc & U), gp(ibm.grad(p));
    volScalarField cont("Rcont", fvc::div(phi));
    if (ibm.moving()) cont += ibm.massSource();

    volVectorField Rvisc("Rvisc", rv), Rconv("Rconv", rc), Rpres("Rpres", gp);
    forAll(alpha, i)
    {
        if (alpha[i] <= 0)
        {
            Rvisc[i] = Rconv[i] = Rpres[i] = Zero; cont[i] = 0; continue;
        }
        Rvisc[i] = rv[i]/alpha[i];
        Rconv[i] = rc[i]/alpha[i] + gpEx[i];
        Rpres[i] = gp[i]/alpha[i] - gpEx[i];
        cont[i] /= alpha[i];
    }
    Rvisc.write(); Rconv.write(); Rpres.write(); cont.write();

    // Continuity split: face-flux error against 3-point Gauss over the wet
    // segment of each face (contiguous, centroid Cf + wetFaceOffset), and the
    // wall source against the one that closes the exact face fluxes
    scalarField faceErrCell(mesh.nCells(), 0);
    {
        auto uEx = [&](const point& x)
        {
            const scalar X = x.x() - c0.x(), Y = x.y() - c0.y();
            const scalar r = Foam::sqrt(X*X + Y*Y), ut = A*r + B/r;
            return vector(-ut*Y/r, ut*X/r, 0);
        };
        const scalar dz = mesh.bounds().span().z();
        volScalarField faceErr("RcontFace", 0*cont), exactSum("exactFaceSum", 0*cont);
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();
        const scalar gx[3] = {-Foam::sqrt(0.6), 0, Foam::sqrt(0.6)};
        const scalar gw[3] = {5.0/18, 8.0/18, 5.0/18};
        forAll(own, f)
        {
            const scalar th = ibm.theta()[f];
            if (th <= 0) continue;
            const vector& Sf = mesh.Sf()[f];
            const vector nf(Sf/mag(Sf));
            const vector t(vector(0, 0, 1) ^ nf);
            const scalar halfL = 0.5*th*mag(Sf)/dz;
            const point xc(mesh.Cf()[f] + ibm.wetFaceOffset()[f]);
            scalar ex = 0;
            for (int k = 0; k < 3; ++k) ex += gw[k]*(uEx(xc + gx[k]*halfL*t) & Sf);
            ex *= th;
            const scalar err = phi[f] - ex;
            faceErr[own[f]] += err; faceErr[nei[f]] -= err;
            exactSum[own[f]] += ex; exactSum[nei[f]] -= ex;
        }
        const scalarField& V = mesh.V();
        volScalarField wallErr("RcontWall", 0*cont);
        forAll(alpha, i)
        {
            if (alpha[i] <= 0) { faceErr[i] = 0; wallErr[i] = 0; continue; }
            const scalar S = ibm.moving() ? ibm.massSource()[i]*V[i] : 0;
            // exact: sum_f ex + S_ideal = 0
            wallErr[i] = (S + exactSum[i])/(alpha[i]*V[i]);
            faceErr[i] /= alpha[i]*V[i];
        }
        faceErr.write(); wallErr.write();
        faceErrCell = faceErr.primitiveField();
    }
    // Residuals of the merged equations: fluid-volume averages over groups
    if (ibm.linking())
    {
        const scalarField& V = mesh.V();
        const scalarField& rf = faceErrCell;
        scalarField sW(mesh.nCells(), 0), sC(sW), sF(sW);
        vectorField sV(mesh.nCells(), Zero);
        forAll(grp, i)
        {
            if (grp[i] < 0) continue;
            const scalar w = alpha[i]*V[i];
            sW[grp[i]] += w; sC[grp[i]] += w*cont[i];
            sF[grp[i]] += w*rf[i]; sV[grp[i]] += w*Rvisc[i];
        }
        volScalarField gC("RcontGroup", cont), gF("RcontFaceGroup", cont);
        gF.primitiveFieldRef() = rf;
        volVectorField gV("RviscGroup", Rvisc);
        forAll(grp, i)
        {
            if (grp[i] < 0) continue;
            gC[i] = sC[grp[i]]/sW[grp[i]];
            gF[i] = sF[grp[i]]/sW[grp[i]];
            gV[i] = sV[grp[i]]/sW[grp[i]];
        }
        gC.write(); gF.write(); gV.write();
    }

    volVectorField Xv
    (
        IOobject("Xvalue", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh, dimensionedVector(dimLength, Zero)
    );
    Xv.primitiveFieldRef() = xv; Xv.write();
    volScalarField al
    (
        IOobject("ibmAlpha", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        ibm.alpha()
    );
    al.write();
    Info<< "secondOrder " << ibm.secondOrder() << ", written to " << runTime.timeName() << nl << "End" << endl;
    return 0;
}
