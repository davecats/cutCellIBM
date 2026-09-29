/*---------------------------------------------------------------------------*\
Application
    ibmOrderTest

Description
    Manufactured-solution tests of the cut-cell operators of the library
    (doc/secondOrderStudy.md), on the interior of a disk of radius R centred
    at "centre" (define the body with "invert true"):

      LD  -lap T = -1 (T = -(R^2 - r^2)/4),               T = 0 on the wall
      LN  -lap T = -(16 r^2 - 8 R^2) (T = r^4 - 2 R^2 r^2), dT/dn = 0
      CD  u.grad T - D lap T = q, u = omega (-y, x), T = (R^2 - r^2) x y, T = 0
      CN  same, T = (R^2 - r^2)^2 x y, dT/dn = 0
      GL  gradient operator alone, p = G.x (G = (0.7 -0.4 0)): the cut-cell
          gradient ibm.grad(p) V is compared with alpha V G; with the
          second-order treatment it is exact up to the sub-sampling error
          of alpha (the exact integral is V_geo G)
      GQ  gradient operator alone, p = X^2 + 0.5 X Y - 0.8 Y^2 (X, Y about
          the centre), compared with alpha V grad p(x_c)

    system/orderTestDict: problem, omega, D, centre, R.
    The equation is assembled with the library operators (theta-weighted
    Laplacian with the selected snGrad scheme, Dirichlet wall coefficient,
    blanking, flux through the wet part of the faces evaluated at the
    wet-face centroid), the explicit corrections are iterated to convergence,
    and the alpha-weighted L2 and Linf errors at the fluid centroids are
    printed.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cutCellGeometry.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary dict
    (
        IOobject("orderTestDict", runTime.system(), mesh,
                 IOobject::MUST_READ, IOobject::NO_WRITE)
    );
    const word problem(dict.get<word>("problem"));
    const scalar omega(dict.getOrDefault<scalar>("omega", 0));
    const scalar D(dict.getOrDefault<scalar>("D", 1));
    const point centre(dict.get<point>("centre"));
    const scalar R(dict.get<scalar>("R"));
    const bool dirichlet = (problem == "LD" || problem == "CD");
    const bool convective = (problem == "CD" || problem == "CN");

    cutCellGeometry ibm(mesh);

    if (problem == "GL" || problem == "GQ")
    {
        // Pressure given at the fluid centroids, gradient operator alone
        const vectorField xc
        (
            mesh.C().primitiveField() + ibm.centroidOffset().primitiveField()
        );
        volScalarField p
        (
            IOobject("p", runTime.timeName(), mesh),
            mesh, dimensionedScalar(dimless, Zero), "zeroGradient"
        );
        vectorField gex(mesh.nCells());
        const vector G(0.7, -0.4, 0);
        forAll(xc, celli)
        {
            const scalar X = xc[celli].x() - centre.x();
            const scalar Y = xc[celli].y() - centre.y();
            if (problem == "GL")
            {
                p[celli] = G & xc[celli];
                gex[celli] = G;
            }
            else
            {
                p[celli] = X*X + 0.5*X*Y - 0.8*Y*Y;
                gex[celli] = vector(2*X + 0.5*Y, 0.5*X - 1.6*Y, 0);
            }
        }
        p.correctBoundaryConditions();

        // grad() returns the volume integral over the cell volume V; the
        // fluid part holds alpha V
        const volVectorField g(ibm.grad(p));
        const scalarField& alpha = ibm.alpha().primitiveField();
        const scalar gRef = gMax(mag(gex));
        scalar L2cut = 0, LinfCut = 0, L2all = 0, LinfAll = 0;
        label nCut = 0, nAll = 0;
        forAll(alpha, celli)
        {
            if (alpha[celli] <= 0) continue;
            const scalar e = mag(g[celli]/alpha[celli] - gex[celli])/gRef;
            L2all += sqr(e); LinfAll = max(LinfAll, e); ++nAll;
            if (ibm.wallCells().test(celli))
            {
                L2cut += sqr(e); LinfCut = max(LinfCut, e); ++nCut;
            }
        }
        reduce(nCut, sumOp<label>()); reduce(nAll, sumOp<label>());
        L2cut = Foam::sqrt(returnReduce(L2cut, sumOp<scalar>())/max(nCut, 1));
        L2all = Foam::sqrt(returnReduce(L2all, sumOp<scalar>())/max(nAll, 1));
        reduce(LinfCut, maxOp<scalar>()); reduce(LinfAll, maxOp<scalar>());

        // Against the volume of the geometric fluid polyhedron: for a linear
        // pressure the scheme is exact if grad(p) V = V_geo G, with
        // the volume of the geometric fluid polyhedron from the same
        // operator, V_geo = [grad(x) V]_x = [grad(y) V]_y
        scalar LinfGeo = 0, LinfVol = 0, L2Geo = 0;
        label nNoLsq = 0, nGeo = 0;
        {
            volScalarField px(p), py(p);
            forAll(xc, celli)
            {
                px[celli] = xc[celli].x();
                py[celli] = xc[celli].y();
            }
            px.correctBoundaryConditions();
            py.correctBoundaryConditions();
            const volVectorField gx(ibm.grad(px)), gy(ibm.grad(py));
            const volVectorField lx(ibm.lsqGrad(px));
            forAll(alpha, celli)
            {
                if (!ibm.wallCells().test(celli)) continue;
                if (mag(lx[celli].x() - 1) > 1e-6)
                {
                    // too few fluid neighbours: no least-squares gradient,
                    // the scheme is not exact there
                    ++nNoLsq;
                    continue;
                }
                const scalar vGeo = gx[celli].x();     // V_geo/V
                LinfVol = max(LinfVol, mag(gy[celli].y() - vGeo)/vGeo);
                const scalar e = mag(g[celli] - vGeo*gex[celli])/(vGeo*gRef);
                LinfGeo = max(LinfGeo, e);
                L2Geo += sqr(e);
                ++nGeo;
            }
            reduce(LinfGeo, maxOp<scalar>());
            reduce(LinfVol, maxOp<scalar>());
            reduce(nNoLsq, sumOp<label>());
            reduce(nGeo, sumOp<label>());
            L2Geo = Foam::sqrt(returnReduce(L2Geo, sumOp<scalar>())/max(nGeo, 1));
        }

        Info<< nl << "orderTest " << problem
            << " secondOrder=" << ibm.secondOrder()
            << " nCells=" << mesh.globalData().nTotalCells()
            << " iterations=0"
            << " L2=" << L2cut << " Linf=" << LinfCut
            << " allCells: L2=" << L2all << " Linf=" << LinfAll
            << " (relative error of the cut-cell gradient, wall cells)";
        if (ibm.secondOrder())
        {
            Info<< nl << "    against V_geo grad p(x_c): L2=" << L2Geo
                << " Linf=" << LinfGeo
                << " (V_geo from grad(x) and grad(y) agree to " << LinfVol
                << "), " << nNoLsq
                << " wall cells without a least-squares gradient excluded";
        }
        Info<< nl << endl;
        Info<< "End\n" << endl;
        return 0;
    }

    volScalarField T
    (
        IOobject("T", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
        mesh
    );

    // Exact solution and source at the fluid centroids
    const vectorField xc(mesh.C().primitiveField() + ibm.centroidOffset().primitiveField());
    scalarField Tex(mesh.nCells()), q(mesh.nCells());
    forAll(xc, celli)
    {
        const scalar X = xc[celli].x() - centre.x(), Y = xc[celli].y() - centre.y();
        const scalar r2 = X*X + Y*Y;
        if (problem == "LD")      { Tex[celli] = -(R*R - r2)/4; q[celli] = -1; }
        else if (problem == "LN") { Tex[celli] = r2*r2 - 2*R*R*r2; q[celli] = -(16*r2 - 8*R*R); }
        else if (problem == "CD")
        {
            Tex[celli] = (R*R - r2)*X*Y;
            q[celli] = omega*(X*X - Y*Y)*(R*R - r2) - D*(-12*X*Y);
        }
        else if (problem == "CN")
        {
            Tex[celli] = sqr(R*R - r2)*X*Y;
            q[celli] = omega*sqr(R*R - r2)*(X*X - Y*Y) - D*X*Y*(32*r2 - 24*R*R);
        }
        else FatalErrorInFunction << "unknown problem " << problem << exit(FatalError);
    }
    // -D lap T + conv = q: for LD/LN the source of -lap T is q (with D = 1)

    const scalarField& V = mesh.V();
    const scalarField& alpha = ibm.alpha().primitiveField();
    volScalarField source
    (
        IOobject("source", runTime.timeName(), mesh),
        mesh, dimensionedScalar(T.dimensions()/dimTime, Zero), "zeroGradient"
    );
    source.primitiveFieldRef() = q;
    if (!dirichlet)
    {
        // pure Neumann: make the source compatible with the discrete fluid volume
        const scalar mean = gSum(q*alpha*V)/gSum(alpha*V);
        source.primitiveFieldRef() -= mean;
    }

    // Flux of the rotation field through the wet part of every face,
    // evaluated at the wet-face centroid (exact for the linear field)
    surfaceScalarField phi
    (
        IOobject("phi", runTime.timeName(), mesh),
        mesh, dimensionedScalar(dimVolume/dimTime, Zero)
    );
    if (convective)
    {
        auto uAt = [&](const point& x)
        {
            return vector(-omega*(x.y() - centre.y()), omega*(x.x() - centre.x()), 0);
        };
        const surfaceVectorField& Sf = mesh.Sf();
        forAll(phi, facei)
        {
            const point x(mesh.Cf()[facei] + ibm.wetFaceOffset()[facei]);
            phi[facei] = ibm.theta()[facei]*(uAt(x) & Sf[facei]);
        }
        forAll(phi.boundaryField(), patchi)
        {
            fvsPatchScalarField& pp = phi.boundaryFieldRef()[patchi];
            const fvsPatchScalarField& tp = ibm.theta().boundaryField()[patchi];
            const fvsPatchVectorField& wp = ibm.wetFaceOffset().boundaryField()[patchi];
            const fvsPatchVectorField& Sfp = Sf.boundaryField()[patchi];
            const fvsPatchVectorField& Cfp = mesh.Cf().boundaryField()[patchi];
            forAll(pp, i)
            {
                pp[i] = tp[i]*(uAt(Cfp[i] + wp[i]) & Sfp[i]);
            }
        }
    }

    const dimensionedScalar DT("D", dimViscosity, D);
    volScalarField wallCoeff(ibm.wallCoeff());
    if (!dirichlet) wallCoeff = dimensionedScalar(wallCoeff.dimensions(), Zero);

    label pRefCell = 0;
    if (!dirichlet)
    {
        forAll(alpha, celli) { if (alpha[celli] > 0) { pRefCell = celli; break; } }
    }

    scalar change = GREAT;
    label nIter = 0;
    for (; nIter < 60 && change > 1e-11; ++nIter)
    {
        volScalarField Told(T);
        fvScalarMatrix TEqn
        (
            fvm::div(phi, T)
          - fvm::laplacian(DT*ibm.theta(), T)
          + fvm::Sp(DT*wallCoeff, T)
          + fvm::Sp(DT*ibm.blankCoeff(), T)
         ==
            source*ibm.alpha()
        );
        if (!dirichlet)
        {
            TEqn.setReference(pRefCell, 0.0);
        }
        TEqn.solve();
        change = gMax(mag(T.primitiveField() - Told.primitiveField()));
        // static geometry: the explicit corrections change only through T
        if (!ibm.secondOrder()) { ++nIter; break; }
    }

    scalarField err(T.primitiveField() - Tex);
    if (!dirichlet)
    {
        err -= gSum(err*alpha*V)/gSum(alpha*V);
    }
    scalar L2 = 0, Linf = 0;
    forAll(err, celli)
    {
        if (alpha[celli] > 0)
        {
            L2 += sqr(err[celli])*alpha[celli]*V[celli];
            Linf = max(Linf, mag(err[celli]));
        }
    }
    L2 = Foam::sqrt(returnReduce(L2, sumOp<scalar>())/gSum(alpha*V));
    reduce(Linf, maxOp<scalar>());

    Info<< nl << "orderTest " << problem
        << " secondOrder=" << ibm.secondOrder()
        << " nCells=" << mesh.globalData().nTotalCells()
        << " iterations=" << nIter
        << " L2=" << L2 << " Linf=" << Linf << nl << endl;

    T.write();
    Info<< "End\n" << endl;
    return 0;
}
