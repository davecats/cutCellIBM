/*---------------------------------------------------------------------------*\
    cutCellIBM: wall shear of the immersed no-slip walls in turbulent flow
\*---------------------------------------------------------------------------*/

#include "cutCellGeometry.H"

namespace
{

// Spalding's law, y+ = u+ + (exp(k u+) - 1 - k u+ - (k u+)^2/2 - (k u+)^3/6)/E,
// written with Re_y = u+ y+ = |u| y/nu as g(u+) = u+ + f(u+) - Re_y/u+ = 0.
// g is increasing, negative for small u+ and positive at u+ = sqrt(Re_y) (the
// viscous-sublayer root), so bisection on (0, sqrt(Re_y)] always converges.
Foam::scalar spaldingUPlus(const Foam::scalar Rey)
{
    using namespace Foam;
    const scalar kappa = 0.41, E = 9.8;
    if (Rey <= 0) return 0;

    auto g = [&](const scalar up)
    {
        const scalar ku = kappa*up;
        if (ku > 50) return GREAT;   // exp would overflow: far beyond the root
        return up + (::exp(ku) - 1 - ku - 0.5*sqr(ku) - pow3(ku)/6)/E - Rey/up;
    };

    scalar lo = 0, hi = ::sqrt(Rey);
    if (g(hi) <= 0) return hi;
    for (label it = 0; it < 100 && hi - lo > 1e-12*hi; ++it)
    {
        const scalar mid = 0.5*(lo + hi);
        if (g(mid) > 0) hi = mid; else lo = mid;
    }
    return 0.5*(lo + hi);
}

// Inverse: u+ of Spalding's law at a given y+ (g(u+) = u+ + f(u+) - y+ is
// increasing, with the root in [0, y+] since f >= 0)
Foam::scalar spaldingUPlusOfYPlus(const Foam::scalar yp)
{
    using namespace Foam;
    const scalar kappa = 0.41, E = 9.8;
    if (yp <= 0) return 0;

    auto g = [&](const scalar up)
    {
        const scalar ku = kappa*up;
        if (ku > 50) return GREAT;
        return up + (::exp(ku) - 1 - ku - 0.5*sqr(ku) - pow3(ku)/6)/E - yp;
    };

    scalar lo = 0, hi = yp;
    for (label it = 0; it < 100 && hi - lo > 1e-12*hi; ++it)
    {
        const scalar mid = 0.5*(lo + hi);
        if (g(mid) > 0) hi = mid; else lo = mid;
    }
    return 0.5*(lo + hi);
}

} // End anonymous namespace


Foam::scalar Foam::cutCellGeometry::wallFunctionUTau
(
    const volVectorField& U,
    const volScalarField& nu,
    const label celli
) const
{
    const label bodyi = cellBody_[celli];
    if
    (
        bodyi < 0
     || bodies_[bodyi].wall() != ibmBody::NOSLIP
     || bodies_[bodyi].wallFunction() != ibmBody::SPALDING
    )
    {
        return 0;
    }

    // Sample: the cell itself at dWall; a linked cell's master (the merged
    // cell's value) at dWall, the merged centroid's distance; a thin cell's
    // wall-function partner at the distance of its own centroid
    label sample = celli;
    scalar y = dWall_[celli];
    if (linkMaster_.size() == mesh_.nCells() && linkMaster_[celli] >= 0)
    {
        sample = linkMaster_[celli];
    }
    else if (wfPartner_.size() == mesh_.nCells() && wfPartner_[celli] >= 0)
    {
        sample = wfPartner_[celli];
        y = centroidWallDist_[sample];
    }

    const vector n(Sw_[celli]/max(Awall_[celli], VSMALL));
    vector slip(U[sample]);
    if (moving_)
    {
        slip -= UbPtr_()[sample];
    }
    const scalar magUt = mag(slip - (slip & n)*n);
    const scalar up = spaldingUPlus(magUt*y/nu[celli]);
    return (up > SMALL) ? magUt/up : 0;
}


Foam::tmp<Foam::volScalarField> Foam::cutCellGeometry::wallViscosity
(
    const volVectorField& U,
    const volScalarField& nu
) const
{
    tmp<volScalarField> tnuW(new volScalarField("ibmWallViscosity", nu));
    scalarField& nuW = tnuW.ref().primitiveFieldRef();

    for (const label celli : wallCells_)
    {
        const scalar uTau = wallFunctionUTau(U, nu, celli);
        if (uTau <= 0) continue;

        // nu_w |u_t|/dWall = u_tau^2 on the cell's own wall
        vector slip(U[celli]);
        if (moving_)
        {
            slip -= UbPtr_()[celli];
        }
        const vector n(Sw_[celli]/max(Awall_[celli], VSMALL));
        const scalar magUt = max(mag(slip - (slip & n)*n), VSMALL);
        nuW[celli] = max(nu[celli], sqr(uTau)*dWall_[celli]/magUt);
    }

    tnuW.ref().correctBoundaryConditions();
    return tnuW;
}


void Foam::cutCellGeometry::correctWallFaceViscosity
(
    const volVectorField& U,
    const volScalarField& nu,
    surfaceScalarField& nuF
) const
{
    if (wallFunctionAlpha_ <= 0)
    {
        return;
    }

    const vectorField& Sf = mesh_.faceAreas();
    const scalarField& deltaCoeffs =
        secondOrder_
      ? cutDeltaCoeffs_.primitiveField()
      : mesh_.nonOrthDeltaCoeffs().primitiveField();

    forAll(wfPartner_, celli)
    {
        const label n = wfPartner_[celli];
        if (n < 0) continue;
        const label facei = wfFace_[celli];

        const scalar uTau = wallFunctionUTau(U, nu, celli);
        if (uTau <= 0) continue;
        const scalar nuc = nu[celli];
        const scalar dUp =
            spaldingUPlusOfYPlus(centroidWallDist_[n]*uTau/nuc)
          - spaldingUPlusOfYPlus(centroidWallDist_[celli]*uTau/nuc);
        if (dUp <= SMALL) continue;

        const scalar nuLaw = max(nuc, uTau/(deltaCoeffs[facei]*dUp));
        const vector nw(Sw_[celli]/max(Awall_[celli], VSMALL));
        const scalar w = sqr(nw & Sf[facei]/max(mag(Sf[facei]), VSMALL));
        nuF[facei] = (1 - w)*nuF[facei] + w*theta_[facei]*nuLaw;
    }
}
