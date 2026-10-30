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

} // End anonymous namespace


Foam::tmp<Foam::volScalarField> Foam::cutCellGeometry::wallViscosity
(
    const volVectorField& U,
    const volScalarField& nu
) const
{
    tmp<volScalarField> tnuW(new volScalarField("ibmWallViscosity", nu));
    scalarField& nuW = tnuW.ref().primitiveFieldRef();

    const vectorField& Sw = Sw_.primitiveField();
    for (const label celli : wallCells_)
    {
        const label bodyi = cellBody_[celli];
        if
        (
            bodyi < 0
         || bodies_[bodyi].wall() != ibmBody::NOSLIP
         || bodies_[bodyi].wallFunction() != ibmBody::SPALDING
        )
        {
            continue;
        }

        // Tangential slip relative to the wall
        vector slip(U[celli]);
        if (moving_)
        {
            slip -= UbPtr_()[celli];
        }
        const vector n(Sw[celli]/max(Awall_[celli], VSMALL));
        const scalar magUt = mag(slip - (slip & n)*n);
        const scalar y = dWall_[celli];
        const scalar Rey = magUt*y/nu[celli];

        const scalar up = spaldingUPlus(Rey);
        if (up > SMALL)
        {
            nuW[celli] = max(nu[celli], nu[celli]*Rey/sqr(up));
        }
    }

    tnuW.ref().correctBoundaryConditions();
    return tnuW;
}
