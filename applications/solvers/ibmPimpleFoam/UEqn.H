// Momentum equation on the cut mesh (see ibmSimpleFoam/UEqn.H); the time
// derivative is weighted with the fluid volume fraction alpha

tmp<volScalarField> tnuEff(turbulence->nuEff());
const volScalarField& nuEff = tnuEff();

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(ibm.alpha(), U)
  + fvm::div(phi, U)
  - fvm::laplacian(ibm.thetaInterpolate(nuEff), U)
  + fvm::Sp(nuEff*ibm.noSlipCoeff(), U)
  + fvm::Sp(nuEff*ibm.blankCoeff(), U)
 ==
    fvOptions(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -ibm.grad(p));

    fvOptions.correct(U);
}
