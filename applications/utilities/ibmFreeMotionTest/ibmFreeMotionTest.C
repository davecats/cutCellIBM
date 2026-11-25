/*---------------------------------------------------------------------------*\
    cutCellIBM: ibmFreeMotionTest

    Checks of the rigid-body integrator of freeMotion without a flow
    (component test, doc/developerNotes.md):

    A  torque-free axisymmetric top (J = diag(1, 1, 2), omega(0) = (0.1, 0,
       1)), all six directions free, explicit and implicit coupling (the
       implicit one with 30 correctors per step, converged): the angular
       momentum is conserved (to round-off), the transverse body-frame
       component of omega turns at (I3 - I1)/I1 omega3 = 1 (Euler's
       equations) and keeps its magnitude 0.1 up to the damping of implicit
       Euler, first order in dt (as the kinetic energy); the phase error is
       of higher order (below 1e-3 after one period)
    B  one free axis, isotropic body, a constant moment: omega = M t/I
       exactly (the single-axis scheme)
    C  a fixed axis (rotation (1 0 1)): its component stays zero
    D  translation in all directions under a constant force: U = F t/m
       exactly

    Prints the measured quantities and PASS/FAIL per check.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "freeMotion.H"
#include "IStringStream.H"
#include "mathematicalConstants.H"

using namespace Foam;

namespace
{
    dictionary motionDict(const string& s)
    {
        IStringStream is(s);
        return dictionary(is);
    }

    // omega in the frame of the body at rest
    vector bodyOmega(const freeMotion& fm)
    {
        return fm.orientation(0).invTransform(fm.omega());
    }

    bool check(const word& name, const bool ok, label& nFail)
    {
        Info<< "    " << name << ": " << (ok ? "PASS" : "FAIL") << nl;
        if (!ok) ++nFail;
        return ok;
    }
}


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::noBanner();
    argList args(argc, argv, false, false, false);

    const scalar pi = constant::mathematical::pi;
    label nFail = 0;

    // A: torque-free axisymmetric top
    Info<< "A: torque-free axisymmetric top, J = diag(1 1 2), omega(0) = "
        << "(0.1 0 1), T = 2 pi" << nl;
    for (const word coupling : {"explicit", "implicit"})
    {
        List<scalar> phaseErr, ampErr;
        for (const scalar dt : {0.01, 0.005, 0.0025})
        {
            freeMotion fm
            (
                motionDict
                (
                    "type free; density 1; coupling " + coupling + ";"
                    " translation (1 1 1); rotation (1 1 1);"
                    " angularVelocity (0.1 0 1); addedMass 0.5;"
                    " volume 1; inertiaTensor (1 0 0 1 0 2);"
                )
            );
            fm.setMassProperties(1, symmTensor::I);
            const vector L0(fm.globalInertia(fm.orientation(0)) & fm.omega());
            const scalar E0 = 0.5*(fm.omega() & L0);
            const label n = label(std::round(2*pi/dt));
            scalar dL = 0;
            scalar phase = 0;
            vector wbPrev(bodyOmega(fm));
            for (label k = 0; k < n; ++k)
            {
                fm.predict(dt, Zero, Zero);
                if (coupling == "implicit")
                {
                    for (label c = 0; c < 30; ++c) fm.correct(dt, Zero, Zero);
                }
                const vector L
                (
                    fm.globalInertia(fm.orientation(0)) & fm.omega()
                );
                dL = max(dL, mag(L - L0)/mag(L0));
                const vector wb(bodyOmega(fm));
                // accumulated turn of the transverse body-frame component
                phase += std::atan2
                (
                    wbPrev.x()*wb.y() - wbPrev.y()*wb.x(),
                    wbPrev.x()*wb.x() + wbPrev.y()*wb.y()
                );
                wbPrev = wb;
            }
            const vector L(fm.globalInertia(fm.orientation(0)) & fm.omega());
            const scalar E = 0.5*(fm.omega() & L);
            const vector wb(bodyOmega(fm));
            phaseErr.append(mag(phase - n*dt));
            ampErr.append
            (
                mag(Foam::sqrt(sqr(wb.x()) + sqr(wb.y())) - 0.1)
            );
            Info<< "    " << coupling << " dt " << dt
                << ": max |L - L0|/|L0| " << dL
                << ", phase " << phase << " (exact " << n*dt << ")"
                << ", body omega_3 " << wb.z()
                << ", |omega_t| " << Foam::sqrt(sqr(wb.x()) + sqr(wb.y()))
                << ", E/E0 - 1 " << E/E0 - 1 << nl;
            check("A " + coupling + " angular momentum conserved", dL < 1e-10, nFail);
        }
        const scalar p1 = Foam::log(ampErr[0]/ampErr[1])/Foam::log(2.0);
        const scalar p2 = Foam::log(ampErr[1]/ampErr[2])/Foam::log(2.0);
        Info<< "    " << coupling << " amplitude error order " << p1 << " "
            << p2 << ", phase errors " << phaseErr << nl;
        check
        (
            "A " + coupling + " amplitude first order, phase",
            p2 > 0.9 && p2 < 1.1 && max(phaseErr) < 1e-3,
            nFail
        );
    }

    // B: one free axis, constant moment
    {
        freeMotion fm
        (
            motionDict
            (
                "type free; density 2; coupling explicit; translation (0 0 0);"
                " rotation (0 0 1); volume 1; inertia 0.5;"
            )
        );
        fm.setMassProperties(1, symmTensor::I);
        const scalar dt = 0.01, M = 0.3;
        for (label k = 0; k < 100; ++k) fm.predict(dt, Zero, vector(0, 0, M));
        const scalar exact = M*1.0/(2*0.5);
        Info<< "B: one free axis, constant moment: omega_z " << fm.omega().z()
            << " (exact " << exact << "), other components "
            << fm.omega().x() << " " << fm.omega().y() << nl;
        check
        (
            "B single axis",
            mag(fm.omega().z() - exact) < 1e-12
         && fm.omega().x() == 0 && fm.omega().y() == 0,
            nFail
        );
    }

    // C: a fixed axis
    {
        freeMotion fm
        (
            motionDict
            (
                "type free; density 1; coupling implicit; translation (1 1 1);"
                " rotation (1 0 1); angularVelocity (0.2 0.5 1); volume 1;"
                " inertiaTensor (1 0.1 0 1.5 0.2 2);"
            )
        );
        fm.setMassProperties(1, symmTensor::I);
        scalar wyMax = 0;
        for (label k = 0; k < 200; ++k)
        {
            fm.predict(0.01, Zero, vector(0.1, 0.2, -0.1));
            for (label c = 0; c < 3; ++c)
            {
                fm.correct(0.01, Zero, vector(0.1, 0.2, -0.1));
            }
            wyMax = max(wyMax, mag(fm.omega().y()));
        }
        Info<< "C: rotation (1 0 1), moment about all axes: max |omega_y| "
            << wyMax << ", omega " << fm.omega() << nl;
        check("C fixed axis", wyMax == 0, nFail);
    }

    // D: translation under a constant force
    {
        freeMotion fm
        (
            motionDict
            (
                "type free; density 3; coupling explicit; volume 0.5;"
                " velocity (0 0 0);"
            )
        );
        fm.setMassProperties(0.5, symmTensor::I);
        const vector F(0.3, -0.6, 0.9);
        const scalar dt = 0.01;
        for (label k = 0; k < 100; ++k) fm.predict(dt, F, Zero);
        const vector exact(F*1.0/(3*0.5));
        Info<< "D: translation, constant force: U " << fm.velocity(0)
            << " (exact " << exact << ")" << nl;
        check("D translation", mag(fm.velocity(0) - exact) < 1e-12, nFail);
    }

    Info<< nl << (nFail ? "FAILED " : "ALL PASS ") << nFail << nl;
    return nFail ? 1 : 0;
}
