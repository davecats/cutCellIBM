#include "ibmScalarTransport.H"
#include "fvm.H"
#include "fvc.H"
#include "Time.H"

namespace Foam
{
    defineTypeNameAndDebug(ibmScalarTransport, 0);
}


Foam::ibmScalarTransport::ibmScalarTransport
(
    const cutCellGeometry& ibm,
    const surfaceScalarField& phi,
    fv::options& fvOptions,
    const word& name,
    const dictionary& dict
)
:
    ibm_(ibm),
    mesh_(ibm.mesh()),
    phi_(phi),
    fvOptions_(fvOptions),
    T_
    (
        IOobject
        (
            name,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    DT_("DT", dimViscosity, dict),
    q_
    (
        IOobject("ibmScalarSource:" + name, mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar(T_.dimensions()/dimTime, Zero),
        "zeroGradient"
    ),
    bodyFixed_(ibm.bodies().size(), false),
    bodyValue_(ibm.bodies().size(), Zero),
    wallCoeff_
    (
        IOobject("ibmScalarWallCoeff:" + name, mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar(dimless/dimArea, Zero),
        "zeroGradient"
    ),
    Twall_
    (
        IOobject("ibmScalarWallValue:" + name, mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar(T_.dimensions(), Zero),
        "zeroGradient"
    ),
    steady_(false),
    content0_(0),
    injected_(0),
    wallFlux_(0),
    sweptDefect_(0),
    alphaPrev_(ibm.alpha().primitiveField())
{
    Info<< "Passive scalar " << T_.name() << ": DT = " << DT_.value() << endl;

    // Wall condition per body: from "walls" of this scalar, else from the
    // body's own scalarWall entry, else zero gradient
    const dictionary* walls = dict.findDict("walls");
    forAll(ibm_.bodies(), bodyi)
    {
        const dictionary& bd = ibm_.bodies()[bodyi].dict();
        const dictionary* wd =
            walls ? walls->findDict(ibm_.bodies()[bodyi].name()) : nullptr;
        const word type
        (
            wd
          ? wd->get<word>("type")
          : bd.getOrDefault<word>("scalarWall", "zeroGradient")
        );
        if (type == "fixedValue")
        {
            bodyFixed_[bodyi] = true;
            bodyValue_[bodyi] =
                wd ? wd->get<scalar>("value") : bd.get<scalar>("scalarWallValue");
        }
        else if (type != "zeroGradient")
        {
            FatalIOErrorInFunction(wd ? *wd : bd)
                << "wall type of scalar " << T_.name()
                << " must be fixedValue or zeroGradient"
                << exit(FatalIOError);
        }
        Info<< "    body " << ibm_.bodies()[bodyi].name() << ": " << type;
        if (bodyFixed_[bodyi]) Info<< " " << bodyValue_[bodyi];
        Info<< endl;
    }

    // Circular (2-D) or spherical source region
    if (dict.found("source"))
    {
        const dictionary& sd = dict.subDict("source");
        const point c(sd.get<point>("centre"));
        const scalar r(sd.get<scalar>("radius"));
        const scalar rate(sd.get<scalar>("rate"));
        const Vector<label>& solD = mesh_.solutionD();
        label n = 0;
        forAll(q_, celli)
        {
            vector d(mesh_.C()[celli] - c);
            for (direction i = 0; i < 3; ++i)
            {
                if (solD[i] < 0) d[i] = 0;
            }
            if (mag(d) <= r)
            {
                q_[celli] = rate;
                ++n;
            }
        }
        reduce(n, sumOp<label>());
        Info<< "    source: rate " << rate << " in " << n << " cells around "
            << c << ", radius " << r << endl;
    }

    {
        const word scheme
        (
            mesh_.ddtScheme("ddt(" + ibm_.alpha().name() + ',' + T_.name() + ')')
        );
        steady_ = (scheme == "steadyState");
    }

    updateWall();
    content0_ = content();
    Info<< "    initial fluid content " << content0_ << nl << endl;
}


void Foam::ibmScalarTransport::updateWall()
{
    scalarField& wc = wallCoeff_.primitiveFieldRef();
    scalarField& tw = Twall_.primitiveFieldRef();
    wc = 0;
    tw = 0;

    const labelList& cellBody = ibm_.cellBody();
    const scalarField& geomCoeff = ibm_.wallCoeff().primitiveField();

    forAll(cellBody, celli)
    {
        const label bodyi = cellBody[celli];
        if (bodyi >= 0 && bodyFixed_[bodyi])
        {
            tw[celli] = bodyValue_[bodyi];
            wc[celli] = geomCoeff[celli];
        }
    }
    wallCoeff_.correctBoundaryConditions();
    Twall_.correctBoundaryConditions();
}


Foam::scalar Foam::ibmScalarTransport::content() const
{
    return gSum(ibm_.alpha().primitiveField()*mesh_.V()*T_.primitiveField());
}


void Foam::ibmScalarTransport::solve()
{
    const scalarField& V = mesh_.V();
    const scalarField& alpha = ibm_.alpha().primitiveField();
    const scalar dt = mesh_.time().deltaTValue();

    // The body may have moved: re-mask the wall, and account for the
    // difference between the geometric old fluid volume and the swept one
    // the equation uses (zero for static bodies)
    updateWall();
    if (ibm_.moving())
    {
        const scalarField& alphaSwept = ibm_.alpha().oldTime().primitiveField();
        sweptDefect_ +=
            gSum((alphaPrev_ - alphaSwept)*V*T_.primitiveField());
    }

    // Sources and sinks of constant/fvOptions listing this field, kept
    // for the budget (they are diagonal: Su in the source, Sp in the diag)
    fvScalarMatrix srcEqn(fvOptions_(T_));

    fvScalarMatrix TEqn
    (
        fvm::ddt(ibm_.alpha(), T_)
      + fvm::div(phi_, T_)
      - fvm::laplacian(DT_*ibm_.theta(), T_)
      + fvm::Sp(DT_*wallCoeff_, T_)
      + fvm::Sp(DT_*ibm_.blankCoeff(), T_)
     ==
        DT_*(wallCoeff_ + ibm_.blankCoeff())*Twall_
      + q_*ibm_.alpha()
      + srcEqn
    );

    TEqn.relax();
    fvOptions_.constrain(TEqn);
    TEqn.solve();
    fvOptions_.correct(T_);

    // Budget
    scalar sourceRate = gSum(q_.primitiveField()*alpha*V);
    if (srcEqn.hasDiag())
    {
        sourceRate += gSum(srcEqn.diag()*T_.primitiveField());
    }
    sourceRate -= gSum(srcEqn.source());
    const scalar wallRate =
        DT_.value()
       *gSum
        (
            wallCoeff_.primitiveField()*V
           *(T_.primitiveField() - Twall_.primitiveField())
        );

    scalar maxSolid = 0;
    for (const label celli : ibm_.solidCells())
    {
        maxSolid = max(maxSolid, mag(T_[celli] - Twall_[celli]));
    }
    reduce(maxSolid, maxOp<scalar>());

    const scalar C = content();

    if (steady_)
    {
        Info<< "scalar " << T_.name() << ": fluid content " << C
            << ", source rate " << sourceRate
            << ", wall flux rate " << wallRate
            << " (balance " << sourceRate - wallRate
            << "), max |T - Twall| in the solid " << maxSolid << endl;
    }
    else
    {
        injected_ += dt*sourceRate;
        wallFlux_ += dt*wallRate;
        Info<< "scalar " << T_.name() << ": fluid content " << C
            << ", injected " << injected_
            << ", through fixed-value walls " << wallFlux_
            << ", swept-volume defect " << sweptDefect_
            << ", budget error "
            << C - content0_ - injected_ + wallFlux_ + sweptDefect_
            << ", max |T - Twall| in the solid " << maxSolid << endl;
    }

    alphaPrev_ = alpha;
}
