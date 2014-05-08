{
    surfaceScalarField voidFractionf(fvc::interpolate(voidFraction));
    surfaceScalarField volumeFractionf(scalar(1)-voidFractionf);

    volScalarField rAU1(1.0/U1Eqn.A());
    surfaceScalarField rAU1f(1.0/fvc::interpolate(U1Eqn.A()));

    volVectorField HbyA1("HbyA1", U1);
    HbyA1 = rAU1*U1Eqn.H();

    surfaceScalarField phiHbyA1
    (
        "phiHbyA1",
        (fvc::interpolate(HbyA1) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU1, U1, phi1)
    );

    //adjustPhi(phiHbyA1, U1, p);

    surfaceScalarField phiHbyA2
    (
        fvc::interpolate(U2) & mesh.Sf()
    );

    surfaceScalarField phiHbyA
    (
        voidFractionf*phiHbyA1 + volumeFractionf*phiHbyA2
    );

    surfaceScalarField Dp
    (
        "Dp",
        voidFractionf*rAU1f
    );

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(Dp, p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi1 = phiHbyA1 - (fvc::snGrad(p)*mesh.magSf())*rAU1f;
//            phi1 = phiHbyA1 - pEqn.flux()/voidFractionf;
            phi2 = phiHbyA2;
            phi = voidFractionf*phi1 + volumeFractionf*phi2;
        }
    }

    U1 = HbyA1 - rAU1*fvc::grad(p);
    U1.correctBoundaryConditions();

    U = voidFraction*U1 + (scalar(1)-voidFraction)*U2;

    #include "continuityErrs.H"
}
