fvVectorMatrix U1Eqn(U1, U1.dimensions()*dimVol/dimTime);

{
    U1Eqn =
    (
        fvm::ddt(U1)
      + fvm::div(phi1, U1)
      - fvm::Sp(fvc::div(phi1), U1)
      - fvm::laplacian(nu, U1)
      - 1/(rho*voidFraction)*kinematicCloud.SU(U1)
    );

    solve(U1Eqn == -fvc::grad(p));

}
