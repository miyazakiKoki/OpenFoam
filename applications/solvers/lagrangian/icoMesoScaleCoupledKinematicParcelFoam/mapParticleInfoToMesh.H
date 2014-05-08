// Mapping average particle velocity
U2.internalField() = mesh.lookupObject<volVectorField>
(
    kinematicCloud.name() + "Up"
);

//U2 = fvc::average(fvc::interpolate(U2));

U2.correctBoundaryConditions();

// Mapping void fraction
voidFraction.internalField() = scalar(1)
    - mesh.lookupObject<volScalarField>
    (
        kinematicCloud.name() + "Theta"
    )
    .internalField();

//voidFraction = fvc::average(fvc::interpolate(voidFraction));

voidFraction.correctBoundaryConditions();
