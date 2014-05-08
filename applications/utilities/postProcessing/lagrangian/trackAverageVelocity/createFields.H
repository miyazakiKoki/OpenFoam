IOdictionary propsDict
(
    IOobject
    (
        "trackAverageVelocityProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED
    )
);

const word cloudName(propsDict.lookup("cloudName"));

scalar startTime(readScalar(propsDict.lookup("startTime")));

// outputMode: compositeFile, filePerTrack
//word outputmode(propsDict.lookup("outputMode"))

//label maxPositions(readLabel(propsDict.lookup("maxPositions")));

// outputFormat: raw, vtk
//word outputFormat(propsDict.lookup("outputFormat"));
word setFormat(propsDict.lookupOrDefault<word>("setFormat", "vtk"));




//**********************************************************//
// copy from icoUncoupledWetCollidingParcelFoam createFields.H


Info<< "\nReading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar rhoInfValue
(
    transportProperties.lookup("rhoInf")
);



volScalarField rhoInf
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    rhoInfValue
);


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//#include "createPhi.H"
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(U) & mesh.Sf()
);



singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);


volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    laminarTransport.nu()*rhoInfValue
);

word kinematicCloudName("kinematicCloud");
args.optionReadIfPresent("cloudName", kinematicCloudName);
/*
Info<< "Constructing kinematicCloud " << kinematicCloudName << endl;
basicWetCollidingCloud kinematicCloud
(
    kinematicCloudName,
    rhoInf,
    U,
    mu,
    g
);

IOobject Hheader
(
    "H",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
);

autoPtr<volVectorField> HPtr;

if (Hheader.headerOk())
{
    Info<< "\nReading field H\n" << endl;

    HPtr.reset(new volVectorField (Hheader, mesh));
}

IOobject HdotGradHheader
(
    "HdotGradH",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
);

autoPtr<volVectorField> HdotGradHPtr;

if (HdotGradHheader.headerOk())
{
    Info<< "Reading field HdotGradH" << endl;

    HdotGradHPtr.reset(new volVectorField(HdotGradHheader, mesh));
}

#include "createNonInertialFrameFields.H"
*/

//**********************************************************************************//
//parameters for this application//

volVectorField Uave
(
    IOobject
    (
        "Uave",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    Foam::vector(0,0,0)
);

