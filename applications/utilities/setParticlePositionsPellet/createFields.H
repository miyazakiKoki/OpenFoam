   Info<< "Reading setParticlePositionsPelletDict\n" << endl;

    IOdictionary setParticlePositionsPelletDict
    (
        IOobject
        (
            "setParticlePositionsPelletDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
