   Info<< "Reading setParticlePositionsDict\n" << endl;

    IOdictionary setParticlePositionsDict
    (
        IOobject
        (
            "setParticlePositionsDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
