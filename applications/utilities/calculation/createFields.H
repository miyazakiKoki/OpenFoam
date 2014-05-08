   Info<< "Reading setParticlePositionsDict\n" << endl;

    IOdictionary calculationDict
    (
        IOobject
        (
            "calculationDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
