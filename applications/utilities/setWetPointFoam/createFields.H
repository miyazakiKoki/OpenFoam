   Info<< "Reading setWetPointDict\n" << endl;

    IOdictionary setWetPointDict
    (
        IOobject
        (
            "setWetPointDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
