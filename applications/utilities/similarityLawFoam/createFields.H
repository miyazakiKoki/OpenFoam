   Info<< "Reading similarityLawDict\n" << endl;

    IOdictionary similarityLawDict
    (
        IOobject
        (
            "similarityLawDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
