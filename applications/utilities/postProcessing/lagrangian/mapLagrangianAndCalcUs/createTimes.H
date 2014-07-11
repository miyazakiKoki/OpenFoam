    Info<< "\nCreate databases as time" << endl;

    HashTable<string> srcOptions(args.options());
    srcOptions.erase("case");
    srcOptions.insert("case", fileName(rootDirSource/caseDirSource));

    argList argsSrc(args, srcOptions, false, false, false);

    Time runTime(Time::controlDictName, argsSrc);
    
    HashTable<string> tgtOptions(args.options());
    tgtOptions.erase("case");
    tgtOptions.insert("case", fileName(rootDirTarget/caseDirTarget));
    
    argList argsTgt(args, tgtOptions, false, false, false);

    Time runTimeTarget(Time::controlDictName, argsTgt);
