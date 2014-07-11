/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    mapLagrangeAndCalcUs

Description
    Map lagrangian data to different mesh and calculate average particle velocity
    in each cell.
    
    Syntax: mapLagrangianAndCalcUs sourceCase targetCase

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "PhaseIncompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "fixedFluxPressureFvPatchScalarField.H"

#ifdef MPPIC
    #include "basicKinematicMPPICCloud.H"
    #define basicKinematicTypeCloud basicKinematicMPPICCloud
#else
    #include "basicKinematicCollidingCloud.H"
    #define basicKinematicTypeCloud basicKinematicCollidingCloud
#endif

int main(int argc, char *argv[])
{
    argList::validArgs.append("sourceCase");
    argList::validArgs.append("targetCase");

    timeSelector::addOptions();
    #include "addRegionOption.H"

    argList::addOption
    (
        "cloudName",
        "name",
        "specify alternative cloud name. default is 'kinematicCloud'"
    );
    
    argList args(argc, argv);

    const fileName casePath = args[1];
    const fileName rootDirSource = casePath.path();
    const fileName caseDirSource = casePath.name();
    
    Info<< "Source: " << rootDirSource << " " << caseDirSource << endl;
    word sourceRegion = fvMesh::defaultRegion;
    
    const fileName targetPath = args[2];
    const fileName rootDirTarget = targetPath.path();
    const fileName caseDirTarget = targetPath.name();

    Info<< "Target: " << rootDirTarget << " " << caseDirTarget << endl;
    word targetRegion = fvMesh::defaultRegion;

    #include "createTimes.H"
    
    instantList timeDirs = timeSelector::select0(runTime, args);
    
    #include "createNamedMesh.H"
    
    Foam::fvMesh meshTarget
    (
        Foam::IOobject
        (
            regionName,
            runTimeTarget.timeName(),
            runTimeTarget,
            Foam::IOobject::MUST_READ
        )
    );

    #include "readGravitationalAcceleration.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    volVectorField UsTarget
    (
        IOobject
        (
            "UsTarget",
            runTimeTarget.timeName(),
            meshTarget,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        meshTarget,
        dimensionedVector("zeroVector", dimVelocity, vector::zero)
    );
    
    labelList npInCell
    (
        meshTarget.C().size(),
        0
    );
        
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        runTimeTarget.setTime(timeDirs[timeI], timeI);
        
        Info<< "Time = " << runTime.timeName() << endl;

        #include "createFields.H"
        
        forAllConstIter (basicKinematicTypeCloud, kinematicCloud, iter)
        {
            label cellI = -1;
            label tetFaceI = -1;
            label tetPtI = -1;
            
            meshTarget.findCellFacePt(iter().position(), cellI, tetFaceI, tetPtI);
            UsTarget[cellI] += iter().U();
            npInCell[cellI]++;
        }
    }
    
    forAll (UsTarget, celli)
    {
        if (npInCell[celli] > 0)
        {
            UsTarget[celli] /= npInCell[celli];
        }
        else
        {
            UsTarget[celli] = vector::zero;
        }
    }
    
    UsTarget.write();
 
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
