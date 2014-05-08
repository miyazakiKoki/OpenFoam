/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    particleTracks

Description
    Generates a VTK file of particle tracks for cases that were computed using
    a tracked-parcel-type cloud.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Cloud.H"
#include "IOdictionary.H"
#include "fvMesh.H"
#include "Time.H"
#include "timeSelector.H"
#include "OFstream.H"
//#include "passiveParticleCloud.H"
#include "basicWetCollidingCloud.H"
#include "writer.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"

    #include "setRootCase.H"

    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"
    #include "createFields.H"
    #include "readGravitationalAcceleration.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    fileName vtkPath(runTime.path()/"VTK");
    mkDir(vtkPath);

    Info<< "Scanning times to determine track data for cloud " << cloudName
        << nl << endl;

    List<label> iCount(U.size(),0);
    
    List<vector> particleAPosition(timeDirs.size(),vector::zero);
    List<vector> particleBPosition(timeDirs.size(),vector::zero);
    scalar particleDiameter;


   

    labelList maxIds(Pstream::nProcs(), -1);
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        Info<< "    Reading particle positions" << endl;
        //passiveParticleCloud myCloud(mesh, cloudName);

        basicWetCollidingCloud myCloud
        (
            cloudName,
            rhoInf,
            U,
            mu,
            g
        );
        label counter = 0;
        Info<< "    Read " << returnReduce(myCloud.size(), sumOp<label>())
            << " particles" << endl;

        //forAllConstIter(passiveParticleCloud, myCloud, iter)
        forAllConstIter(basicWetCollidingCloud, myCloud, iter)
        {
            particleDiameter = iter().d();
            if(counter==0)
            {
                particleAPosition[timeI]=iter().position();
            }else{
                particleBPosition[timeI]=iter().position();
            }
            counter++;
        }
    }
    
    List<scalar> particlesDistance(timeDirs.size());
    List<scalar> particleSeparationDistance(timeDirs.size());
    forAll(particlesDistance,i)
    {
        particlesDistance[i]= mag(particleAPosition[i]-particleBPosition[i]);
        //particlesDistance[i]= particleAPosition[i].x()-particleBPosition[i].x();
        particleSeparationDistance[i]=particlesDistance[i]-particleDiameter;
        
    }

    
    forAll(particleSeparationDistance,i)
    {
    Info<<"particleSeparationDistance"<<particleSeparationDistance[i]<<endl;
    }
    
    
    particleSeparationDistance[0]=particleSeparationDistance[1];
    
    fileName rootName(runTime.path());
    fileName fName("contactTime");
    OFstream os(rootName/fName);
    Info<< "Writing to the output file"<<endl;

        forAll(particleSeparationDistance,i)
        {
            if(particleSeparationDistance[i]<0)
            {
            os << "Time"<<timeDirs[i]<<endl;
            os << "particleSeparationDistance"<<particleSeparationDistance[i]<<endl;  
            }else{
            os << "No Contact"<<endl;
            }          
        }
    fileName rootNameB(runTime.path());
    fileName fNameB("particleSeparationDistance");
    OFstream osB(rootNameB/fNameB);
    Info<< "Writing to the output file"<<endl;

        forAll(particleSeparationDistance,i)
        {
            
            osB << "particleSeparationDistance"<<particleSeparationDistance[i]<<endl;  
                     
        }
    
    

    
    


/*
    fileName rootName(runTime.path());
    fileName fName("position");
    OFstream os(rootName/fName);
    Info<< "Writing to the output file"<<endl;

    forAll(particleAPosition,i)
    {
    	os <<"particleAPosition"<<particleAPosition[i]<<endl;
    }
    fileName rootNameB(runTime.path());
    fileName fNameB("positionB");
    OFstream osB(rootNameB/fNameB);
    Info<< "Writing to the output file"<<endl;

    forAll(particleBPosition,i)
    {
    	osB <<"particleBPosition"<<particleBPosition[i]<<endl;
    }
*/


    return 0;
}


// ************************************************************************* //
