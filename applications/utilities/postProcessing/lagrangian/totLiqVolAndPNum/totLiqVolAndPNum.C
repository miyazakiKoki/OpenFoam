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

    List<scalar> totLiqVol(timeDirs.size(),0.0);
    List<scalar> PNumber(timeDirs.size(),0.0);
    List<scalar> DropPNumber(timeDirs.size(),0.0);
    List<scalar> DropPWeight(timeDirs.size(),0.0);
    List<scalar> DropLVol(timeDirs.size(),0.0);
    List<scalar> DropLWeight(timeDirs.size(),0.0);


   

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

        forAllConstIter(basicWetCollidingCloud, myCloud, iter)
        {
            totLiqVol[timeI] += iter().Vliq();
        }
        PNumber[timeI] = myCloud.size();
    }

    forAll(DropPNumber, i)
    {
        DropPNumber[i] = initParticleNumber - PNumber[i];
        DropPWeight[i] = DropPNumber[i] * (4/3)*3.14159*pow(particleDiameter/2,3)*particleDensity;
        DropLVol[i] = totLiqVol[0]-totLiqVol[i];
        DropLWeight[i] = DropLVol[i] * liquidDensity;
    }


    
    fileName rootName(runTime.path());
    fileName fName("totLiqVolInDrum");
    OFstream os(rootName/fName);
    Info<< "Writing to the output file"<<endl;

        forAll(timeDirs, i)
        {
            os<<timeDirs[i]<<","<<totLiqVol[i]<<endl;
                     
        }

    fileName rootName1(runTime.path());
    fileName fName1("totParticleNumberInDrum");
    OFstream os1(rootName1/fName1);
    Info<< "Writing to the output file"<<endl;

        forAll(timeDirs, i)
        {
            os1<<timeDirs[i]<<","<<PNumber[i]<<endl;
                     
        }

    fileName rootName2(runTime.path());
    fileName fName2("totParticleNumberDropped");
    OFstream os2(rootName2/fName2);
    Info<< "Writing to the output file"<<endl;

        forAll(timeDirs, i)
        {
            os2<<timeDirs[i]<<","<<DropPNumber[i]<<endl;
                     
        }

    fileName rootName3(runTime.path());
    fileName fName3("totLiquidVolumeDropped");
    OFstream os3(rootName3/fName3);
    Info<< "Writing to the output file"<<endl;

        forAll(timeDirs, i)
        {
            os3<<timeDirs[i]<<","<<DropLVol<<endl;
                     
        }

    fileName rootName4(runTime.path());
    fileName fName4("totLiquidWeightDropped");
    OFstream os4(rootName4/fName4);
    Info<< "Writing to the output file"<<endl;

        forAll(timeDirs, i)
        {
            os4<<timeDirs[i]<<","<<DropLWeight<<endl;
                     
        }

    
        

    


    return 0;
}


// ************************************************************************* //
