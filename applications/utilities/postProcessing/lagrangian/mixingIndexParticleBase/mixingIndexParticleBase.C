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


    Info<< "Scanning times to determine track data for cloud " << cloudName
        << nl << endl;
    Info<< "start Time="<<startTime<<endl;

    List<scalar> totLiqV(timeDirs.size(),0.0);
    List<scalar> Vave (timeDirs.size());
    List<scalar> mixingIndex(timeDirs.size(),0.0);
    List<scalar> S (timeDirs.size(),0.0);
    scalar S0;



    labelList maxIds(Pstream::nProcs(), -1);
    forAll(timeDirs, timeI)
    {

        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        
        if(timeI >= 1)
        {
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

            Info<< "    Read " << returnReduce(myCloud.size(), sumOp<label>())
            << " particles" << endl;

            Info<<"value of myCloud.size() is " << myCloud.size() << endl;
            
            forAllConstIter(basicWetCollidingCloud, myCloud, iter)
            {
                totLiqV[timeI] += iter().Vliq(); 
            }

            Vave[timeI] = totLiqV[timeI]/myCloud.size();

            forAllConstIter(basicWetCollidingCloud, myCloud, iter)
            {
                S[timeI] += pow((iter().Vliq()-Vave[timeI]),2);
            }

        }
    }

    forAll(timeDirs, timeI)
    {

        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        
        if(timeI == 1)
        {
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

            Info<< "    Read " << returnReduce(myCloud.size(), sumOp<label>())
            << " particles" << endl;

            Info<<"value of myCloud.size() is " << myCloud.size() << endl;

            Vave[timeI] = totLiqV[timeI]/myCloud.size();

            //intentionally make dummyVliq in time0 from time1 to calculate S0
            List<scalar> dummyVliq(myCloud.size(),0.0);
            int k;
            forAllConstIter(basicWetCollidingCloud, myCloud, iter)
            {   
                dummyVliq[k] = iter().Vliq();
                k+=1;
            }

            //finding maximum value of Vliq
            int i;
            scalar maxVliq = dummyVliq[0];
            for( i =1; i<myCloud.size(); i++)
            {
                if(maxVliq < dummyVliq[i])
                {
                    maxVliq = dummyVliq[i];
                }
            }
            // rewriting the value of Vliq. if Vliq>maxVliq/2, then maxVliq, else 0
            k=0;
            forAllConstIter(basicWetCollidingCloud, myCloud, iter)
            {
                if(iter().Vliq() > maxVliq/2)
                {
                    dummyVliq[k] = maxVliq;
                }else
                {
                    dummyVliq[k] = 0;
                }
                k+=1;
            }

            forAll(dummyVliq,k)
            {
                S0 += pow((dummyVliq[k]-Vave[timeI]),2);
            }

        }
    }

    S[0] = S0;

    forAll(mixingIndex,timeI)
    {
        mixingIndex[timeI] = (S[0]-S[timeI])/S[0];
    }

    fileName rootName(runTime.path());
    fileName fName("mixingIndexParticleBase");
    OFstream os(rootName/fName);
    Info<< "Writing to the output file"<<endl;

    forAll(mixingIndex,i)
    {
        os<<timeDirs[i]<<","<<mixingIndex[i]<<endl;               
    }

    return 0;
}




// ************************************************************************* //
