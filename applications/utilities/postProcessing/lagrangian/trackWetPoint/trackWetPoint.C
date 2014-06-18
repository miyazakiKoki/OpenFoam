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
    Info<< "end Time="<<endTime<<endl;

    //List<scalar> totLiqV(timeDirs.size(),0.0);
    //List<scalar> Vave (timeDirs.size());
    List<vector> wetPointPosition(timeDirs.size(),vector::zero);
    List<vector> wetPointVector(timeDirs.size(),vector::zero);
    List<vector> wetPointPosition1(timeDirs.size(),vector::zero);
    List<vector> wetPointVector1(timeDirs.size(),vector::zero);
    List<vector> wetPointPosition2(timeDirs.size(),vector::zero);
    List<vector> wetPointVector2(timeDirs.size(),vector::zero);
    //scalar nCell;
    //nCell = U.size();

    labelList maxIds(Pstream::nProcs(), -1);
    forAll(timeDirs, timeI)
    {

        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        //if(runTime.value() >= startTime)
        if( (runTime.value()>=startTime) && (runTime.value()<=endTime) )
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
            //List<scalar> VCell(U.size(),0.0);
            forAllConstIter(basicWetCollidingCloud, myCloud, iter)
            {
                //totLiqV[timeI] += iter().Vliq(); 
                wetPointPosition[timeI] = iter().liquidPositions()[0];
                wetPointVector[timeI] = iter().liquidPositionVectors()[0];
                wetPointPosition1[timeI] = iter().liquidPositions()[1];
                wetPointVector1[timeI] = iter().liquidPositionVectors()[1];
                wetPointPosition2[timeI] = iter().liquidPositions()[2];
                wetPointVector2[timeI] = iter().liquidPositionVectors()[2];
                //label cellI = iter().cell();
                //VCell[cellI] += iter().Vliq(); 
            }

            
        }else{
            Info<< "pass" << endl;
        }
    }




    fileName rootName(runTime.path());
    fileName fName("trackWetPoint");
    OFstream os(rootName/fName);
    Info<< "Writing to the output file"<<endl;

    forAll(wetPointPosition,i)
    {
        os<<timeDirs[i]<<","<<"liquidPosition"<<","<<wetPointPosition[i]<<","<<"liquidPositionVector"<<","<<wetPointVector[i]<<","<<"liquidPosition1"<<","<<wetPointPosition1[i]<<","<<"liquidPositionVector1"<<","<<wetPointVector1[i]<<"liquidPosition2"<<","<<wetPointPosition2[i]<<","<<"liquidPositionVector2"<<","<<wetPointVector2[i]<<endl;
       
    }

    return 0;
}




// ************************************************************************* //
