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
    Info<< "start Time="<<startTime<<endl;
    List<label> iCount(U.size(),0);



    labelList maxIds(Pstream::nProcs(), -1);
    forAll(timeDirs, timeI)
    {

        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        if(runTime.value() >= startTime)
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

            //forAllConstIter(passiveParticleCloud, myCloud, iter)
            forAllConstIter(basicWetCollidingCloud, myCloud, iter)
            {


                //vector particlePosition1 = iter().position();
                //vector particleVelocity = iter().U();
                //label particleNumber= myCloud.size();
                //particlePosition[iter]=iter().position();

                /*Info <<"particle position " << particlePosition1<<endl;
                Info <<"particle velocity " << particleVelocity<<endl;
                Info <<"particle number " << particlenumber<<endl;*/

                label cellI = iter().cell();
                iCount[cellI]++;
                Uave[cellI] += iter().U();



            }
        }else{
            Info<< "pass" << endl;
        }
    }

    forAll(Uave,cellI)
    {
        if(iCount[cellI]==0)
        {
            Uave[cellI]=vector::zero;
        }else{
            Uave[cellI]/=iCount[cellI];
        }
    }

    Uave.write();

/*
    fileName rootName(runTime.path());
    fileName fName("Uave");
    OFstream os(rootName/fName);
    Info<< "Writing to the output file"<<endl;

    forAll(Uave,i)
    {
    	os <<"cellNumber"<<i<<"Uave"<<Uave[i]<<endl;
    }
*/
    return 0;
}


// ************************************************************************* //
