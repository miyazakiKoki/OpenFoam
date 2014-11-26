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
    List<scalar> S (timeDirs.size(),0.0);
    
    List<scalar> liqVSTDEV(timeDirs.size());
    List<scalar> CV(timeDirs.size());



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

            liqVSTDEV[timeI] = Foam::pow(S[timeI]/myCloud.size(),0.5);
            
            Info<<"value of S is "<<S[timeI]<<endl;
            Info<<"value of liqVSTDEV is "<<liqVSTDEV[timeI]<<endl;
            Info<<"value of Vave is "<<Vave[timeI]<<endl;
        
            CV[timeI] = liqVSTDEV[timeI] / Vave[timeI];
            Info<<"value of CV is "<<CV[timeI]<<endl;
            
        }

    }

    CV[0] = CV[1];


    fileName rootName(runTime.path());
    fileName fName("CV%");
    OFstream os(rootName/fName);
    Info<< "Writing to the output file"<<endl;

    forAll(CV,i)
    {
        os<<timeDirs[i]<<","<<CV[i]<<endl;               
    }

    return 0;
}




// ************************************************************************* //
