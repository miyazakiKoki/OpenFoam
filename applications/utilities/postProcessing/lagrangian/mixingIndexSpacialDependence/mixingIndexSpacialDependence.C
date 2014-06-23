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


    List<scalar> Vave (timeDirs.size());
    List<scalar> mixingIndex(timeDirs.size(),0.0);
    scalar nCell;
    nCell = U.size();
    List<scalar> alphaR(timeDirs.size(),0.0);
    List<scalar> VcellTot(timeDirs.size(),0.0);
    List<scalar> VliqTot(timeDirs.size(),0.0);
    List<scalar> S(timeDirs.size(),0.0);
    


    labelList maxIds(Pstream::nProcs(), -1);
    forAll(timeDirs, timeI)
    {

        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        if(timeI > 1)
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
            List<scalar> VliqCell(U.size(),0.0);
            List<scalar> nPCell(U.size(),0.0);
            List<scalar> alphaI(U.size(),0.0);
            
            /*forAll(mesh.V(),celli)
            {
                Info << mesh.V()[celli] << endl;
            }*/
            
            forAllConstIter(basicWetCollidingCloud, myCloud, iter)
            {
                VliqTot[timeI] += iter().Vliq(); 
                label cellI = iter().cell();

                // number of particles in cellI
                nPCell[cellI] += 1;

                // amount of liquid in cellI
                VliqCell[cellI] += iter().Vliq(); 
                
               
            }

            


            forAll(nPCell, cellI)
            {
                if(nPCell[cellI] != 0)
                {
                    // Sum of cell volume of that particles exists
                    VcellTot[timeI] += mesh.V()[cellI];
                }
            }
            //Info<<"VcellTot["<<timeI<<"]="<<VcellTot[timeI]<<endl;

            alphaR[timeI] = VliqTot[timeI]/VcellTot[timeI] ;


            forAll(VliqCell, cellI)
            {
                if(nPCell[cellI] != 0)
                {
                    alphaI[cellI] = VliqCell[cellI]/mesh.V()[cellI] ;
                    //Info<<"alphaI["<<cellI<<"]="<<alphaI[cellI]<<endl;
                    //Info<<"mesh.V()["<<cellI<<"]="<<mesh.V()[cellI]<<endl;
                    S[timeI] += pow(alphaI[cellI]-alphaR[timeI],2)*mesh.V()[cellI];
                }
            }
            //Info<< "Value of S["<<timeI<<"]is "<<S[timeI]<<endl;

        }
        else if(timeI == 1)
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

            List<scalar> VliqCell(U.size(),0.0);
            List<scalar> nPCell(U.size(),0.0);
            List<scalar> alphaI(U.size(),0.0);
            List<scalar> alphaI0(U.size(),0.0);
            List<scalar> VliqTime0(myCloud.size(),0.0);
            
            int k=0;
            //writing Vliq by artificially making initial state to calculate S[0]------------
            //read the value of Vliq for each particle. VliqTime0 is array to rewrite the value of Vliq for Time0
            forAllConstIter(basicWetCollidingCloud, myCloud, iter)
            {   
                VliqTime0[k] = iter().Vliq();
                k+=1;
            }

            //finding maximum value of Vliq
            int i;
            scalar maxVliq = VliqTime0[0];
            for( i =1; i<myCloud.size(); i++)
            {
                if(maxVliq < VliqTime0[i])
                {
                    maxVliq = VliqTime0[i];
                }
            }
            // rewriting the value of Vliq. if Vliq>maxVliq/2, then maxVliq, else 0
            k=0;
            forAllConstIter(basicWetCollidingCloud, myCloud, iter)
            {
                if(iter().Vliq() > maxVliq/2)
                {
                    VliqTime0[k] = maxVliq;
                }else
                {
                    VliqTime0[k] = 0;
                }
                k+=1;
            }
            //Info<<"Value of VliqTime0 is"<<endl;
            //forAll(VliqTime0,i)
            //{
            //    Info<<VliqTime0[i]<<endl;
            //}

            k=0;
            List<scalar> VliqCell0(U.size(),0.0);
            forAllConstIter(basicWetCollidingCloud, myCloud, iter)
            {
                VliqTot[timeI] += iter().Vliq(); 
                label cellI = iter().cell();

                // number of particles in cellI
                nPCell[cellI] += 1;

                

                // amount of liquid in cellI
                VliqCell[cellI] += iter().Vliq(); 
                VliqCell0[cellI] += VliqTime0[k];
                k+=1;
                  
            }


            


            forAll(nPCell, cellI)
            {
                if(nPCell[cellI] != 0)
                {
                    // Sum of cell volume of that particles exists
                    VcellTot[timeI] += mesh.V()[cellI];
                }
            }
            //Info<<"VcellTot["<<timeI<<"]="<<VcellTot[timeI]<<endl;
            alphaR[timeI] = VliqTot[timeI]/VcellTot[timeI] ;
            


            forAll(VliqCell, cellI)
            {
                if(nPCell[cellI] != 0)
                {
                    alphaI[cellI] = VliqCell[cellI]/mesh.V()[cellI] ;
                    S[timeI] += pow(alphaI[cellI]-alphaR[timeI],2)*mesh.V()[cellI];
                    alphaI0[cellI] = VliqCell0[cellI]/mesh.V()[cellI] ;
                    S[0] += pow(alphaI0[cellI]-alphaR[timeI],2)*mesh.V()[cellI];
                }
            }
            

        }            
        
    }
    
    VcellTot[0] = VcellTot[1];

    forAll(S,timeI)
    {
        S[timeI] /= VcellTot[timeI];
        //Info << "S["<<timeI<<"] = "<<S[timeI]<<endl;
    }
    
    
    forAll(mixingIndex,timeI)
    {
        mixingIndex[timeI] = (S[0]-S[timeI])/S[0];
        //Info << "mixingIndex["<<timeI<<"] = "<<mixingIndex[timeI]<<endl;
    }


    fileName rootName(runTime.path());
    fileName fName("mixingIndexSpacialDependence");
    OFstream os(rootName/fName);
    Info<< "Writing to the output file"<<endl;

    forAll(mixingIndex,i)
    {
        os<<timeDirs[i]<<","<<mixingIndex[i]<<endl;               
    }

    return 0;
}




// ************************************************************************* //
