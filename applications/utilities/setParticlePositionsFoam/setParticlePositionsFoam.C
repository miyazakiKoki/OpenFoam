/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    laplacianFoam

Description
    Solves a simple Laplace equation, e.g. for thermal diffusion in a solid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    Info << "\n setting particle positions \n" << endl;

    word positioningType
    (
        setParticlePositionsDict.lookup("positioningType")
    );

    if(positioningType == "autoCubic")
    {
        dictionary autoCubicPositionSubDict
        (
            setParticlePositionsDict.subDict("autoCubicPositionSubDict")
        );

        scalar particleDiameter
        (
            readScalar(autoCubicPositionSubDict.lookup("particleDiameter"))
        );

        label particleNumber
        (
            readScalar(autoCubicPositionSubDict.lookup("particleNumber"))
        );

        scalar boundaryLengthX
        (
            readScalar(autoCubicPositionSubDict.lookup("boundaryLengthX"))
        );

        scalar boundaryLengthY
        (
            readScalar(autoCubicPositionSubDict.lookup("boundaryLengthY"))
        );

        scalar boundaryLengthZ
        (
            readScalar(autoCubicPositionSubDict.lookup("boundaryLengthZ"))
        );

        
        //the center of the cube is (0 0 0)//

        //number of max particles in x direction//
        label maxParticX = boundaryLengthX/particleDiameter;
        //number of max particles in y direction//
        label maxParticY = boundaryLengthY/particleDiameter;
        //number of max particles in z direction//
        label maxParticZ = boundaryLengthZ/particleDiameter;

        //the first position of the particle//
        scalar firstPositionX = -boundaryLengthX/2;
        scalar firstPositionY = -boundaryLengthY/2;
        scalar firstPositionZ = -boundaryLengthZ/2;

        //space between the first (and the last) particle and boundary wall in x direction//
        scalar spaceX = (boundaryLengthX-maxParticX*particleDiameter)/2;
        //space between the first (and the last) particle and boundary wall in y direction//
        scalar spaceY = (boundaryLengthY-maxParticY*particleDiameter)/2;
        //space between the first (and the last) particle and boundary wall in z direction//
        scalar spaceZ = (boundaryLengthZ-maxParticZ*particleDiameter)/2;

        Info<<"number of max particles in x direction is "<<maxParticX<<endl;
        Info<<"number of max particles in y direction is "<<maxParticY<<endl;
        Info<<"number of max particles in z direction is "<<maxParticZ<<endl;


        //Info << "runTime.rootPath(): " << runTime.rootPath() << endl;
        //Info << "runTime.caseName(): " << runTime.caseName() << endl;
        //Info << "runTime.path(): " << runTime.path() << endl;

        fileName rootName(runTime.path());
        fileName fName("particlePositions");   
        OFstream os(rootName/fName);
        Info<< "Writing to the output file"<<endl;

        List<vector> veclist;
        veclist.setSize(particleNumber);
        //counter//
        label l1=0;

        for(label k=0;k<maxParticZ;k++)
        {
            for(label j=0;j<maxParticX;j++)
            {
                for(label i=0;i<maxParticY;i++)
                {
                    l1=l1+1;
                    if(l1<=particleNumber)
                    {
                        veclist[l1-1].z()=spaceZ+particleDiameter/2+particleDiameter*k+firstPositionZ;
                        veclist[l1-1].x()=spaceX+particleDiameter/2+particleDiameter*j+firstPositionX;
                        veclist[l1-1].y()=spaceY+particleDiameter/2+particleDiameter*i+firstPositionY;

                        os <<"("<<veclist[l1-1].x()<<"  "<<veclist[l1-1].y()<<"  "<<veclist[l1-1].z()<<")"<<endl;
                    }
                }
            }
        }


        Info<< "End\n" << endl;
    }
    return 0;
}


// ************************************************************************* //
