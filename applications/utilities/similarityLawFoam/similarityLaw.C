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
    vector cubePointMin
    (
        similarityLawDict.lookup("cubePointMin")
    );
    
    vector cubePointMax
    (
        similarityLawDict.lookup("cubePointMax")
    );
    
    label xDivisionNumber
    (
        readScalar(similarityLawDict.lookup("xDivisionNumber"))
    );
    
    label yDivisionNumber
    (
        readScalar(similarityLawDict.lookup("yDivisionNumber"))
    );
    
    label zDivisionNumber
    (
        readScalar(similarityLawDict.lookup("zDivisionNumber"))
    );
    
    //***********Calculating the length of the mesh defined by cubePointMin and cubePointMax************//
    scalar cubeXLength = mag(cubePointMax.x()-cubePointMin.x());
    scalar cubeYLength = mag(cubePointMax.y()-cubePointMin.y());
    scalar cubeZLength = mag(cubePointMax.z()-cubePointMin.z());
     
    
    //***********Calculating the length of the mesh devided by each division number************//
    scalar meshXLength = cubeXLength/xDivisionNumber;
    scalar meshYLength = cubeYLength/yDivisionNumber;
    scalar meshZLength = cubeZLength/zDivisionNumber;
    
    //***********making bound box for each division. This mesh is to judge whether the particle exist in the mesh*************//
    //****number of bound box needed******//
    label nbb=xDivisionNumber*yDivisionNumber*zDivisionNumber;
    
    //******making boundbox arrays*******//
    boundBox bb[nbb];
    
    //counter//
    label l=0;
    for(label k=0;k<xDivisionNumber;k++)
        {
            for(label j=0;j<zDivisionNumber;j++)
            {
                for(label i=0;i<yDivisionNumber;i++)
                {
                    l=l+1;
                    //minPoint of each boundBox//
                    scalar meshXPointMin=cubePointMin.x()+k*meshXLength;
                    scalar meshYPointMin=cubePointMin.y()+i*meshYLength;
                    scalar meshZPointMin=cubePointMin.z()+j*meshZLength;
                    
                    //maxPoint of each boundBox//
                    scalar meshXPointMax=cubePointMax.x()+(k+1)*meshXLength;
                    scalar meshYPointMax=cubePointMax.y()+(i+1)*meshYLength;
                    scalar meshZPointMax=cubePointMax.z()+(j+1)*meshZLength;
                    
                    //minPoint vector for mesh for boundBox//
                    vector meshPointMin (meshXPointMin,meshYPointMin,meshYPointMin);
                    //maxPoint vector for mesh for boundBox//
                    vector meshPointMax (meshXPointMax,meshYPointMax,meshYPointMax);
                    
                    
                    //making boundBox//
                    bb[l-1]=boundBox(meshPointMin,meshPointMax); 
                }
            }
        }
        
    
    
    
    
    
    
    
    
    
    
    


    return 0;
}


// ************************************************************************* //
