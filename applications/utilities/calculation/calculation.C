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


    Info << "\n calculating \n" << endl;

    word positioningType
    (
        calculationDict.lookup("positioningType")
    );

    if(positioningType == "autoCubic")
    {
        dictionary autoCubicPositionSubDict
        (
            calculationDict.subDict("autoCubicPositionSubDict")
        );

        scalar particleDiameter
        (
            readScalar(autoCubicPositionSubDict.lookup("particleDiameter"))
        );

        vector liquidPositionVectors
        (
            autoCubicPositionSubDict.lookup("liquidPositionVectors")
        );

        scalar stepFraction
        (
            readScalar(autoCubicPositionSubDict.lookup("stepFraction"))
        );

        vector particlePosition
        (
            autoCubicPositionSubDict.lookup("particlePosition")
        );
        
        vector omega
        (
            autoCubicPositionSubDict.lookup("omega")
        );


        vector liquidPositionVelocity;
        liquidPositionVelocity = (particleDiameter/2)*omega^liquidPositionVectors;
        vector newliquidPosition;
        newliquidPosition = liquidPositionVelocity*stepFraction+liquidPositionVectors+particlePosition;
        vector newliquidPositionVector;
        newliquidPositionVector = newliquidPosition-particlePosition;
        vector a(1,0,0);
        vector b;
        b = omega^a;
        
        

        fileName rootName(runTime.path());
        fileName fName("calculation");   
        OFstream os(rootName/fName);
        Info<< "Writing to the output file"<<endl;

        

        os <<"newliquidPosition"<<"("<<newliquidPosition.x()<<"  "<<newliquidPosition.y()<<"  "<<newliquidPosition.z()<<")"<<endl;
        os <<"newLiquidPositionVector"<<"("<<newliquidPositionVector.x()<<"  "<<newliquidPositionVector.y()<<"  "<<newliquidPositionVector.z()<<")"<<endl;
        os <<"omega^(1,0,0)="<<b<<endl;
                  


        Info<< "End\n" << endl;
    }
    return 0;
}


// ************************************************************************* //
