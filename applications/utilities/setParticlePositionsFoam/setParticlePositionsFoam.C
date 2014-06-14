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
double GetRandom(double min, double max);

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
        fileName fName("autoCubicParticlePositions");   
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


    //-----------------------------------------detail Cubic------------------------------------//
    // make cubic with minimum point and max point

    if(positioningType == "detailCubic")
    {

        dictionary detailCubicSubDict
        (
            setParticlePositionsDict.subDict("detailCubicSubDict")
        );

        scalar particleDiameter
        (
            readScalar(detailCubicSubDict.lookup("particleDiameter"))
        );

        label particleNumber
        (
            readScalar(detailCubicSubDict.lookup("particleNumber"))
        );

        vector boundaryPointMin
        (
            detailCubicSubDict.lookup("boundaryPointMin")
        );

        vector boundaryPointMax
        (
            detailCubicSubDict.lookup("boundaryPointMax")
        );

        vector pileDirection
        (
            detailCubicSubDict.lookup("pileDirection")
        );

        

        scalar boundaryLengthX = boundaryPointMax.x()-boundaryPointMin.x();
        scalar boundaryLengthY = boundaryPointMax.y()-boundaryPointMin.y();
        scalar boundaryLengthZ = boundaryPointMax.z()-boundaryPointMin.z();


        //number of max particles in x direction//
        label maxParticX = boundaryLengthX/particleDiameter;
        //number of max particles in y direction//
        label maxParticY = boundaryLengthY/particleDiameter;
        //number of max particles in z direction//
        label maxParticZ = boundaryLengthZ/particleDiameter;

        //space between the first (and the last) particle and boundary wall in x direction//
        scalar spaceX = (boundaryLengthX-maxParticX*particleDiameter)/2;
        //space between the first (and the last) particle and boundary wall in y direction//
        scalar spaceY = (boundaryLengthY-maxParticY*particleDiameter)/2;
        //space between the first (and the last) particle and boundary wall in z direction//
        scalar spaceZ = (boundaryLengthZ-maxParticZ*particleDiameter)/2;
 
        //the first position of the particle//
        scalar firstPositionX = boundaryPointMin.x()+spaceX+particleDiameter/2;
        scalar firstPositionY = boundaryPointMin.y()+spaceY+particleDiameter/2;
        scalar firstPositionZ = boundaryPointMin.z()+spaceZ+particleDiameter/2;



        Info<<"number of max particles in x direction is "<<maxParticX<<endl;
        Info<<"number of max particles in y direction is "<<maxParticY<<endl;
        Info<<"number of max particles in z direction is "<<maxParticZ<<endl;


        //Info << "runTime.rootPath(): " << runTime.rootPath() << endl;
        //Info << "runTime.caseName(): " << runTime.caseName() << endl;
        //Info << "runTime.path(): " << runTime.path() << endl;

        fileName rootName(runTime.path());
        fileName fName("detailCubicParticlePositions");   
        OFstream os(rootName/fName);
        Info<< "Writing to the output file"<<endl;

        vector x  (1,0,0);
        vector y  (0,1,0);
        vector z  (0,0,1);


        List<vector> veclist;
        veclist.setSize(particleNumber);

        if(pileDirection == z)
        {
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
                            veclist[l1-1].z()=particleDiameter*k+firstPositionZ;
                            veclist[l1-1].x()=particleDiameter*j+firstPositionX;
                            veclist[l1-1].y()=particleDiameter*i+firstPositionY;
    
                            os <<"("<<veclist[l1-1].x()<<"  "<<veclist[l1-1].y()<<"  "<<veclist[l1-1].z()<<")"<<endl;
                        }
                    }
                }
            }
        }
        else if(pileDirection == x)
        {
            //counter//
            label l1=0;

            for(label j=0;j<maxParticX;j++)
            {
                for(label k=0;k<maxParticZ;k++)
                {
                    for(label i=0;i<maxParticY;i++)
                    {
                        l1=l1+1;
                        if(l1<=particleNumber)
                        {
                            veclist[l1-1].z()=particleDiameter*k+firstPositionZ;
                            veclist[l1-1].x()=particleDiameter*j+firstPositionX;
                            veclist[l1-1].y()=particleDiameter*i+firstPositionY;
    
                            os <<"("<<veclist[l1-1].x()<<"  "<<veclist[l1-1].y()<<"  "<<veclist[l1-1].z()<<")"<<endl;
                        }
                    }
                }
            }
        }
        else if(pileDirection == y)
        {
            //counter//
            label l1=0;

            for(label i=0;i<maxParticY;i++)
            {
                for(label k=0;k<maxParticZ;k++)
                {
                    for(label j=0;j<maxParticX;j++)
                    {
                        l1=l1+1;
                        if(l1<=particleNumber)
                        {
                            veclist[l1-1].z()=particleDiameter*k+firstPositionZ;
                            veclist[l1-1].x()=particleDiameter*j+firstPositionX;
                            veclist[l1-1].y()=particleDiameter*i+firstPositionY;
    
                            os <<"("<<veclist[l1-1].x()<<"  "<<veclist[l1-1].y()<<"  "<<veclist[l1-1].z()<<")"<<endl;
                        }
                    }
                }
            }
        }



        Info<< "End\n" << endl;
    }

    //-----------------------------------------random cylinder-----------------------------
    // make cylinder shaped. particles are allocated randomly

    if(positioningType == "randomCylinder")
    {   
        dictionary randomCylinderSubDict
        (
            setParticlePositionsDict.subDict("randomCylinderSubDict")
        );

        scalar particleDiameter
        (
            readScalar(randomCylinderSubDict.lookup("particleDiameter"))
        );

        label particleNumber
        (
            readLabel(randomCylinderSubDict.lookup("particleNumber"))
        );

        scalar radius
        (
            readScalar(randomCylinderSubDict.lookup("cylinderRadius"))
        );

        scalar length
        (
            readScalar(randomCylinderSubDict.lookup("cylinderLength"))
        );

        vector origin
        (
            randomCylinderSubDict.lookup("origin")
        );

        label trialNum
        (
            readLabel(randomCylinderSubDict.lookup("trialNumber"))
        );

        label densityR
        (
            readLabel(randomCylinderSubDict.lookup("densityRadius"))
        );

        label densityT
        (
            readLabel(randomCylinderSubDict.lookup("densityTheta"))
        );

        label densityL
        (
            readLabel(randomCylinderSubDict.lookup("densityLength"))
        );

        vector axisDirection
        (
            randomCylinderSubDict.lookup("axisDirection")
        );


        // Min and Max of each variable were multiplied by densityR,T,L respectively to increase the resolution, since random function I made can only output integer.
        // by doing so, it can output decimal value.
        scalar thetaMin = 0*densityT;
        scalar thetaMax = 360 * densityT;
        scalar rMin = 0*densityR;
        scalar rMax = (radius-particleDiameter/2) * densityR;
        scalar lMin = -(length-particleDiameter)*densityL/2;
        scalar lMax = (length-particleDiameter)*densityL/2;

        //dummy list to store particles in cylindrical coordinates x is r, y is theta, z is length

        List<vector> dummyCylindic(trialNum);

        //dummy list to store particles in rectangula coordinates

        List<vector> dummyRectangula(trialNum);

        //List for particles settled in cylindrical coordinates

        List<vector> cylindic(particleNumber);

        //List for particles transformed into rectangular coordinates

        List<vector> rectangula(particleNumber);


        vector x (1,0,0);
        vector y (0,1,0);
        vector z (0,0,1);

        if(y == axisDirection)
        {


            // the first particle is stored in cylindic and rectangula list

            cylindic[0].x() = GetRandom(rMin,rMax)/densityR;
            cylindic[0].y() = GetRandom(thetaMin,thetaMax)/densityT;
            cylindic[0].z() = GetRandom(lMin,lMax)/densityL;

            Info << "r is" << cylindic[0].x() << endl;
            Info << "theta is" << cylindic[0].y() << endl;
            Info << "length is " << cylindic[0].z() << endl;


            // cylindic[0] is transormed into rectangular coordinates
            //inside sin is written in radian

            rectangula[0].x() = cylindic[0].x()*Foam::sin(constant::mathematical::pi*cylindic[0].y()/180);
            rectangula[0].y() = cylindic[0].z();
            rectangula[0].z() = cylindic[0].x()*Foam::cos(constant::mathematical::pi*cylindic[0].y()/180);

            Info << "x is" << rectangula[0].x() << endl;
            Info << "y is" << rectangula[0].y() << endl;
            Info << "z is " << rectangula[0].z() << endl;
        
            label l = 0;

            forAll(dummyCylindic,i)
            {

                if( l< particleNumber-1 )
                {
                    dummyCylindic[i].x() = GetRandom(rMin,rMax)/densityR;
                    dummyCylindic[i].y() = GetRandom(thetaMin,thetaMax)/densityT;
                    dummyCylindic[i].z() = GetRandom(lMin,lMax)/densityL;
                    Info << " dummyCylindic[ "<<i<<"].x()" << dummyCylindic[i].x() << endl;
                    Info << " dummyCylindic[ "<<i<<"].y()" << dummyCylindic[i].y() << endl;
                    Info << " dummyCylindic[ "<<i<<"].z()" << dummyCylindic[i].z() << endl;

                    dummyRectangula[i].x() = dummyCylindic[i].x()*Foam::sin(constant::mathematical::pi*dummyCylindic[i].y()/180);
                    dummyRectangula[i].y() = dummyCylindic[i].z();
                    dummyRectangula[i].z() = dummyCylindic[i].x()*Foam::cos(constant::mathematical::pi*dummyCylindic[i].y()/180);
    
                    label counter = 0;
                    label t;
                    for(t=0; t<=l; t++)
                    {
                        scalar S;
                        S = mag(dummyRectangula[i]-rectangula[t]);
                        if(S>particleDiameter)
                        {
                            counter += 1;   
                        }else{
                            counter += 0;
                        }
                    }
             
                    if(counter == l+1)
                    {
                        l = l+1;
                        rectangula[l] = dummyRectangula[i];
                    }
                }
            } 
        

            fileName rootName(runTime.path());
            fileName fName("yAxisRandomCylinder");   
            OFstream os(rootName/fName);
            Info<< "Writing to the output file"<<endl;

            label k;
            for(k=0; k<=l; k++)
            {
                rectangula[k] += origin ;
            }

            for(k=0; k<=l; k++)
            {
                os <<"("<<rectangula[k].x()<<"  "<<rectangula[k].y()<<"  "<<rectangula[k].z()<<")"<<endl;
            }
        
            Info<< "End\n" << endl;
        }
        else if(x == axisDirection)
        {

            // the first particle is stored in cylindic and rectangula list

            cylindic[0].x() = GetRandom(rMin,rMax)/densityR;
            cylindic[0].y() = GetRandom(thetaMin,thetaMax)/densityT;
            cylindic[0].z() = GetRandom(lMin,lMax)/densityL;

            Info << "r is" << cylindic[0].x() << endl;
            Info << "theta is" << cylindic[0].y() << endl;
            Info << "length is " << cylindic[0].z() << endl;


            // cylindic[0] is transormed into rectangular coordinates
            //inside sin is written in radian

            rectangula[0].x() = cylindic[0].z();
            rectangula[0].y() = cylindic[0].x()*Foam::sin(constant::mathematical::pi*cylindic[0].y()/180);
            rectangula[0].z() = cylindic[0].x()*Foam::cos(constant::mathematical::pi*cylindic[0].y()/180);

            Info << "x is" << rectangula[0].x() << endl;
            Info << "y is" << rectangula[0].y() << endl;
            Info << "z is " << rectangula[0].z() << endl;
        
            label l = 0;

            forAll(dummyCylindic,i)
            {

                if( l< particleNumber-1 )
                {
                    dummyCylindic[i].x() = GetRandom(rMin,rMax)/densityR;
                    dummyCylindic[i].y() = GetRandom(thetaMin,thetaMax)/densityT;
                    dummyCylindic[i].z() = GetRandom(lMin,lMax)/densityL;
                    Info << " dummyCylindic[ "<<i<<"].x()" << dummyCylindic[i].x() << endl;
                    Info << " dummyCylindic[ "<<i<<"].y()" << dummyCylindic[i].y() << endl;
                    Info << " dummyCylindic[ "<<i<<"].z()" << dummyCylindic[i].z() << endl;



                    dummyRectangula[i].x() = dummyCylindic[i].z();
                    dummyRectangula[i].y() = dummyCylindic[i].x()*Foam::sin(constant::mathematical::pi*dummyCylindic[i].y()/180);
                    dummyRectangula[i].z() = dummyCylindic[i].x()*Foam::cos(constant::mathematical::pi*dummyCylindic[i].y()/180);
    
                    label counter = 0;
                    label t;
                    for(t=0; t<=l; t++)
                    {
                        scalar S;
                        S = mag(dummyRectangula[i]-rectangula[t]);
                        if(S>particleDiameter)
                        {
                            counter += 1;   
                        }else{
                            counter += 0;
                        }
                    }
             
                    if(counter == l+1)
                    {
                        l = l+1;
                        rectangula[l] = dummyRectangula[i];
                    }
                }
            } 
        

            fileName rootName(runTime.path());
            fileName fName("xAxisRandomCylinder");   
            OFstream os(rootName/fName);
            Info<< "Writing to the output file"<<endl;

            label k;
            for(k=0; k<=l; k++)
            {
                rectangula[k] += origin ;
            }

            for(k=0; k<=l; k++)
            {
                os <<"("<<rectangula[k].x()<<"  "<<rectangula[k].y()<<"  "<<rectangula[k].z()<<")"<<endl;
            }
        
            Info<< "End\n" << endl;

        }
        else if(z == axisDirection)
        {

            // the first particle is stored in cylindic and rectangula list

            cylindic[0].x() = GetRandom(rMin,rMax)/densityR;
            cylindic[0].y() = GetRandom(thetaMin,thetaMax)/densityT;
            cylindic[0].z() = GetRandom(lMin,lMax)/densityL;

            Info << "r is" << cylindic[0].x() << endl;
            Info << "theta is" << cylindic[0].y() << endl;
            Info << "length is " << cylindic[0].z() << endl;


            // cylindic[0] is transormed into rectangular coordinates
            //inside sin is written in radian

            rectangula[0].x() = cylindic[0].x()*Foam::sin(constant::mathematical::pi*cylindic[0].y()/180);
            rectangula[0].y() = cylindic[0].x()*Foam::cos(constant::mathematical::pi*cylindic[0].y()/180);
            rectangula[0].z() = cylindic[0].z();

            Info << "x is" << rectangula[0].x() << endl;
            Info << "y is" << rectangula[0].y() << endl;
            Info << "z is " << rectangula[0].z() << endl;
        
            label l = 0;

            forAll(dummyCylindic,i)
            {

                if( l< particleNumber-1 )
                {
                    dummyCylindic[i].x() = GetRandom(rMin,rMax)/densityR;
                    dummyCylindic[i].y() = GetRandom(thetaMin,thetaMax)/densityT;
                    dummyCylindic[i].z() = GetRandom(lMin,lMax)/densityL;
                    Info << " dummyCylindic[ "<<i<<"].x()" << dummyCylindic[i].x() << endl;
                    Info << " dummyCylindic[ "<<i<<"].y()" << dummyCylindic[i].y() << endl;
                    Info << " dummyCylindic[ "<<i<<"].z()" << dummyCylindic[i].z() << endl;

                    dummyRectangula[i].x() = dummyCylindic[i].x()*Foam::sin(constant::mathematical::pi*dummyCylindic[i].y()/180);
                    dummyRectangula[i].y() = dummyCylindic[i].x()*Foam::cos(constant::mathematical::pi*dummyCylindic[i].y()/180);
                    dummyRectangula[i].z() = dummyCylindic[i].z();
    
                    label counter = 0;
                    label t;
                    for(t=0; t<=l; t++)
                    {
                        scalar S;
                        S = mag(dummyRectangula[i]-rectangula[t]);
                        if(S>particleDiameter)
                        {
                            counter += 1;   
                        }else{
                            counter += 0;
                        }
                    }
             
                    if(counter == l+1)
                    {
                        l = l+1;
                        rectangula[l] = dummyRectangula[i];
                    }
                }
            } 
        

            fileName rootName(runTime.path());
            fileName fName("zAxisRandomCylinder");   
            OFstream os(rootName/fName);
            Info<< "Writing to the output file"<<endl;

            label k;
            for(k=0; k<=l; k++)
            {
                rectangula[k] += origin ;
            }

            for(k=0; k<=l; k++)
            {
                os <<"("<<rectangula[k].x()<<"  "<<rectangula[k].y()<<"  "<<rectangula[k].z()<<")"<<endl;
            }
        
            Info<< "End\n" << endl;

        }
    }

    if(positioningType == "packingCubic")
    {

        dictionary packingCubicSubDict
        (
            setParticlePositionsDict.subDict("packingCubicSubDict")
        );

        scalar particleDiameter
        (
            readScalar(packingCubicSubDict.lookup("particleDiameter"))
        );

        label particleNumber
        (
            readScalar(packingCubicSubDict.lookup("particleNumber"))
        );

        vector boundaryPointMin
        (
            packingCubicSubDict.lookup("boundaryPointMin")
        );

        vector boundaryPointMax
        (
            packingCubicSubDict.lookup("boundaryPointMax")
        );

        vector pileDirection
        (
            packingCubicSubDict.lookup("pileDirection")
        );

        

        scalar boundaryLengthX = boundaryPointMax.x()-boundaryPointMin.x();
        scalar boundaryLengthY = boundaryPointMax.y()-boundaryPointMin.y();
        scalar boundaryLengthZ = boundaryPointMax.z()-boundaryPointMin.z();

        vector x  (1,0,0);
        vector y  (0,1,0);
        vector z  (0,0,1);

        List<vector> veclist;
        veclist.setSize(particleNumber);

        if(pileDirection == z)
        {


            //number of max particles in x direction in odd row and even row when z axis is even//
            label maxParticXEvenWhenZEven = boundaryLengthX/particleDiameter;
            label maxParticXOddWhenZEven  = (boundaryLengthX-particleDiameter/2)/particleDiameter;
            //number of max particles in x direction in odd row and even row when z axis is odd//
            label maxParticXEvenWhenZOdd  = (boundaryLengthX-particleDiameter/2)/particleDiameter;
            label maxParticXOddWhenZOdd   = boundaryLengthX/particleDiameter;
        
            //number of max particles in y direction when NumZrow is Even and Odd//
            //label NumYRowEven = (boundaryLengthY/Foam::sin(constant::mathematical::pi/3))/particleDiameter;
            scalar yL = particleDiameter * Foam::sin(constant::mathematical::pi/3);
            label NumYRowEven = (boundaryLengthY - particleDiameter)/yL + 1;
            label NumYRowOdd  = (boundaryLengthY-particleDiameter/2*Foam::tan(constant::mathematical::pi/6)-particleDiameter)/yL +1;
            //label NumYRowOdd  = ((boundaryLengthY-particleDiameter/2)/Foam::sin(constant::mathematical::pi/3))/particleDiameter;
        
            //number of max particles in z direction//
            label NumZRow = (boundaryLengthZ/Foam::sin(constant::mathematical::pi/3))/particleDiameter;

            //Particle number difference between maxParticleXOdd and maxParticXeven
            label NdiffX = maxParticXEvenWhenZEven-maxParticXOddWhenZEven;
            label NdiffY = NumYRowEven-NumYRowOdd;

            //space between the first (and the last) particle and boundary wall in x direction//
            scalar spaceX = (boundaryLengthX-(maxParticXEvenWhenZEven*particleDiameter+(-NdiffX+1)*particleDiameter/2))/2;
            //space between the first (and the last) particle and boundary wall in y direction//
            //scalar spaceY = Foam::sin(constant::mathematical::pi/3)*(boundaryLengthY/Foam::sin(constant::mathematical::pi/3)-NumYRow*particleDiameter)/2+NdiffY*((particleDiameter/2)/Foam::sin(constant::mathematical::pi/6));
            scalar spaceY = (boundaryLengthY - ((NumYRowEven-1)*yL+particleDiameter+(-NdiffY+1)*particleDiameter/2*Foam::tan(constant::mathematical::pi/6)))/2;
            //space between the first (and the last) particle and boundary wall in z direction//
            scalar spaceZ = Foam::sin(constant::mathematical::pi/3)*(boundaryLengthZ/Foam::sin(constant::mathematical::pi/3)-NumZRow*particleDiameter)/2;
 
            //the first position of the particle//
            //the first Position of X and y when z axis is even
            scalar firstPositionXEven = boundaryPointMin.x()+spaceX+particleDiameter/2;
            scalar firstPositionYEven = boundaryPointMin.y()+spaceY+particleDiameter/2;
            //the first Position of X and y when z axis is odd
            scalar firstPositionXOdd = firstPositionXEven + particleDiameter*Foam::cos(constant::mathematical::pi/3);
            scalar firstPositionYOdd = firstPositionYEven + particleDiameter/2*Foam::tan(constant::mathematical::pi/6);
            scalar firstPositionZ = boundaryPointMin.z()+spaceZ+particleDiameter/2;



            Info<<"number of max particles in x direction (Y even Z even) is " << maxParticXEvenWhenZEven<<endl;
            Info<<"number of max particles in x direction (Y odd Z even) is "<<maxParticXOddWhenZEven<<endl;
            Info<<"number of rows in y direction (Even) is "<<NumYRowEven<<endl;
            Info<<"number of rows in y direction (Odd) is "<<NumYRowOdd<<endl;
            Info<<"number of rows z in direction is "<<NumZRow<<endl;
            Info << " spaceX = " << spaceX<<endl;
            Info << " spaceY = " << spaceY<<endl;
            Info << " spaceZ = " << spaceZ<<endl;
            Info << "NdiffX " << NdiffX;


            //Info << "runTime.rootPath(): " << runTime.rootPath() << endl;
            //Info << "runTime.caseName(): " << runTime.caseName() << endl;
            //Info << "runTime.path(): " << runTime.path() << endl;

            fileName rootName(runTime.path());
            fileName fName("packingCubicParticlePositions");   
            OFstream os(rootName/fName);
            Info<< "Writing to the output file"<<endl;


            //counter//
            label l1=0;
            label n=0;

            for(label k=0;k<NumZRow;k++)
            {
                if(k%2 == 0)
                {
                    for(label j=0;j<NumYRowEven;j++)
                    {
                        if(j%2 == 0)
                        {
                            for(label i=0;i<maxParticXEvenWhenZEven;i++)
                            {
                        
                        
                                l1=l1+1;
                                if(l1<=particleNumber)
                                {
                                    veclist[l1-1].z()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionZ;
                                    veclist[l1-1].x()=particleDiameter*i+firstPositionXEven;
                                    veclist[l1-1].y()=(j*particleDiameter*Foam::sin(constant::mathematical::pi/3))+firstPositionYEven;
                                    n = n+1;
                            
                                    os <<"("<<veclist[l1-1].x()<<"  "<<veclist[l1-1].y()<<"  "<<veclist[l1-1].z()<<")"<<endl;
                                }
                        
                            }
                        }
                        else
                        {
                            for(label i=0;i<maxParticXOddWhenZEven;i++)
                            {
                        
                        
                                l1=l1+1;
                                if(l1<=particleNumber)
                                {
                                    veclist[l1-1].z()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionZ;
                                    veclist[l1-1].x()=particleDiameter/2+particleDiameter*i+firstPositionXEven;
                                    veclist[l1-1].y()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*(j-1)+particleDiameter*Foam::sin(constant::mathematical::pi/3)+firstPositionYEven;
                                    n = n+1;
                            
                                    os <<"("<<veclist[l1-1].x()<<"  "<<veclist[l1-1].y()<<"  "<<veclist[l1-1].z()<<")"<<endl;
                                }
                        
                            }
                        }
                    }
                }
                else
                {
                    for(label j=0;j<NumYRowOdd;j++)
                    {
                        if(j%2 == 0)
                        {
                            for(label i=0;i<maxParticXEvenWhenZOdd;i++)
                            {
                        
                        
                                l1=l1+1;
                                if(l1<=particleNumber)
                                {
                                    veclist[l1-1].z()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionZ;
                                    veclist[l1-1].x()=particleDiameter*i+firstPositionXOdd;
                                    veclist[l1-1].y()=(j*particleDiameter*Foam::sin(constant::mathematical::pi/3))+firstPositionYOdd;
                                    n = n+1;
                            
                                    os <<"("<<veclist[l1-1].x()<<"  "<<veclist[l1-1].y()<<"  "<<veclist[l1-1].z()<<")"<<endl;
                                }
                        
                            }
                        }
                        else
                        {
                            for(label i=0;i<maxParticXOddWhenZOdd;i++)
                            {
                        
                        
                                l1=l1+1;
                                if(l1<=particleNumber)
                                {
                                    veclist[l1-1].z()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionZ;
                                    veclist[l1-1].x()=particleDiameter/2+particleDiameter*(i-1)+firstPositionXOdd;
                                    //veclist[l1-1].x()=particleDiameter*i+firstPositionXOdd;
                                    veclist[l1-1].y()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*(j-1)+particleDiameter*Foam::sin(constant::mathematical::pi/3)+firstPositionYOdd;
                                    n = n+1;
                            
                                    os <<"("<<veclist[l1-1].x()<<"  "<<veclist[l1-1].y()<<"  "<<veclist[l1-1].z()<<")"<<endl;
                                }
                        
                            }
                        }
                    }
                }
            }
            Info << "max Particle that can be in this box is " << l1 << endl;
            Info <<" number of particles in this box is " << n << endl;
        }
        else if(pileDirection == x)
        {
            //number of max particles in z direction in odd row and even row when x axis is even//
            label maxParticZEvenWhenXEven = boundaryLengthZ/particleDiameter;
            label maxParticZOddWhenXEven  = (boundaryLengthZ-particleDiameter/2)/particleDiameter;
            //number of max particles in x direction in odd row and even row when z axis is odd//
            label maxParticZEvenWhenXOdd  = (boundaryLengthZ-particleDiameter/2)/particleDiameter;
            label maxParticZOddWhenXOdd   = boundaryLengthZ/particleDiameter;
        
            //number of max particles in y direction when NumZrow is Even and Odd//
            //label NumYRowEven = (boundaryLengthY/Foam::sin(constant::mathematical::pi/3))/particleDiameter;
            scalar yL = particleDiameter * Foam::sin(constant::mathematical::pi/3);
            label NumYRowEven = (boundaryLengthY - particleDiameter)/yL + 1;
            label NumYRowOdd  = (boundaryLengthY-particleDiameter/2*Foam::tan(constant::mathematical::pi/6)-particleDiameter)/yL +1;
            //label NumYRowOdd  = ((boundaryLengthY-particleDiameter/2)/Foam::sin(constant::mathematical::pi/3))/particleDiameter;
        
            //number of max particles in z direction//
            label NumXRow = (boundaryLengthX/Foam::sin(constant::mathematical::pi/3))/particleDiameter;

            //Particle number difference between maxParticleXOdd and maxParticXeven
            label NdiffZ = maxParticZEvenWhenXEven-maxParticZOddWhenXEven;
            label NdiffY = NumYRowEven-NumYRowOdd;

            //space between the first (and the last) particle and boundary wall in z direction//
            scalar spaceZ = (boundaryLengthZ-(maxParticZEvenWhenXEven*particleDiameter+(-NdiffZ+1)*particleDiameter/2))/2;
            //space between the first (and the last) particle and boundary wall in y direction//
            //scalar spaceY = Foam::sin(constant::mathematical::pi/3)*(boundaryLengthY/Foam::sin(constant::mathematical::pi/3)-NumYRow*particleDiameter)/2+NdiffY*((particleDiameter/2)/Foam::sin(constant::mathematical::pi/6));
            scalar spaceY = (boundaryLengthY - ((NumYRowEven-1)*yL+particleDiameter+(-NdiffY+1)*particleDiameter/2*Foam::tan(constant::mathematical::pi/6)))/2;
            //space between the first (and the last) particle and boundary wall in z direction//
            scalar spaceX = Foam::sin(constant::mathematical::pi/3)*(boundaryLengthX/Foam::sin(constant::mathematical::pi/3)-NumXRow*particleDiameter)/2;
 
            //the first position of the particle//
            //the first Position of z and y when x axis is even
            scalar firstPositionZEven = boundaryPointMin.z()+spaceZ+particleDiameter/2;
            scalar firstPositionYEven = boundaryPointMin.y()+spaceY+particleDiameter/2;
            //the first Position of z and y when x axis is odd
            scalar firstPositionZOdd = firstPositionZEven + particleDiameter*Foam::cos(constant::mathematical::pi/3);
            scalar firstPositionYOdd = firstPositionYEven + particleDiameter/2*Foam::tan(constant::mathematical::pi/6);
            scalar firstPositionX = boundaryPointMin.x()+spaceX+particleDiameter/2;



            Info<<"number of max particles in z direction (Y even X even) is " << maxParticZEvenWhenXEven<<endl;
            Info<<"number of max particles in z direction (Y odd x even) is "<<maxParticZOddWhenXEven<<endl;
            Info<<"number of rows in y direction (Even) is "<<NumYRowEven<<endl;
            Info<<"number of rows in y direction (Odd) is "<<NumYRowOdd<<endl;
            Info<<"number of rows X in direction is "<<NumXRow<<endl;
            Info << " spaceZ = " << spaceZ<<endl;
            Info << " spaceY = " << spaceY<<endl;
            Info << " spaceX = " << spaceX<<endl;
            Info << "NdiffZ " << NdiffZ;


            //Info << "runTime.rootPath(): " << runTime.rootPath() << endl;
            //Info << "runTime.caseName(): " << runTime.caseName() << endl;
            //Info << "runTime.path(): " << runTime.path() << endl;

            fileName rootName(runTime.path());
            fileName fName("packingCubicParticlePositions");   
            OFstream os(rootName/fName);
            Info<< "Writing to the output file"<<endl;


            //counter//
            label l1=0;
            label n=0;

            for(label k=0;k<NumXRow;k++)
            {
                if(k%2 == 0)
                {
                    for(label j=0;j<NumYRowEven;j++)
                    {
                        if(j%2 == 0)
                        {
                            for(label i=0;i<maxParticZEvenWhenXEven;i++)
                            {
                        
                        
                                l1=l1+1;
                                if(l1<=particleNumber)
                                {
                                    veclist[l1-1].x()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionX;
                                    veclist[l1-1].z()=particleDiameter*i+firstPositionZEven;
                                    veclist[l1-1].y()=(j*particleDiameter*Foam::sin(constant::mathematical::pi/3))+firstPositionYEven;
                                    n = n+1;
                            
                                    os <<"("<<veclist[l1-1].x()<<"  "<<veclist[l1-1].y()<<"  "<<veclist[l1-1].z()<<")"<<endl;
                                }
                        
                            }
                        }
                        else
                        {
                            for(label i=0;i<maxParticZOddWhenXEven;i++)
                            {
                        
                        
                                l1=l1+1;
                                if(l1<=particleNumber)
                                {
                                    veclist[l1-1].x()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionX;
                                    veclist[l1-1].z()=particleDiameter/2+particleDiameter*i+firstPositionZEven;
                                    veclist[l1-1].y()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*(j-1)+particleDiameter*Foam::sin(constant::mathematical::pi/3)+firstPositionYEven;
                                    n = n+1;
                            
                                    os <<"("<<veclist[l1-1].x()<<"  "<<veclist[l1-1].y()<<"  "<<veclist[l1-1].z()<<")"<<endl;
                                }
                        
                            }
                        }
                    }
                }
                else
                {
                    for(label j=0;j<NumYRowOdd;j++)
                    {
                        if(j%2 == 0)
                        {
                            for(label i=0;i<maxParticZEvenWhenXOdd;i++)
                            {
                        
                        
                                l1=l1+1;
                                if(l1<=particleNumber)
                                {
                                    veclist[l1-1].x()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionX;
                                    veclist[l1-1].z()=particleDiameter*i+firstPositionZOdd;
                                    veclist[l1-1].y()=(j*particleDiameter*Foam::sin(constant::mathematical::pi/3))+firstPositionYOdd;
                                    n = n+1;
                            
                                    os <<"("<<veclist[l1-1].x()<<"  "<<veclist[l1-1].y()<<"  "<<veclist[l1-1].z()<<")"<<endl;
                                }
                        
                            }
                        }
                        else
                        {
                            for(label i=0;i<maxParticZOddWhenXOdd;i++)
                            {
                        
                        
                                l1=l1+1;
                                if(l1<=particleNumber)
                                {
                                    veclist[l1-1].x()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionX;
                                    veclist[l1-1].z()=particleDiameter/2+particleDiameter*(i-1)+firstPositionZOdd;
                                    //veclist[l1-1].x()=particleDiameter*i+firstPositionXOdd;
                                    veclist[l1-1].y()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*(j-1)+particleDiameter*Foam::sin(constant::mathematical::pi/3)+firstPositionYOdd;
                                    n = n+1;
                            
                                    os <<"("<<veclist[l1-1].x()<<"  "<<veclist[l1-1].y()<<"  "<<veclist[l1-1].z()<<")"<<endl;
                                }
                        
                            }
                        }
                    }
                }
            }
            Info << "max Particle that can be in this box is " << l1 << endl;
            Info <<" number of particles in this box is " << n << endl;
        }
        else if(pileDirection == y)
        {
            //number of max particles in x direction in odd row and even row when y axis is even//
            label maxParticXEvenWhenYEven = boundaryLengthX/particleDiameter;
            label maxParticXOddWhenYEven  = (boundaryLengthX-particleDiameter/2)/particleDiameter;
            //number of max particles in x direction in odd row and even row when y axis is odd//
            label maxParticXEvenWhenYOdd  = (boundaryLengthX-particleDiameter/2)/particleDiameter;
            label maxParticXOddWhenYOdd   = boundaryLengthX/particleDiameter;
        
            //number of max particles in z direction when NumYrow is Even and Odd//
            //label NumZRowEven = (boundaryLengthZ/Foam::sin(constant::mathematical::pi/3))/particleDiameter;
            scalar zL = particleDiameter * Foam::sin(constant::mathematical::pi/3);
            label NumZRowEven = (boundaryLengthZ - particleDiameter)/zL + 1;
            label NumZRowOdd  = (boundaryLengthZ-particleDiameter/2*Foam::tan(constant::mathematical::pi/6)-particleDiameter)/zL +1;
            //label NumZRowOdd  = ((boundaryLengthZ-particleDiameter/2)/Foam::sin(constant::mathematical::pi/3))/particleDiameter;
        
            //number of max particles in y direction//
            label NumYRow = (boundaryLengthY/Foam::sin(constant::mathematical::pi/3))/particleDiameter;

            //Particle number difference between maxParticleXOdd and maxParticXeven
            label NdiffX = maxParticXEvenWhenYEven-maxParticXOddWhenYEven;
            label NdiffZ = NumZRowEven-NumZRowOdd;

            //space between the first (and the last) particle and boundary wall in x direction//
            scalar spaceX = (boundaryLengthX-(maxParticXEvenWhenYEven*particleDiameter+(-NdiffX+1)*particleDiameter/2))/2;
            //space between the first (and the last) particle and boundary wall in z direction//
            //scalar spaceZ = Foam::sin(constant::mathematical::pi/3)*(boundaryLengthZ/Foam::sin(constant::mathematical::pi/3)-NumZRow*particleDiameter)/2+NdiffZ*((particleDiameter/2)/Foam::sin(constant::mathematical::pi/6));
            scalar spaceZ = (boundaryLengthZ - ((NumZRowEven-1)*zL+particleDiameter+(-NdiffZ+1)*particleDiameter/2*Foam::tan(constant::mathematical::pi/6)))/2;
            //space between the first (and the last) particle and boundary wall in y direction//
            scalar spaceY = Foam::sin(constant::mathematical::pi/3)*(boundaryLengthY/Foam::sin(constant::mathematical::pi/3)-NumYRow*particleDiameter)/2;
 
            //the first position of the particle//
            //the first Position of X and z when y axis is even
            scalar firstPositionXEven = boundaryPointMin.x()+spaceX+particleDiameter/2;
            scalar firstPositionZEven = boundaryPointMin.z()+spaceZ+particleDiameter/2;
            //the first Position of X and z when y axis is odd
            scalar firstPositionXOdd = firstPositionXEven + particleDiameter*Foam::cos(constant::mathematical::pi/3);
            scalar firstPositionZOdd = firstPositionZEven + particleDiameter/2*Foam::tan(constant::mathematical::pi/6);
            scalar firstPositionY = boundaryPointMin.y()+spaceY+particleDiameter/2;



            Info<<"number of max particles in x direction (Z even Y even) is " << maxParticXEvenWhenYEven<<endl;
            Info<<"number of max particles in x direction (Z odd Y even) is "<<maxParticXOddWhenYEven<<endl;
            Info<<"number of rows in z direction (Even) is "<<NumZRowEven<<endl;
            Info<<"number of rows in z direction (Odd) is "<<NumZRowOdd<<endl;
            Info<<"number of rows y in direction is "<<NumYRow<<endl;
            Info << " spaceX = " << spaceX<<endl;
            Info << " spaceZ = " << spaceZ<<endl;
            Info << " spaceY = " << spaceY<<endl;
            Info << "NdiffX " << NdiffX;


            //Info << "runTime.rootPath(): " << runTime.rootPath() << endl;
            //Info << "runTime.caseName(): " << runTime.caseName() << endl;
            //Info << "runTime.path(): " << runTime.path() << endl;

            fileName rootName(runTime.path());
            fileName fName("packingCubicParticlePositions");   
            OFstream os(rootName/fName);
            Info<< "Writing to the output file"<<endl;


            //counter//
            label l1=0;
            label n=0;

            for(label k=0;k<NumYRow;k++)
            {
                if(k%2 == 0)
                {
                    for(label j=0;j<NumZRowEven;j++)
                    {
                        if(j%2 == 0)
                        {
                            for(label i=0;i<maxParticXEvenWhenYEven;i++)
                            {
                        
                        
                                l1=l1+1;
                                if(l1<=particleNumber)
                                {
                                    veclist[l1-1].y()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionY;
                                    veclist[l1-1].x()=particleDiameter*i+firstPositionXEven;
                                    veclist[l1-1].z()=(j*particleDiameter*Foam::sin(constant::mathematical::pi/3))+firstPositionZEven;
                                    n = n+1;
                            
                                    os <<"("<<veclist[l1-1].x()<<"  "<<veclist[l1-1].y()<<"  "<<veclist[l1-1].z()<<")"<<endl;
                                }
                        
                            }
                        }
                        else
                        {
                            for(label i=0;i<maxParticXOddWhenYEven;i++)
                            {
                        
                        
                                l1=l1+1;
                                if(l1<=particleNumber)
                                {
                                    veclist[l1-1].y()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionY;
                                    veclist[l1-1].x()=particleDiameter/2+particleDiameter*i+firstPositionXEven;
                                    veclist[l1-1].z()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*(j-1)+particleDiameter*Foam::sin(constant::mathematical::pi/3)+firstPositionZEven;
                                    n = n+1;
                            
                                    os <<"("<<veclist[l1-1].x()<<"  "<<veclist[l1-1].y()<<"  "<<veclist[l1-1].z()<<")"<<endl;
                                }
                        
                            }
                        }
                    }
                }
                else
                {
                    for(label j=0;j<NumZRowOdd;j++)
                    {
                        if(j%2 == 0)
                        {
                            for(label i=0;i<maxParticXEvenWhenYOdd;i++)
                            {
                        
                        
                                l1=l1+1;
                                if(l1<=particleNumber)
                                {
                                    veclist[l1-1].y()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionY;
                                    veclist[l1-1].x()=particleDiameter*i+firstPositionXOdd;
                                    veclist[l1-1].z()=(j*particleDiameter*Foam::sin(constant::mathematical::pi/3))+firstPositionZOdd;
                                    n = n+1;
                            
                                    os <<"("<<veclist[l1-1].x()<<"  "<<veclist[l1-1].y()<<"  "<<veclist[l1-1].z()<<")"<<endl;
                                }
                        
                            }
                        }
                        else
                        {
                            for(label i=0;i<maxParticXOddWhenYOdd;i++)
                            {
                        
                        
                                l1=l1+1;
                                if(l1<=particleNumber)
                                {
                                    veclist[l1-1].y()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionY;
                                    veclist[l1-1].x()=particleDiameter/2+particleDiameter*(i-1)+firstPositionXOdd;
                                    //veclist[l1-1].x()=particleDiameter*i+firstPositionXOdd;
                                    veclist[l1-1].z()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*(j-1)+particleDiameter*Foam::sin(constant::mathematical::pi/3)+firstPositionZOdd;
                                    n = n+1;
                            
                                    os <<"("<<veclist[l1-1].x()<<"  "<<veclist[l1-1].y()<<"  "<<veclist[l1-1].z()<<")"<<endl;
                                }
                        
                            }
                        }
                    }
                }
            }
            Info << "max Particle that can be in this box is " << l1 << endl;
            Info <<" number of particles in this box is " << n << endl;
        }



        Info<< "End\n" << endl;
    }

    if(positioningType == "packingCylinder")
    {

        dictionary packingCylinderSubDict
        (
            setParticlePositionsDict.subDict("packingCylinderSubDict")
        );

        scalar particleDiameter
        (
            readScalar(packingCylinderSubDict.lookup("particleDiameter"))
        );

        label particleNumber
        (
            readScalar(packingCylinderSubDict.lookup("particleNumber"))
        );

        scalar cylinderRadius
        (
            readScalar(packingCylinderSubDict.lookup("cylinderRadius"))
        );

        scalar cylinderLength
        (
            readScalar(packingCylinderSubDict.lookup("cylinderLength"))
        );

        vector axisDirection
        (
            packingCylinderSubDict.lookup("axisDirection")
        );
        
        vector origin
        (
            packingCylinderSubDict.lookup("origin")
        );

        

        vector x  (1,0,0);
        vector y  (0,1,0);
        vector z  (0,0,1);

        List<vector> veclist;
        veclist.setSize(particleNumber);


        if(axisDirection == z)
        {
            //Boundary is taken little larger than the cylinder
            scalar boundaryLengthX = cylinderRadius*2+particleDiameter;
            scalar boundaryLengthY = cylinderRadius*2+particleDiameter;
            scalar boundaryLengthZ = cylinderLength+particleDiameter;
            vector boundaryPointMin;
            boundaryPointMin.x() = origin.x()-boundaryLengthX/2;
            boundaryPointMin.y() = origin.y()-boundaryLengthY/2;
            boundaryPointMin.z() = origin.z()-boundaryLengthZ/2;


            //number of max particles in x direction in odd row and even row when z axis is even//
            label maxParticXEvenWhenZEven = boundaryLengthX/particleDiameter;
            label maxParticXOddWhenZEven  = (boundaryLengthX-particleDiameter/2)/particleDiameter;
            //number of max particles in x direction in odd row and even row when z axis is odd//
            label maxParticXEvenWhenZOdd  = (boundaryLengthX-particleDiameter/2)/particleDiameter;
            label maxParticXOddWhenZOdd   = boundaryLengthX/particleDiameter;
        
            //number of max particles in y direction when NumZrow is Even and Odd//
            //label NumYRowEven = (boundaryLengthY/Foam::sin(constant::mathematical::pi/3))/particleDiameter;
            scalar yL = particleDiameter * Foam::sin(constant::mathematical::pi/3);
            label NumYRowEven = (boundaryLengthY - particleDiameter)/yL + 1;
            label NumYRowOdd  = (boundaryLengthY-particleDiameter/2*Foam::tan(constant::mathematical::pi/6)-particleDiameter)/yL +1;
            //label NumYRowOdd  = ((boundaryLengthY-particleDiameter/2)/Foam::sin(constant::mathematical::pi/3))/particleDiameter;
        
            //number of max particles in z direction//
            label NumZRow = (boundaryLengthZ/Foam::sin(constant::mathematical::pi/3))/particleDiameter;

            //Particle number difference between maxParticleXOdd and maxParticXeven
            label NdiffX = maxParticXEvenWhenZEven-maxParticXOddWhenZEven;
            label NdiffY = NumYRowEven-NumYRowOdd;

            //space between the first (and the last) particle and boundary wall in x direction//
            scalar spaceX = (boundaryLengthX-(maxParticXEvenWhenZEven*particleDiameter+(-NdiffX+1)*particleDiameter/2))/2;
            //space between the first (and the last) particle and boundary wall in y direction//
            //scalar spaceY = Foam::sin(constant::mathematical::pi/3)*(boundaryLengthY/Foam::sin(constant::mathematical::pi/3)-NumYRow*particleDiameter)/2+NdiffY*((particleDiameter/2)/Foam::sin(constant::mathematical::pi/6));
            scalar spaceY = (boundaryLengthY - ((NumYRowEven-1)*yL+particleDiameter+(-NdiffY+1)*particleDiameter/2*Foam::tan(constant::mathematical::pi/6)))/2;
            //space between the first (and the last) particle and boundary wall in z direction//
            scalar spaceZ = Foam::sin(constant::mathematical::pi/3)*(boundaryLengthZ/Foam::sin(constant::mathematical::pi/3)-NumZRow*particleDiameter)/2;
 
            //the first position of the particle//
            //the first Position of X and y when z axis is even
            scalar firstPositionXEven = boundaryPointMin.x()+spaceX+particleDiameter/2;
            scalar firstPositionYEven = boundaryPointMin.y()+spaceY+particleDiameter/2;
            //the first Position of X and y when z axis is odd
            scalar firstPositionXOdd = firstPositionXEven + particleDiameter*Foam::cos(constant::mathematical::pi/3);
            scalar firstPositionYOdd = firstPositionYEven + particleDiameter/2*Foam::tan(constant::mathematical::pi/6);
            scalar firstPositionZ = boundaryPointMin.z()+spaceZ+particleDiameter/2;



            Info<<"number of max particles in x direction (Y even Z even) is " << maxParticXEvenWhenZEven<<endl;
            Info<<"number of max particles in x direction (Y odd Z even) is "<<maxParticXOddWhenZEven<<endl;
            Info<<"number of rows in y direction (Even) is "<<NumYRowEven<<endl;
            Info<<"number of rows in y direction (Odd) is "<<NumYRowOdd<<endl;
            Info<<"number of rows z in direction is "<<NumZRow<<endl;
            Info << " spaceX = " << spaceX<<endl;
            Info << " spaceY = " << spaceY<<endl;
            Info << " spaceZ = " << spaceZ<<endl;
            Info << "NdiffX " << NdiffX;

            //Dummy list to store particle position of rectangular coordinates
            label num = maxParticXEvenWhenZEven*NumYRowEven*NumZRow;
            List<vector> dummylist;
            dummylist.setSize(num);


            //Info << "runTime.rootPath(): " << runTime.rootPath() << endl;
            //Info << "runTime.caseName(): " << runTime.caseName() << endl;
            //Info << "runTime.path(): " << runTime.path() << endl;

            fileName rootName(runTime.path());
            fileName fName("packingCylinderParticlePositions");   
            OFstream os(rootName/fName);
            Info<< "Writing to the output file"<<endl;


            //counter//
            label l1=0;
            label n=0;

            for(label k=0;k<NumZRow;k++)
            {
                if(k%2 == 0)
                {
                    for(label j=0;j<NumYRowEven;j++)
                    {
                        if(j%2 == 0)
                        {
                            for(label i=0;i<maxParticXEvenWhenZEven;i++)
                            {
                        
                        
                                l1=l1+1;
                                
                                dummylist[l1-1].z()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionZ;
                                dummylist[l1-1].x()=particleDiameter*i+firstPositionXEven;
                                dummylist[l1-1].y()=(j*particleDiameter*Foam::sin(constant::mathematical::pi/3))+firstPositionYEven;
                                n = n+1;
                            
                                //os <<"("<<dummylist[l1-1].x()<<"  "<<dummylist[l1-1].y()<<"  "<<dummylist[l1-1].z()<<")"<<endl;
                                
                        
                            }
                        }
                        else
                        {
                            for(label i=0;i<maxParticXOddWhenZEven;i++)
                            {
                        
                        
                                l1=l1+1;
                                
                                dummylist[l1-1].z()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionZ;
                                dummylist[l1-1].x()=particleDiameter/2+particleDiameter*i+firstPositionXEven;
                                dummylist[l1-1].y()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*(j-1)+particleDiameter*Foam::sin(constant::mathematical::pi/3)+firstPositionYEven;
                                n = n+1;
                            
                                //os <<"("<<dummylist[l1-1].x()<<"  "<<dummylist[l1-1].y()<<"  "<<dummylist[l1-1].z()<<")"<<endl;
                                
                        
                            }
                        }
                    }
                }
                else
                {
                    for(label j=0;j<NumYRowOdd;j++)
                    {
                        if(j%2 == 0)
                        {
                            for(label i=0;i<maxParticXEvenWhenZOdd;i++)
                            {
                        
                        
                                l1=l1+1;
                                
                                dummylist[l1-1].z()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionZ;
                                dummylist[l1-1].x()=particleDiameter*i+firstPositionXOdd;
                                dummylist[l1-1].y()=(j*particleDiameter*Foam::sin(constant::mathematical::pi/3))+firstPositionYOdd;
                                n = n+1;
                            
                                //os <<"("<<dummylist[l1-1].x()<<"  "<<dummylist[l1-1].y()<<"  "<<dummylist[l1-1].z()<<")"<<endl;
                                
                        
                            }
                        }
                        else
                        {
                            for(label i=0;i<maxParticXOddWhenZOdd;i++)
                            {
                        
                        
                                l1=l1+1;
                                
                                dummylist[l1-1].z()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionZ;
                                dummylist[l1-1].x()=particleDiameter/2+particleDiameter*(i-1)+firstPositionXOdd;
                                //dummylist[l1-1].x()=particleDiameter*i+firstPositionXOdd;
                                dummylist[l1-1].y()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*(j-1)+particleDiameter*Foam::sin(constant::mathematical::pi/3)+firstPositionYOdd;
                                n = n+1;
                            
                                //os <<"("<<dummylist[l1-1].x()<<"  "<<dummylist[l1-1].y()<<"  "<<dummylist[l1-1].z()<<")"<<endl;
                                
                        
                            }
                        }
                    }
                }
            }
            Info << "max Particle that can be in this box is " << l1 << endl;
            Info <<" number of particles in this box is " << n << endl;

            
            //Check whether particle is inside the cylinder
            
            scalar LengthfromAxis;
            scalar AxisDirecLength;
            label k=0;
            label i=0;

            for(i=0; i<n; i++)
            {
                LengthfromAxis = Foam::sqrt(pow(dummylist[i].x()-origin.x(),2)+pow(dummylist[i].y()-origin.y(),2));
                AxisDirecLength = Foam::mag(dummylist[i].z()-origin.z());
                if(  ( (LengthfromAxis+particleDiameter/2)<=cylinderRadius )&&( (AxisDirecLength+particleDiameter/2)<=(cylinderLength/2) )  )
                {   
                    
                    if( k < particleNumber)
                    {   
                        k = k+1;
                        veclist[k-1] = dummylist[i];
                        os <<"("<<veclist[k-1].x()<<"  "<<veclist[k-1].y()<<"  "<<veclist[k-1].z()<<")"<<endl;
                    }
                }
            }
            Info << "number of particles in this cylinder is " << k << endl;
                


        }
        else if(axisDirection == x)
        {

            //Boundary is taken little larger than the cylinder
            scalar boundaryLengthZ = cylinderRadius*2+particleDiameter;
            scalar boundaryLengthY = cylinderRadius*2+particleDiameter;
            scalar boundaryLengthX = cylinderLength+particleDiameter;
            vector boundaryPointMin;
            boundaryPointMin.x() = origin.x()-boundaryLengthX/2;
            boundaryPointMin.y() = origin.y()-boundaryLengthY/2;
            boundaryPointMin.z() = origin.z()-boundaryLengthZ/2;

            //number of max particles in z direction in odd row and even row when x axis is even//
            label maxParticZEvenWhenXEven = boundaryLengthZ/particleDiameter;
            label maxParticZOddWhenXEven  = (boundaryLengthZ-particleDiameter/2)/particleDiameter;
            //number of max particles in x direction in odd row and even row when z axis is odd//
            label maxParticZEvenWhenXOdd  = (boundaryLengthZ-particleDiameter/2)/particleDiameter;
            label maxParticZOddWhenXOdd   = boundaryLengthZ/particleDiameter;
        
            //number of max particles in y direction when NumZrow is Even and Odd//
            //label NumYRowEven = (boundaryLengthY/Foam::sin(constant::mathematical::pi/3))/particleDiameter;
            scalar yL = particleDiameter * Foam::sin(constant::mathematical::pi/3);
            label NumYRowEven = (boundaryLengthY - particleDiameter)/yL + 1;
            label NumYRowOdd  = (boundaryLengthY-particleDiameter/2*Foam::tan(constant::mathematical::pi/6)-particleDiameter)/yL +1;
            //label NumYRowOdd  = ((boundaryLengthY-particleDiameter/2)/Foam::sin(constant::mathematical::pi/3))/particleDiameter;
        
            //number of max particles in z direction//
            label NumXRow = (boundaryLengthX/Foam::sin(constant::mathematical::pi/3))/particleDiameter;

            //Particle number difference between maxParticleXOdd and maxParticXeven
            label NdiffZ = maxParticZEvenWhenXEven-maxParticZOddWhenXEven;
            label NdiffY = NumYRowEven-NumYRowOdd;

            //space between the first (and the last) particle and boundary wall in z direction//
            scalar spaceZ = (boundaryLengthZ-(maxParticZEvenWhenXEven*particleDiameter+(-NdiffZ+1)*particleDiameter/2))/2;
            //space between the first (and the last) particle and boundary wall in y direction//
            //scalar spaceY = Foam::sin(constant::mathematical::pi/3)*(boundaryLengthY/Foam::sin(constant::mathematical::pi/3)-NumYRow*particleDiameter)/2+NdiffY*((particleDiameter/2)/Foam::sin(constant::mathematical::pi/6));
            scalar spaceY = (boundaryLengthY - ((NumYRowEven-1)*yL+particleDiameter+(-NdiffY+1)*particleDiameter/2*Foam::tan(constant::mathematical::pi/6)))/2;
            //space between the first (and the last) particle and boundary wall in z direction//
            scalar spaceX = Foam::sin(constant::mathematical::pi/3)*(boundaryLengthX/Foam::sin(constant::mathematical::pi/3)-NumXRow*particleDiameter)/2;
 
            //the first position of the particle//
            //the first Position of z and y when x axis is even
            scalar firstPositionZEven = boundaryPointMin.z()+spaceZ+particleDiameter/2;
            scalar firstPositionYEven = boundaryPointMin.y()+spaceY+particleDiameter/2;
            //the first Position of z and y when x axis is odd
            scalar firstPositionZOdd = firstPositionZEven + particleDiameter*Foam::cos(constant::mathematical::pi/3);
            scalar firstPositionYOdd = firstPositionYEven + particleDiameter/2*Foam::tan(constant::mathematical::pi/6);
            scalar firstPositionX = boundaryPointMin.x()+spaceX+particleDiameter/2;



            Info<<"number of max particles in z direction (Y even X even) is " << maxParticZEvenWhenXEven<<endl;
            Info<<"number of max particles in z direction (Y odd x even) is "<<maxParticZOddWhenXEven<<endl;
            Info<<"number of rows in y direction (Even) is "<<NumYRowEven<<endl;
            Info<<"number of rows in y direction (Odd) is "<<NumYRowOdd<<endl;
            Info<<"number of rows X in direction is "<<NumXRow<<endl;
            Info << " spaceZ = " << spaceZ<<endl;
            Info << " spaceY = " << spaceY<<endl;
            Info << " spaceX = " << spaceX<<endl;
            Info << "NdiffZ " << NdiffZ;


            //Info << "runTime.rootPath(): " << runTime.rootPath() << endl;
            //Info << "runTime.caseName(): " << runTime.caseName() << endl;
            //Info << "runTime.path(): " << runTime.path() << endl;

            //Dummy list to store particle position of rectangular coordinates
            label num = maxParticZEvenWhenXEven*NumYRowEven*NumXRow;
            List<vector> dummylist;
            dummylist.setSize(num);

            fileName rootName(runTime.path());
            fileName fName("packingCylinderParticlePositions");   
            OFstream os(rootName/fName);
            Info<< "Writing to the output file"<<endl;


            //counter//
            label l1=0;
            label n=0;

            for(label k=0;k<NumXRow;k++)
            {
                if(k%2 == 0)
                {
                    for(label j=0;j<NumYRowEven;j++)
                    {
                        if(j%2 == 0)
                        {
                            for(label i=0;i<maxParticZEvenWhenXEven;i++)
                            {
                        
                        
                                l1=l1+1;
                                
                                    dummylist[l1-1].x()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionX;
                                    dummylist[l1-1].z()=particleDiameter*i+firstPositionZEven;
                                    dummylist[l1-1].y()=(j*particleDiameter*Foam::sin(constant::mathematical::pi/3))+firstPositionYEven;
                                    n = n+1;
                            
                                    //os <<"("<<dummylist[l1-1].x()<<"  "<<dummylist[l1-1].y()<<"  "<<dummylist[l1-1].z()<<")"<<endl;
                                
                        
                            }
                        }
                        else
                        {
                            for(label i=0;i<maxParticZOddWhenXEven;i++)
                            {
                        
                        
                                l1=l1+1;
                                
                                    dummylist[l1-1].x()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionX;
                                    dummylist[l1-1].z()=particleDiameter/2+particleDiameter*i+firstPositionZEven;
                                    dummylist[l1-1].y()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*(j-1)+particleDiameter*Foam::sin(constant::mathematical::pi/3)+firstPositionYEven;
                                    n = n+1;
                            
                                    //os <<"("<<dummylist[l1-1].x()<<"  "<<dummylist[l1-1].y()<<"  "<<dummylist[l1-1].z()<<")"<<endl;
                                
                        
                            }
                        }
                    }
                }
                else
                {
                    for(label j=0;j<NumYRowOdd;j++)
                    {
                        if(j%2 == 0)
                        {
                            for(label i=0;i<maxParticZEvenWhenXOdd;i++)
                            {
                        
                        
                                l1=l1+1;
                                
                                    dummylist[l1-1].x()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionX;
                                    dummylist[l1-1].z()=particleDiameter*i+firstPositionZOdd;
                                    dummylist[l1-1].y()=(j*particleDiameter*Foam::sin(constant::mathematical::pi/3))+firstPositionYOdd;
                                    n = n+1;
                            
                                    //os <<"("<<dummylist[l1-1].x()<<"  "<<dummylist[l1-1].y()<<"  "<<dummylist[l1-1].z()<<")"<<endl;
                                
                        
                            }
                        }
                        else
                        {
                            for(label i=0;i<maxParticZOddWhenXOdd;i++)
                            {
                        
                        
                                l1=l1+1;
                                
                                    dummylist[l1-1].x()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionX;
                                    dummylist[l1-1].z()=particleDiameter/2+particleDiameter*(i-1)+firstPositionZOdd;
                                    //dummylist[l1-1].x()=particleDiameter*i+firstPositionXOdd;
                                    dummylist[l1-1].y()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*(j-1)+particleDiameter*Foam::sin(constant::mathematical::pi/3)+firstPositionYOdd;
                                    n = n+1;
                            
                                    //os <<"("<<dummylist[l1-1].x()<<"  "<<dummylist[l1-1].y()<<"  "<<dummylist[l1-1].z()<<")"<<endl;
                                
                        
                            }
                        }
                    }
                }
            }
            Info << "max Particle that can be in this box is " << l1 << endl;
            Info <<" number of particles in this box is " << n << endl;

            //Check whether particle is inside the cylinder
            
            scalar LengthfromAxis;
            scalar AxisDirecLength;
            label k=0;
            label i=0;

            for(i=0; i<n; i++)
            {
                LengthfromAxis = Foam::sqrt(pow(dummylist[i].z()-origin.z(),2)+pow(dummylist[i].y()-origin.y(),2));
                AxisDirecLength = Foam::mag(dummylist[i].x()-origin.x());
                if(  ( (LengthfromAxis+particleDiameter/2)<=cylinderRadius )&&( (AxisDirecLength+particleDiameter/2)<=(cylinderLength/2) )  )
                {   
                    
                    if( k < particleNumber)
                    {   
                        k = k+1;
                        veclist[k-1] = dummylist[i];
                        os <<"("<<veclist[k-1].x()<<"  "<<veclist[k-1].y()<<"  "<<veclist[k-1].z()<<")"<<endl;
                    }
                }
            }
            Info << "number of particles in this cylinder is " << k << endl;

        }
        else if(axisDirection == y)
        {

            //Boundary is taken little larger than the cylinder
            scalar boundaryLengthZ = cylinderRadius*2+particleDiameter;
            scalar boundaryLengthX = cylinderRadius*2+particleDiameter;
            scalar boundaryLengthY = cylinderLength+particleDiameter;
            vector boundaryPointMin;
            boundaryPointMin.x() = origin.x()-boundaryLengthX/2;
            boundaryPointMin.y() = origin.y()-boundaryLengthY/2;
            boundaryPointMin.z() = origin.z()-boundaryLengthZ/2;

            //number of max particles in x direction in odd row and even row when y axis is even//
            label maxParticXEvenWhenYEven = boundaryLengthX/particleDiameter;
            label maxParticXOddWhenYEven  = (boundaryLengthX-particleDiameter/2)/particleDiameter;
            //number of max particles in x direction in odd row and even row when y axis is odd//
            label maxParticXEvenWhenYOdd  = (boundaryLengthX-particleDiameter/2)/particleDiameter;
            label maxParticXOddWhenYOdd   = boundaryLengthX/particleDiameter;
        
            //number of max particles in z direction when NumYrow is Even and Odd//
            //label NumZRowEven = (boundaryLengthZ/Foam::sin(constant::mathematical::pi/3))/particleDiameter;
            scalar zL = particleDiameter * Foam::sin(constant::mathematical::pi/3);
            label NumZRowEven = (boundaryLengthZ - particleDiameter)/zL + 1;
            label NumZRowOdd  = (boundaryLengthZ-particleDiameter/2*Foam::tan(constant::mathematical::pi/6)-particleDiameter)/zL +1;
            //label NumZRowOdd  = ((boundaryLengthZ-particleDiameter/2)/Foam::sin(constant::mathematical::pi/3))/particleDiameter;
        
            //number of max particles in y direction//
            label NumYRow = (boundaryLengthY/Foam::sin(constant::mathematical::pi/3))/particleDiameter;

            //Particle number difference between maxParticleXOdd and maxParticXeven
            label NdiffX = maxParticXEvenWhenYEven-maxParticXOddWhenYEven;
            label NdiffZ = NumZRowEven-NumZRowOdd;

            //space between the first (and the last) particle and boundary wall in x direction//
            scalar spaceX = (boundaryLengthX-(maxParticXEvenWhenYEven*particleDiameter+(-NdiffX+1)*particleDiameter/2))/2;
            //space between the first (and the last) particle and boundary wall in z direction//
            //scalar spaceZ = Foam::sin(constant::mathematical::pi/3)*(boundaryLengthZ/Foam::sin(constant::mathematical::pi/3)-NumZRow*particleDiameter)/2+NdiffZ*((particleDiameter/2)/Foam::sin(constant::mathematical::pi/6));
            scalar spaceZ = (boundaryLengthZ - ((NumZRowEven-1)*zL+particleDiameter+(-NdiffZ+1)*particleDiameter/2*Foam::tan(constant::mathematical::pi/6)))/2;
            //space between the first (and the last) particle and boundary wall in y direction//
            scalar spaceY = Foam::sin(constant::mathematical::pi/3)*(boundaryLengthY/Foam::sin(constant::mathematical::pi/3)-NumYRow*particleDiameter)/2;
 
            //the first position of the particle//
            //the first Position of X and z when y axis is even
            scalar firstPositionXEven = boundaryPointMin.x()+spaceX+particleDiameter/2;
            scalar firstPositionZEven = boundaryPointMin.z()+spaceZ+particleDiameter/2;
            //the first Position of X and z when y axis is odd
            scalar firstPositionXOdd = firstPositionXEven + particleDiameter*Foam::cos(constant::mathematical::pi/3);
            scalar firstPositionZOdd = firstPositionZEven + particleDiameter/2*Foam::tan(constant::mathematical::pi/6);
            scalar firstPositionY = boundaryPointMin.y()+spaceY+particleDiameter/2;



            Info<<"number of max particles in x direction (Z even Y even) is " << maxParticXEvenWhenYEven<<endl;
            Info<<"number of max particles in x direction (Z odd Y even) is "<<maxParticXOddWhenYEven<<endl;
            Info<<"number of rows in z direction (Even) is "<<NumZRowEven<<endl;
            Info<<"number of rows in z direction (Odd) is "<<NumZRowOdd<<endl;
            Info<<"number of rows y in direction is "<<NumYRow<<endl;
            Info << " spaceX = " << spaceX<<endl;
            Info << " spaceZ = " << spaceZ<<endl;
            Info << " spaceY = " << spaceY<<endl;
            Info << "NdiffX " << NdiffX;

            //Dummy list to store particle position of rectangular coordinates
            label num = maxParticXEvenWhenYEven*NumZRowEven*NumYRow;
            List<vector> dummylist;
            dummylist.setSize(num);


            //Info << "runTime.rootPath(): " << runTime.rootPath() << endl;
            //Info << "runTime.caseName(): " << runTime.caseName() << endl;
            //Info << "runTime.path(): " << runTime.path() << endl;

            fileName rootName(runTime.path());
            fileName fName("packingCylinderParticlePositions");   
            OFstream os(rootName/fName);
            Info<< "Writing to the output file"<<endl;


            //counter//
            label l1=0;
            label n=0;

            for(label k=0;k<NumYRow;k++)
            {
                if(k%2 == 0)
                {
                    for(label j=0;j<NumZRowEven;j++)
                    {
                        if(j%2 == 0)
                        {
                            for(label i=0;i<maxParticXEvenWhenYEven;i++)
                            {
                        
                        
                                l1=l1+1;
                                
                                    dummylist[l1-1].y()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionY;
                                    dummylist[l1-1].x()=particleDiameter*i+firstPositionXEven;
                                    dummylist[l1-1].z()=(j*particleDiameter*Foam::sin(constant::mathematical::pi/3))+firstPositionZEven;
                                    n = n+1;
                            
                                    //os <<"("<<dummylist[l1-1].x()<<"  "<<dummylist[l1-1].y()<<"  "<<dummylist[l1-1].z()<<")"<<endl;
                                
                        
                            }
                        }
                        else
                        {
                            for(label i=0;i<maxParticXOddWhenYEven;i++)
                            {
                        
                        
                                l1=l1+1;
                                
                                    dummylist[l1-1].y()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionY;
                                    dummylist[l1-1].x()=particleDiameter/2+particleDiameter*i+firstPositionXEven;
                                    dummylist[l1-1].z()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*(j-1)+particleDiameter*Foam::sin(constant::mathematical::pi/3)+firstPositionZEven;
                                    n = n+1;
                            
                                    //os <<"("<<dummylist[l1-1].x()<<"  "<<dummylist[l1-1].y()<<"  "<<dummylist[l1-1].z()<<")"<<endl;
                                
                        
                            }
                        }
                    }
                }
                else
                {
                    for(label j=0;j<NumZRowOdd;j++)
                    {
                        if(j%2 == 0)
                        {
                            for(label i=0;i<maxParticXEvenWhenYOdd;i++)
                            {
                        
                        
                                l1=l1+1;
                                
                                    dummylist[l1-1].y()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionY;
                                    dummylist[l1-1].x()=particleDiameter*i+firstPositionXOdd;
                                    dummylist[l1-1].z()=(j*particleDiameter*Foam::sin(constant::mathematical::pi/3))+firstPositionZOdd;
                                    n = n+1;
                            
                                    //os <<"("<<dummylist[l1-1].x()<<"  "<<dummylist[l1-1].y()<<"  "<<dummylist[l1-1].z()<<")"<<endl;
                                
                        
                            }
                        }
                        else
                        {
                            for(label i=0;i<maxParticXOddWhenYOdd;i++)
                            {
                        
                        
                                l1=l1+1;
                                
                                    dummylist[l1-1].y()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*k+firstPositionY;
                                    dummylist[l1-1].x()=particleDiameter/2+particleDiameter*(i-1)+firstPositionXOdd;
                                    //dummylist[l1-1].x()=particleDiameter*i+firstPositionXOdd;
                                    dummylist[l1-1].z()=(particleDiameter*Foam::sin(constant::mathematical::pi/3))*(j-1)+particleDiameter*Foam::sin(constant::mathematical::pi/3)+firstPositionZOdd;
                                    n = n+1;
                            
                                    //os <<"("<<dummylist[l1-1].x()<<"  "<<dummylist[l1-1].y()<<"  "<<dummylist[l1-1].z()<<")"<<endl;
                                
                        
                            }
                        }
                    }
                }
            }
            Info << "max Particle that can be in this box is " << l1 << endl;
            Info <<" number of particles in this box is " << n << endl;

            //Check whether particle is inside the cylinder
            
            scalar LengthfromAxis;
            scalar AxisDirecLength;
            label k=0;
            label i=0;

            for(i=0; i<n; i++)
            {
                LengthfromAxis = Foam::sqrt(pow(dummylist[i].z()-origin.z(),2)+pow(dummylist[i].x()-origin.x(),2));
                AxisDirecLength = Foam::mag(dummylist[i].y()-origin.y());
                if(  ( (LengthfromAxis+particleDiameter/2)<=cylinderRadius )&&( (AxisDirecLength+particleDiameter/2)<=(cylinderLength/2) )  )
                {   
                    
                    if( k < particleNumber)
                    {   
                        k = k+1;
                        veclist[k-1] = dummylist[i];
                        os <<"("<<veclist[k-1].x()<<"  "<<veclist[k-1].y()<<"  "<<veclist[k-1].z()<<")"<<endl;
                    }
                }
            }
            Info << "number of particles in this cylinder is " << k << endl;

        }



        Info<< "End\n" << endl;
    }



    
    return 0;
}

double GetRandom(double min, double max)
{
    return min + rand()*(max-min+1.0)/(1.0+RAND_MAX);
}


// ************************************************************************* //
