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
        setWetPointDict.lookup("positioningType")
    );

    if(positioningType == "random")
    {
        dictionary randomSubDict
        (
            setWetPointDict.subDict("randomSubDict")
        );

        label wetPointNumber
        (
            readScalar(randomSubDict.lookup("wetPointNumber"))
        );

        scalar radius
        (
            readScalar(randomSubDict.lookup("particleRadius"))
        );

        scalar minEdgeLength
        (
            readScalar(randomSubDict.lookup("minEdgeLength"))
        );

        scalar thetaDensity
        (
            readScalar(randomSubDict.lookup("thetaDensity"))
        );

        scalar phiDensity
        (
            readScalar(randomSubDict.lookup("phiDensity"))
        );

        label trialNum
        (
            readScalar(randomSubDict.lookup("trialNumber"))
        );

        scalar thetaMin = 0;
        scalar thetaMax = 180 * thetaDensity;
        scalar phiMin = 0;
        scalar phiMax = 360 * phiDensity;
        
        //dummy list to store wetPoints in spherical coordinates x is r, y is theta, z is phi
        List<vector> dummySpheric(trialNum);

        //dummy list to store particles in rectangula coordinates
        List<vector> dummyRectangula(trialNum);

        //List for particles settled in spherical coordinates
        List<vector> spheric(wetPointNumber);

        //List for particles transformed into rectangular coordinates
        List<vector> rectangula(wetPointNumber);


        // the first particle is stored in cylindic and rectangula list

        spheric[0].x() = radius;
        spheric[0].y() = GetRandom(thetaMin,thetaMax)/thetaDensity;
        spheric[0].z() = GetRandom(phiMin,phiMax)/phiDensity;

        Info << "r is" << spheric[0].x() << endl;
        Info << "theta is" << spheric[0].y() << endl;
        Info << "phi is " << spheric[0].z() << endl;

        // spheric[0] is transormed into rectangular coordinates
        //inside sin is written in radian

        rectangula[0].x() = radius * Foam::sin(constant::mathematical::pi*spheric[0].y()/180) * Foam::cos(constant::mathematical::pi*spheric[0].z()/180);
        rectangula[0].y() = radius * Foam::sin(constant::mathematical::pi*spheric[0].y()/180) * Foam::sin(constant::mathematical::pi*spheric[0].z()/180);
        rectangula[0].z() = radius * Foam::cos(constant::mathematical::pi*spheric[0].y()/180);

        Info << "x is" << rectangula[0].x() << endl;
        Info << "y is" << rectangula[0].y() << endl;
        Info << "z is " << rectangula[0].z() << endl;

        label l = 0;

        forAll(dummySpheric,i)
        {
            if( l< wetPointNumber-1 )
            {
                dummySpheric[i].x() = radius;
                dummySpheric[i].y() = GetRandom(thetaMin,thetaMax)/thetaDensity;
                dummySpheric[i].z() = GetRandom(phiMin,phiMax)/phiDensity;
                Info << " dummySpheric[ "<<i<<"].x()" << dummySpheric[i].x() << endl;
                Info << " dummySpheric[ "<<i<<"].y()" << dummySpheric[i].y() << endl;
                Info << " dummySpheric[ "<<i<<"].z()" << dummySpheric[i].z() << endl;

                dummyRectangula[i].x() = dummySpheric[i].x() * Foam::sin(constant::mathematical::pi*dummySpheric[i].y()/180) * Foam::cos(constant::mathematical::pi*dummySpheric[i].z()/180);
                dummyRectangula[i].y() = dummySpheric[i].x() * Foam::sin(constant::mathematical::pi*dummySpheric[i].y()/180) * Foam::sin(constant::mathematical::pi*dummySpheric[i].z()/180);
                dummyRectangula[i].z() = dummySpheric[i].x() * Foam::cos(constant::mathematical::pi*dummySpheric[i].y()/180);


                label counter = 0;
                label t;
                for(t=0; t<=l; t++)
                {
                    scalar S;
                    S = mag(dummyRectangula[i]-rectangula[t]);
                    if(S>minEdgeLength)
                    {
                        counter += 1;   
                    }else{
                        counter += 0;
                    }
                }
                Info<<"value of counter is "<< counter << endl;
             
                if(counter == l+1)
                {
                    l = l+1;
                    rectangula[l] = dummyRectangula[i];
                    Info << "x is" << rectangula[l].x() << endl;
                    Info << "y is" << rectangula[l].y() << endl;
                    Info << "z is " << rectangula[l].z() << endl;
                }
            }
        }

        fileName rootName(runTime.path());
        fileName fName("randomWetPoint");   
        OFstream os(rootName/fName);
        Info<< "Writing to the output file"<<endl;

        label k;

        for(k=0; k<=l; k++)
        {
            os <<"("<<rectangula[k].x()<<"  "<<rectangula[k].y()<<"  "<<rectangula[k].z()<<")"<<endl;
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
