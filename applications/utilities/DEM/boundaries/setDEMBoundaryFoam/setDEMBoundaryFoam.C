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
    setDEMBoundaryFoam

Description
    Set boundary conditions used for DEM simulation.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    runTime++;

    Info << "\nChanging boundary conditions\n" << endl;

    word boundaryType
    (
        setDEMBoundaryDict.lookup("boundaryType")
    );

    if(boundaryType == "rotaryDrum")
    {
        dictionary rotaryDrumSubDict
        (
            setDEMBoundaryDict.subDict("rotaryDrumSubDict")
        );

        vector centreOfRotation
        (
            rotaryDrumSubDict.lookup("centreOfRotation")
        );

        vector axisOfRotation
        (
            rotaryDrumSubDict.lookup("axisOfRotation")
        );

        axisOfRotation /= (mag(axisOfRotation)+VSMALL);

        scalar rpm
        (
            readScalar(rotaryDrumSubDict.lookup("rpm"))
        );

        forAll(U.boundaryField(), patchI)
        {
            forAll(U.boundaryField()[patchI], faceI)
            {
                vector facePos = mesh.boundaryMesh()[patchI].faceCentres()[faceI];
                vector posFromCentre = facePos - centreOfRotation;
                vector posFromAxis =
                    posFromCentre - (posFromCentre&axisOfRotation)*axisOfRotation;
                vector faceVel =
                    2.*constant::mathematical::pi*rpm/60. * (posFromAxis^axisOfRotation);
                U.boundaryField()[patchI][faceI] = faceVel;
            }
        }

        U.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
