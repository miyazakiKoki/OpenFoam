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

\*---------------------------------------------------------------------------*/

#include "basicWetCollidingCloud.H"

#include "makeParcelCloudFunctionObjects.H"

// Kinematic
#include "makeParcelForces.H"
#include "makeParcelDispersionModels.H"
#include "makeParcelInjectionModels.H"
#include "makeParcelCollisionModels.H"
#include "makeParcelPatchInteractionModels.H"
#include "makeParcelSurfaceFilmModels.H"

// Include user defined submodels
#include "makeUserDefinedParcelCloudFunctionObjects.H"
#include "makeUserDefinedParcelForces.H"
#include "makeUserDefinedParcelCollisionModels.H"

// Wet
#include "makeParcelWetModels.H"
#include "makeWetParcelInjectionModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeParcelCloudFunctionObjects(basicWetCollidingCloud);

    // Kinematic sub-models
    makeParcelForces(basicWetCollidingCloud);
    makeParcelDispersionModels(basicWetCollidingCloud);
    makeParcelInjectionModels(basicWetCollidingCloud);
    typedef basicWetCollidingCloud::collidingCloudType wetCollidingType;
    makeParcelCollisionModels(wetCollidingType);
    makeParcelPatchInteractionModels(basicWetCollidingCloud);
    makeParcelSurfaceFilmModels(basicWetCollidingCloud);

    // User defined sub-models
    makeUserDefinedParcelCloudFunctionObjects(basicWetCollidingCloud);
    makeUserDefinedParcelForces(basicWetCollidingCloud);
    makeUserDefinedParcelCollisionModels(wetCollidingType);

    // Wet sub-models
    makeParcelWetModels(basicWetCollidingCloud);
    makeWetParcelInjectionModels(basicWetCollidingCloud);

}


// ************************************************************************* //
