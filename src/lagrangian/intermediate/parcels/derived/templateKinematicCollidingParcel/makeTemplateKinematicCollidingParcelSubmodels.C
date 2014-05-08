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

#include "templateKinematicCollidingCloud.H"

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeParcelCloudFunctionObjects(templateKinematicCollidingCloud);

    // Kinematic sub-models
    makeParcelForces(templateKinematicCollidingCloud);
    makeParcelDispersionModels(templateKinematicCollidingCloud);
    makeParcelInjectionModels(templateKinematicCollidingCloud);
    makeParcelCollisionModels(templateKinematicCollidingCloud);
    makeParcelPatchInteractionModels(templateKinematicCollidingCloud);
    makeParcelSurfaceFilmModels(templateKinematicCollidingCloud);

    // User defined sub-models
    makeUserDefinedParcelCloudFunctionObjects(templateKinematicCollidingCloud);
    makeUserDefinedParcelForces(templateKinematicCollidingCloud);

}


// ************************************************************************* //
