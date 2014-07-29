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

#ifndef makeParcelWetModels_H
#define makeParcelWetModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NoWet.H"
#include "Pendular.H"

#include "NoPendular.H"
#include "Lambert.H"
#include "LambertDispersion.H"
#include "LambertDispersionNoWetP.H"

#include "NoPendularWall.H"
#include "LambertWall.H"
#include "LambertWallDispersion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeParcelWetModels(CloudType)                                        \
                                                                              \
    makeWetModel(CloudType);                                                  \
    makeWetModelType(NoWet, CloudType);                                       \
    makeWetModelType(Pendular, CloudType);                                    \
                                                                              \
    makePendularModel(CloudType);                                             \
    makePendularModelType(NoPendular, CloudType);                             \
    makePendularModelType(Lambert, CloudType);                                \
    makePendularModelType(LambertDispersion, CloudType);                      \
    makePendularModelType(LambertDispersionNoWetP, CloudType);                \
                                                                              \
    makePendularWallModel(CloudType);                                         \
    makePendularWallModelType(NoPendularWall, CloudType);                     \
    makePendularWallModelType(LambertWall, CloudType);                        \
    makePendularWallModelType(LambertWallDispersion, CloudType);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
