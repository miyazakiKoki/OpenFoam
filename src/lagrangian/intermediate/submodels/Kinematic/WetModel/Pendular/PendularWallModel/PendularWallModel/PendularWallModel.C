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

#include "PendularWallModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PendularWallModel<CloudType>::PendularWallModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type,
    const scalar& surfaceTension,
    const scalar& contactAngle,
    const scalar& liqFrac,
    const scalar& viscosity,
    const scalar& minSep
)
:
    dict_(dict),
    owner_(owner),
    coeffDict_(dict.subDict(type + "Coeffs")),
    surfaceTension_(surfaceTension),
    contactAngle_(contactAngle),
    liqFrac_(liqFrac),
    viscosity_(viscosity),
    minSep_(minSep)
{}

template<class CloudType>
Foam::PendularWallModel<CloudType>::PendularWallModel
(
    const dictionary& dict,
    CloudType& owner,
    const scalar& surfaceTension,
    const scalar& contactAngle,
    const scalar& liqFrac,
    const scalar& viscosity,
    const scalar& minSep
)
:
    dict_(dict),
    owner_(owner),
    surfaceTension_(surfaceTension),
    contactAngle_(contactAngle),
    liqFrac_(liqFrac),
    viscosity_(viscosity),
    minSep_(minSep)

{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PendularWallModel<CloudType>::~PendularWallModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const CloudType&
Foam::PendularWallModel<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
CloudType&
Foam::PendularWallModel<CloudType>::owner()
{
    return owner_;
}


template<class CloudType>
const Foam::dictionary& Foam::PendularWallModel<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
const Foam::dictionary&
Foam::PendularWallModel<CloudType>::coeffDict() const
{
    return coeffDict_;
}


template<class CloudType>
const Foam::scalar&
Foam::PendularWallModel<CloudType>::surfaceTension() const
{
    return surfaceTension_;
}


template<class CloudType>
const Foam::scalar&
Foam::PendularWallModel<CloudType>::contactAngle() const
{
    return contactAngle_;
}


template<class CloudType>
const Foam::scalar&
Foam::PendularWallModel<CloudType>::liqFrac() const
{
    return liqFrac_;
}

template<class CloudType>
const Foam::scalar&
Foam::PendularWallModel<CloudType>::viscosity() const
{
    return viscosity_;
}

template<class CloudType>
const Foam::scalar&
Foam::PendularWallModel<CloudType>::minSep() const
{
    return minSep_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PendularWallModelNew.C"

// ************************************************************************* //
