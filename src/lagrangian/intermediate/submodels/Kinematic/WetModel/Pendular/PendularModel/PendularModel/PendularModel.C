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

#include "PendularModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PendularModel<CloudType>::PendularModel
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
Foam::PendularModel<CloudType>::PendularModel
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
Foam::PendularModel<CloudType>::~PendularModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const CloudType&
Foam::PendularModel<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
const Foam::dictionary& Foam::PendularModel<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
const Foam::dictionary&
Foam::PendularModel<CloudType>::coeffDict() const
{
    return coeffDict_;
}


template<class CloudType>
const Foam::scalar&
Foam::PendularModel<CloudType>::surfaceTension() const
{
    return surfaceTension_;
}


template<class CloudType>
const Foam::scalar&
Foam::PendularModel<CloudType>::contactAngle() const
{
    return contactAngle_;
}


template<class CloudType>
const Foam::scalar&
Foam::PendularModel<CloudType>::liqFrac() const
{
    return liqFrac_;
}

template<class CloudType>
const Foam::scalar&
Foam::PendularModel<CloudType>::viscosity() const
{
    return viscosity_;
}

template<class CloudType>
const Foam::scalar&
Foam::PendularModel<CloudType>::minSep() const
{
    return minSep_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PendularModelNew.C"

// ************************************************************************* //
