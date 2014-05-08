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

#include "WetModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::WetModel<CloudType>::WetModel(CloudType& owner)
:
    SubModelBase<CloudType>(owner)
{}


template<class CloudType>
Foam::WetModel<CloudType>::WetModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    SubModelBase<CloudType>(owner, dict, typeName, type)
{}


template<class CloudType>
Foam::WetModel<CloudType>::WetModel(WetModel<CloudType>& cm)
:
    SubModelBase<CloudType>(cm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::WetModel<CloudType>::~WetModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::WetModel<CloudType>::controlsWallInteraction() const
{
    notImplemented
    (
        "bool Foam::WetModel<CloudType>::controlsWallInteraction()"
    );
    return false;
}

template<class CloudType>
void Foam::WetModel<CloudType>::bond()
{
    notImplemented("void Foam::WetModel<CloudType>::bond()");
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "WetModelNew.C"

// ************************************************************************* //
