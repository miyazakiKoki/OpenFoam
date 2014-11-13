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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::WetParcel<ParcelType>::WetParcel
(
    const polyMesh& owner,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    ParcelType(owner, position, cellI, tetFaceI, tetPtI),
    Vliq_(0.0),
    partVliq_(),
    liquidPositionVectors_(),
    liquidPositions_(),
    contactList_(),
    previousContactList_()
{}


template<class ParcelType>
inline Foam::WetParcel<ParcelType>::WetParcel
(
    const polyMesh& owner,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const scalar dTarget0,
    const vector& U0,
    const vector& f0,
    const vector& angularMomentum0,
    const vector& torque0,
    const typename ParcelType::constantProperties& constProps,
    const scalar Vliq,
    const List<scalar> partVliq,
    const List<vector> liquidPositionVectors,
    const List<vector> liquidPositions,
    const List<label> contactList,
    const List<label> previousContactList
)
:
    ParcelType
    (
        owner,
        position,
        cellI,
        tetFaceI,
        tetPtI,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        f0,
        angularMomentum0,
        torque0,
        constProps
    ),
    Vliq_(Vliq),
    partVliq_(partVliq),
    liquidPositionVectors_(liquidPositionVectors),
    liquidPositions_(liquidPositions),
    contactList_(contactList),
    previousContactList_(previousContactList)
{}

// * * * * * * * WetParcel Member Functions  * * * * * * * //

template<class ParcelType>
inline Foam::scalar Foam::WetParcel<ParcelType>::Vliq() const
{
    return Vliq_;
}


template<class ParcelType>
inline Foam::scalar& Foam::WetParcel<ParcelType>::Vliq()
{
    return Vliq_;
}

template<class ParcelType>
inline Foam::List<Foam::scalar> Foam::WetParcel<ParcelType>::partVliq() const
{
    return partVliq_;
}

template<class ParcelType>
inline Foam::List<Foam::scalar>& Foam::WetParcel<ParcelType>::partVliq()
{
    return partVliq_;
}

template<class ParcelType>
inline Foam::List<Foam::vector> Foam::WetParcel<ParcelType>::liquidPositions() const
{
    return liquidPositions_;
}

template<class ParcelType>
inline Foam::List<Foam::vector>& Foam::WetParcel<ParcelType>::liquidPositions()
{
    return liquidPositions_;
}

template<class ParcelType>
inline Foam::List<Foam::vector> Foam::WetParcel<ParcelType>::liquidPositionVectors() const
{
    return liquidPositionVectors_;
}

template<class ParcelType>
inline Foam::List<Foam::vector>& Foam::WetParcel<ParcelType>::liquidPositionVectors()
{
    return liquidPositionVectors_;
}

template<class ParcelType>
inline Foam::List<Foam::label> Foam::WetParcel<ParcelType>::contactList() const
{
    return contactList_;
}

template<class ParcelType>
inline Foam::List<Foam::label>& Foam::WetParcel<ParcelType>::contactList()
{
    return contactList_;
}

template<class ParcelType>
inline Foam::List<Foam::label> Foam::WetParcel<ParcelType>::previousContactList() const
{
    return previousContactList_;
}

template<class ParcelType>
inline Foam::List<Foam::label>& Foam::WetParcel<ParcelType>::previousContactList()
{
    return previousContactList_;
}

// ************************************************************************* //
