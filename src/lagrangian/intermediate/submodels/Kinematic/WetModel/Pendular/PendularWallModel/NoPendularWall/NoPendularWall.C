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

#include "NoPendularWall.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoPendularWall<CloudType>::NoPendularWall
(
    const dictionary& dict,
    CloudType& cloud,
    const scalar& surfaceTension,
    const scalar& contactAngle,
    const scalar& liqFrac,
    const scalar& viscosity,
    const scalar& minSep
)
:
    PendularWallModel<CloudType>
    (
        dict,
        cloud,
        surfaceTension,
        contactAngle,
        liqFrac,
        viscosity,
        minSep
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoPendularWall<CloudType>::~NoPendularWall()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::NoPendularWall<CloudType>::pREff
(
    const typename CloudType::parcelType& p
) const
{
    return p.d()/2;
}

template<class CloudType>
void Foam::NoPendularWall<CloudType>::evaluatePendularWall
(
    typename CloudType::parcelType& p,
    const List<point>& flatSitePoints,
    const List<WallSiteData<vector> >& flatSiteData,
    const List<point>& sharpSitePoints,
    const List<WallSiteData<vector> >& sharpSiteData
) const
{}



// ************************************************************************* //
