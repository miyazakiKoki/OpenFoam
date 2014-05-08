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

#include "VoidFractionDivided.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class CloudType>
label Foam::VoidFractionDivided<CloudType>::refinementLevel() const
{
    return refinementLevel_;
}

template<class CloudType>
label Foam::VoidFractionDivided<CloudType>::satParNo() const
{
    return satParNo_;
}

template<class CloudType>
List<vector> Foam::VoidFractionDivided<CloudType>::satPositions() const
{
    return satPositions_;
}

template<class CloudType>
scalar Foam::VoidFractionDivided<CloudType>::maxDistance() const
{
    return maxDistance_;
}

// ************************************************************************* //
