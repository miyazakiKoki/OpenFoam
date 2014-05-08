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

#include "TemplateCollidingParcel.H"
#include "IOstreams.H"
#include "IOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::TemplateCollidingParcel<ParcelType>::propertyList_ =
    Foam::TemplateCollidingParcel<ParcelType>::propertyList();


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::TemplateCollidingParcel<ParcelType>::TemplateCollidingParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields)
{
    if (readFields)
    {
    }

    // Check state of Istream
    is.check
    (
        "TemplateCollidingParcel<ParcelType>::TemplateCollisionParcel"
        "(const polyMesh&, Istream&, bool)"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::TemplateCollidingParcel<ParcelType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParcelType::readFields(c);
}


template<class ParcelType>
template<class CloudType>
void Foam::TemplateCollidingParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const TemplateCollidingParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p);
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const TemplateCollidingParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
