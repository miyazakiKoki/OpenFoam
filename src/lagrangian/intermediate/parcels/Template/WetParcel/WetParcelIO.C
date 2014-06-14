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

#include "WetParcel.H"
#include "IOstreams.H"
#include "IOField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::WetParcel<ParcelType>::propertyList_ =
    Foam::WetParcel<ParcelType>::propertyList();


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::WetParcel<ParcelType>::WetParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    Vliq_(0.0),
    partVliq_(0),
    liquidPositions_(),
    liquidPositionVectors_()
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            Vliq_ = readScalar(is);
            partVliq_ = readList<scalar>(is);
            liquidPositions_ = readList<vector>(is);
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&Vliq_),
                sizeof(Vliq_)
            );
            is.read
            (
                reinterpret_cast<char*>(&partVliq_),
                sizeof(partVliq_)
            );
            is.read
            (
                reinterpret_cast<char*>(&liquidPositions_),
                sizeof(liquidPositions_)
            );
            is.read
            (
                reinterpret_cast<char*>(&liquidPositionVectors_),
                sizeof(liquidPositionVectors_)
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "WetParcel<ParcelType>::WetParcel"
        "(const polyMesh&, Istream&, bool)"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::WetParcel<ParcelType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParcelType::readFields(c);

    IOField<scalar> Vliq(c.fieldIOobject("Vliq", IOobject::MUST_READ));
    c.checkFieldIOobject(c, Vliq);

    IOField<scalarField>  partVliq(c.fieldIOobject("partVliq", IOobject::MUST_READ));
    c.checkFieldIOobject(c, partVliq);

    IOField<vectorField>  liquidPositions(c.fieldIOobject("liquidPositions", IOobject::MUST_READ));
    c.checkFieldIOobject(c, liquidPositions);

    IOField<vectorField>  liquidPositionVectors(c.fieldIOobject("liquidPositionVectors", IOobject::MUST_READ));
    c.checkFieldIOobject(c, liquidPositionVectors);


    label i = 0;

    forAllIter(typename CloudType, c, iter)
    {
        WetParcel<ParcelType>& p = iter();

        p.Vliq_ = Vliq[i];
        p.partVliq_ = partVliq[i];
        p.liquidPositions_ = liquidPositions[i];
        p.liquidPositionVectors_ = liquidPositionVectors[i];

        i++;
    }

}


template<class ParcelType>
template<class CloudType>
void Foam::WetParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    label np =  c.size();

    IOField<scalar> Vliq(c.fieldIOobject("Vliq", IOobject::NO_READ), np);
    IOField<scalarField>  partVliq(c.fieldIOobject("partVliq", IOobject::NO_READ), np);
    IOField<vectorField>  liquidPositions(c.fieldIOobject("liquidPositions", IOobject::NO_READ), np);
    IOField<vectorField>  liquidPositionVectors(c.fieldIOobject("liquidPositionVectors", IOobject::NO_READ), np);

    label i = 0;

    forAllConstIter(typename CloudType, c, iter)
    {
        const WetParcel<ParcelType>& p = iter();

        Vliq[i] = p.Vliq();
        partVliq[i] = p.partVliq();
        liquidPositions[i] = p.liquidPositions();
        liquidPositionVectors[i] = p.liquidPositionVectors();



        i++;
    }

    Vliq.write();
    partVliq.write();
    liquidPositions.write();
    liquidPositionVectors.write();

}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const WetParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.Vliq();
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.partVliq();
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.liquidPositions();
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.liquidPositionVectors();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.Vliq_),
            sizeof(p.Vliq())

        );
        os.write
        (
            reinterpret_cast<const char*>(&p.partVliq_),
            sizeof(p.partVliq())
        );
        os.write
        (
            reinterpret_cast<const char*>(&p.liquidPositions_),
            sizeof(p.liquidPositions())
        );
        os.write
        (
            reinterpret_cast<const char*>(&p.liquidPositionVectors_),
            sizeof(p.liquidPositionVectors())
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const WetParcel<ParcelType>&)"
    );

    return os;
}
/*






template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const WetParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.Vliq();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.Vliq_),
            sizeof(p.Vliq())
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const WetParcel<ParcelType>&)"
    );

    return os;
}
*/

// ************************************************************************* //
