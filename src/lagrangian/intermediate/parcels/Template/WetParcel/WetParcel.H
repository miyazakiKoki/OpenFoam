/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Class
    Foam::WetParcel

Description
    Wrapper around kinematic parcel types to add collision modelling

SourceFiles
    WetParcelI.H
    WetParcel.C
    WetParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef WetParcel_H
#define WetParcel_H

#include "particle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class WetParcel;

// Forward declaration of friend functions

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const WetParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                        Class CollidingParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class WetParcel
:
    public ParcelType
{
protected:

    // Protected members

        //Total liquid volume deposited on particle surface
        scalar Vliq_;
        //Liquid volume deposited on part of the particle surface
        List<scalar> partVliq_;

        // Liquid position vector(Unit) on particle surface
        List<vector> liquidPositions_;
        List<vector> liquidPositionVectors_;

public:

    // Static data members

        //- Runtime type information
        TypeName("WetParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " Vliq"
          + "partVliq"
          + " liquidPositions"
          + " liquidPositionVectors"
        );

public:
    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline WetParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI
        );

        //- Construct from components
        inline WetParcel
        (
            const polyMesh& mesh,
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
            const List<vector> liquidPositions
        );

        //- Construct from Istream
        WetParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        WetParcel(const WetParcel& p);

        //- Construct as a copy
        WetParcel(const WetParcel& p, const polyMesh& mesh);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new WetParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new WetParcel(*this, mesh));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<WetParcel<ParcelType> > operator()(Istream& is) const
            {
                return autoPtr<WetParcel<ParcelType> >
                (
                    new WetParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to the liquid volume
            inline scalar Vliq() const;

            inline List<scalar> partVliq() const;

            inline List<vector> liquidPositionVectors() const;

            inline List<vector> liquidPositions() const;


        // Edit

            //- Return access to the liquid volume
            inline scalar& Vliq();

            inline List<scalar>& partVliq();

            inline List<vector>& liquidPositionVectors();

            inline List<vector>& liquidPositions();

        // Tracking

            //- Move the parcel
            template<class TrackData>
            bool move(TrackData& td, const scalar trackTime);


       // I-O

            //- Read
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType>
            static void writeFields(const CloudType& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const WetParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "WetParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "WetParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
