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

Class
    Foam::LambertWall

Description
    Forces between particles and walls, interacting with a spring,
    slider, damper model

\*---------------------------------------------------------------------------*/

#ifndef LambertWall_H
#define LambertWall_H

#include "PendularWallModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class LambertWall Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class LambertWall
:
    public PendularWallModel<CloudType>
{
    // Private data

        //- Volume factor for determining the equivalent size of a
        //  parcel where nParticles is not 1.  The equivalent size of
        //  the parcel is
        //      parcelEquivVolume = volumeFactor*nParticles*p.volume()
        //  so
        //      parcelEquivD = cbrt(volumeFactor*nParticles)*p.d()
        //  + When volumeFactor = 1, the particles are compressed
        //    together so that the equivalent volume of the parcel is
        //    the sum of the constituent particles
        //  + When volumeFactor = 3*sqrt(2)/pi, the particles are
        //    close packed, but uncompressed.
        //  + When volumeFactor > 3*sqrt(2)/pi, the particles loosely
        //    grouped.
        // 3*sqrt(2)/pi = 1.350474 is the volume factor for close
        // packing, i.e pi/(3*sqrt(2)) is the maximum close packing
        // factor
        scalar volumeFactor_;

        //- Switch to control use of equivalent size particles.  Used
        //  because the calculation can be very expensive.
        bool useEquivalentSize_;


    // Private Member Functions

        //- Calculate the wall interaction for a parcel at a given site
        void evaluatePendularWall
        (
            typename CloudType::parcelType& p,
            const point& site,
            const WallSiteData<vector>& data,
            scalar pREff
        ) const;


public:

    //- Runtime type information
    TypeName("LambertWall");


    // Constructors

        //- Construct from dictionary
        LambertWall
        (
            const dictionary& dict,
            CloudType& cloud,
            const scalar& surfaceTension,
            const scalar& contactAngle,
            const scalar& liqFrac,
            const scalar& viscosity,
            const scalar& minSep
        );


    //- Destructor
    virtual ~LambertWall();


    // Member Functions

        //- Return the volumeFactor
        inline scalar volumeFactor() const
        {
            return volumeFactor_;
        }

        //- Return the effective radius for a particle for the model
        virtual scalar pREff(const typename CloudType::parcelType& p) const;

        //- Calculate the wall interaction for a parcel
        virtual void evaluatePendularWall
        (
            typename CloudType::parcelType& p,
            const List<point>& flatSitePoints,
            const List<WallSiteData<vector> >& flatSiteData,
            const List<point>& sharpSitePoints,
            const List<WallSiteData<vector> >& sharpSiteData
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "LambertWall.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
