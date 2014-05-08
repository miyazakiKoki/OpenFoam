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
    Foam::NoPendularWall

Description
    Forces between particles and walls, interacting with a spring,
    slider, damper model

\*---------------------------------------------------------------------------*/

#ifndef NoPendularWall_H
#define NoPendularWall_H

#include "PendularWallModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class NoPendularWall Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoPendularWall
:
    public PendularWallModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoPendularWall
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
    virtual ~NoPendularWall();


    // Member Functions

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
#   include "NoPendularWall.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
