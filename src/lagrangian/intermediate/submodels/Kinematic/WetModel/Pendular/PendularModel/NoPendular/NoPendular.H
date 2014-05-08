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
    Foam::NoPendular

Description
    Pair forces between particles colliding with a spring, slider, damper model

\*---------------------------------------------------------------------------*/

#ifndef NoPendular_H
#define NoPendular_H

#include "PendularModel.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class NoPendular Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoPendular
:
    public PendularModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoPendular
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
    virtual ~NoPendular();

        //- Calculate the pair interaction between parcels
        virtual void evaluatePendular
        (
            typename CloudType::parcelType& pA,
            typename CloudType::parcelType& pB
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NoPendular.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
