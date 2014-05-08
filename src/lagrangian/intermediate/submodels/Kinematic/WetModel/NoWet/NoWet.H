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
    Foam::NoWet

Description
    Place holder for 'none' option

SourceFiles
    NoWet.C

\*---------------------------------------------------------------------------*/

#ifndef NoWet_H
#define NoWet_H

#include "WetModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class NoWet Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoWet
:
    public WetModel<CloudType>
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        NoWet(const dictionary& dict, CloudType& owner);

        //- Construct copy
        NoWet(NoWet<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<WetModel<CloudType> > clone()
        {
            return autoPtr<WetModel<CloudType> >
            (
                new NoWet<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoWet();


    // Member Functions

        //- Flag to indicate whether model activates wet model
        virtual bool active() const;

        //- Indicates whether model determines wall bonding or not,
        //  used to determine what value to use for wallImpactDistance
        virtual bool controlsWallInteraction() const;

        // Bond function
        virtual void bond();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NoWet.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
