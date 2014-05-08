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
    Foam::ErgunPlusWenYuDragForce

Description
    Drag model based on assumption of solid spheres

\*---------------------------------------------------------------------------*/

#ifndef ErgunPlusWenYuDragForce_H
#define ErgunPlusWenYuDragForce_H

#include "ParticleForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class ErgunPlusWenYuDragForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ErgunPlusWenYuDragForce
:
    public ParticleForce<CloudType>
{
    // Private Member Functions

        //- Drag coefficient multiplied by Reynolds number
        scalar Cd(const scalar Re) const;

        scalar beta
        (
            const typename CloudType::parcelType& p,
            const scalar Re,
            const scalar muc,
            const scalar vfp
        ) const;

        scalar ErgunBeta
        (
            const typename CloudType::parcelType& p,
            const scalar Re,
            const scalar muc,
            const scalar vfp
        ) const;

        scalar WenYuBeta
        (
            const typename CloudType::parcelType& p,
            const scalar Re,
            const scalar muc,
            const scalar vfp
        ) const;

public:

    //- Runtime type information
    TypeName("ErgunPlusWenYuDrag");


    // Constructors

        //- Construct from mesh
        ErgunPlusWenYuDragForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        ErgunPlusWenYuDragForce(const ErgunPlusWenYuDragForce<CloudType>& df);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType> > clone() const
        {
            return autoPtr<ParticleForce<CloudType> >
            (
                new ErgunPlusWenYuDragForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ErgunPlusWenYuDragForce();


    // Member Functions

        // Evaluation

            //- Calculate the coupled force
            virtual forceSuSp calcCoupled
            (
                const typename CloudType::parcelType& p,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ErgunPlusWenYuDragForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
