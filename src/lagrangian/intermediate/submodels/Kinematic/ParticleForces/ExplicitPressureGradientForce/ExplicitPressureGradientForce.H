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
    Foam::ExplicitPressureGradientForce

Description
    Calculates particle gravity force

SourceFiles
    ExplicitPressureGradientForceI.H
    ExplicitPressureGradientForce.C

\*---------------------------------------------------------------------------*/

#ifndef ExplicitPressureGradientForce_H
#define ExplicitPressureGradientForce_H

#include "ParticleForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                       Class ExplicitPressureGradientForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ExplicitPressureGradientForce
:
    public ParticleForce<CloudType>
{
    //- Private Member Functions

        //- Check the dimension of pressure term
        bool isDividedByRho_;


public:

    //- Runtime type information
    TypeName("explicitPressureGradient");


    // Constructors

        //- Construct from mesh
        ExplicitPressureGradientForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        ExplicitPressureGradientForce(const ExplicitPressureGradientForce& gf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType> > clone() const
        {
            return autoPtr<ParticleForce<CloudType> >
            (
                new ExplicitPressureGradientForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ExplicitPressureGradientForce();


    // Member Functions

        // Evaluation

            //- Calculate the non-coupled force
            virtual forceSuSp calcNonCoupled
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

#include "ExplicitPressureGradientForceI.H"

#ifdef NoRepository
    #include "ExplicitPressureGradientForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
