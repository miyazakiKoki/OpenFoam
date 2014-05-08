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
    Foam::CellAverageParticleVelocity

Description
    Creates particle void fraction field on carrier phase

SourceFiles
    CellAverageParticleVelocity.C

\*---------------------------------------------------------------------------*/

#ifndef CellAverageParticleVelocity_H
#define CellAverageParticleVelocity_H

#include "CloudFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class CellAverageParticleVelocity Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CellAverageParticleVelocity
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;


        //- Average particle velocity field
        autoPtr<volVectorField> UpPtr_;

        //- Count total particle volume in cell
        autoPtr<volScalarField> pVolPtr_;

protected:

    // Protected Member Functions

        //- Write post-processing info
        virtual void write();


public:

    //- Runtime type information
    TypeName("cellAverageParticleVelocity");


    // Constructors

        //- Construct from dictionary
        CellAverageParticleVelocity(const dictionary& dict, CloudType& owner);

        //- Construct copy
        CellAverageParticleVelocity(const CellAverageParticleVelocity<CloudType>& vf);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType> > clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType> >
            (
                new CellAverageParticleVelocity<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~CellAverageParticleVelocity();


    // Member Functions

        // Evaluation

            //- Pre-evolve hook
            virtual void preEvolve();

            //- Post-evolve hook
            virtual void postEvolve();

            //- Post-move hook
            virtual void postMove
            (
                const parcelType& p,
                const label cellI,
                const scalar dt,
                const point& position0,
                bool& keepParticle
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CellAverageParticleVelocity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
