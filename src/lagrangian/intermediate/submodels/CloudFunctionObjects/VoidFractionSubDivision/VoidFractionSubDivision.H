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
    Foam::VoidFractionSubDivision

Description
    Creates particle void fraction field on carrier phase

SourceFiles
    VoidFractionSubDivision.C

\*---------------------------------------------------------------------------*/

#ifndef VoidFractionSubDivision_H
#define VoidFractionSubDivision_H

#include "CloudFunctionObject.H"
#include "volFields.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class VoidFractionSubDivision Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class VoidFractionSubDivision
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;


        //- Void fraction field
        autoPtr<volScalarField> thetaPtr_;


protected:

    // Protected Member Functions

        //- Write post-processing info
        virtual void write();

        label refinementLevel_;

        List<vector> elementPosition_;

        scalar annulusThicknessFactor_;

        void recursivelyCalcVoidFraction
        (
            const parcelType& p,
            const label cellI,
            const scalar &dt,
            labelHashSet &lhash
        );

public:

    //- Runtime type information
    TypeName("voidFractionSubDivision");


    // Constructors

        //- Construct from dictionary
        VoidFractionSubDivision(const dictionary& dict, CloudType& owner);

        //- Construct copy
        VoidFractionSubDivision(const VoidFractionSubDivision<CloudType>& vf);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType> > clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType> >
            (
                new VoidFractionSubDivision<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~VoidFractionSubDivision();


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
                bool& keepParticle
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "VoidFractionSubDivision.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
