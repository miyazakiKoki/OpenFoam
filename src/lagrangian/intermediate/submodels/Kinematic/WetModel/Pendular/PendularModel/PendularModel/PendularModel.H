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
    Foam::PendularModel

Description
    Templated pendular interaction class

SourceFiles
    PendularModel.C
    PendularModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef PendularModel_H
#define PendularModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class PendularModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PendularModel
{
    // Private data

        //- The CollisionModel dictionary
        const dictionary& dict_;

        //- Reference to the owner cloud class
        CloudType& owner_;

        //- The coefficients dictionary
        const dictionary coeffDict_;

        //- Surface tension coefficient
        scalar surfaceTension_;

        //- Contact angle
        scalar contactAngle_;

        //- Fraction of the liquid to form one bridge
        scalar liqFrac_;

        //- Viscosity
        scalar viscosity_;

        //- Minimum separation distance divided by effective radius
        scalar minSep_;

public:

    //- Runtime type information
    TypeName("pendularModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        PendularModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner,
            const scalar& surfaceTension,
            const scalar& contactAngle,
            const scalar& liqFrac,
            const scalar& viscosity,
            const scalar& minSep
        ),
        (
            dict,
            owner,
            surfaceTension,
            contactAngle,
            liqFrac,
            viscosity,
            minSep
        )
    );


    // Constructors

        //- Construct from components
        PendularModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type,
            const scalar& surfaceTension,
            const scalar& contactAngle,
            const scalar& liqFrac,
            const scalar& viscosity,
            const scalar& minSep
        );

        //- Construct from components
        PendularModel
        (
            const dictionary& dict,
            CloudType& owner,
            const scalar& surfaceTension,
            const scalar& contactAngle,
            const scalar& liqFrac,
            const scalar& viscosity,
            const scalar& minSep
        );

    //- Destructor
    virtual ~PendularModel();


    //- Selector
    static autoPtr<PendularModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner,
        const scalar& surfaceTension,
        const scalar& contactAngle,
        const scalar& liqFrac,
        const scalar& viscosity,
        const scalar& minSep
    );


    // Access

        //- Return the owner cloud object
        const CloudType& owner() const;

        //- Return the dictionary
        const dictionary& dict() const;

        //- Return the coefficients dictionary
        const dictionary& coeffDict() const;

        //- Calculate the pendular interaction between parcels
        virtual void evaluatePendular
        (
            typename CloudType::parcelType& pA,
            typename CloudType::parcelType& pB
        ) const = 0;

        //- Return the surfaceTension
        const scalar& surfaceTension() const;

        //- Return the contactAngle
        const scalar& contactAngle() const;

        //- Return the liqFrac
        const scalar& liqFrac() const;

        //- Return the viscosity
        const scalar& viscosity() const;

        //- Return the minSep
        const scalar& minSep() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makePendularModel(CloudType)                                          \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(PendularModel<CloudType>, 0);         \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        PendularModel<CloudType>,                                             \
        dictionary                                                            \
    );


#define makePendularModelType(SS, CloudType)                                  \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType>, 0);                    \
                                                                              \
    PendularModel<CloudType>::adddictionaryConstructorToTable<SS<CloudType> > \
        add##SS##CloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PendularModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
