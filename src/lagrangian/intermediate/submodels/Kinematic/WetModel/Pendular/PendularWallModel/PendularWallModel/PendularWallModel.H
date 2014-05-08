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
    Foam::PendularWallModel

Description
    Templated wall interaction class

SourceFiles
    PendularWallModel.C
    PendularWallModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef PendularWallModel_H
#define PendularWallModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "WallSiteData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class PendularWallModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PendularWallModel
{
    // Private data

        //- The CollisionModel dictionary
        const dictionary& dict_;

        //- Reference to the owner cloud class
        CloudType& owner_;

        //- The coefficients dictionary
        const dictionary coeffDict_;

        //- Surface tension coefficient
        const scalar surfaceTension_;

        //- Contact angle
        const scalar contactAngle_;

        //- Fraction of the liquid to form one bridge
        const scalar liqFrac_;

        //- Viscosity
        const scalar viscosity_;

        //- Minimum separation distance divided by effective radius
        const scalar minSep_;

public:

    //- Runtime type information
    TypeName("pendularWallModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        PendularWallModel,
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
        (dict, owner, surfaceTension, contactAngle, liqFrac, viscosity, minSep)
    );


    // Constructors

        //- Construct from components
        PendularWallModel
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
        PendularWallModel
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
    virtual ~PendularWallModel();


    //- Selector
    static autoPtr<PendularWallModel<CloudType> > New
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

        //- Return non-const access to the owner cloud object
        CloudType& owner();

        //- Return the dictionary
        const dictionary& dict() const;

        //- Return the coefficients dictionary
        const dictionary& coeffDict() const;

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


    // Member Functions

        //- Return the effective radius for a particle for the model
        virtual scalar pREff(const typename CloudType::parcelType& p) const = 0;

        //- Calculate the wall interaction for a parcel
        virtual void evaluatePendularWall
        (
            typename CloudType::parcelType& p,
            const List<point>& flatSitePoints,
            const List<WallSiteData<vector> >& flatSiteData,
            const List<point>& sharpSitePoints,
            const List<WallSiteData<vector> >& sharpSiteData
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makePendularWallModel(CloudType)                                          \
                                                                                  \
    defineNamedTemplateTypeNameAndDebug(PendularWallModel<CloudType>, 0);         \
                                                                                  \
    defineTemplateRunTimeSelectionTable                                           \
    (                                                                             \
        PendularWallModel<CloudType>,                                             \
        dictionary                                                                \
    );


#define makePendularWallModelType(SS, CloudType)                                  \
                                                                                  \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType>, 0);                        \
                                                                                  \
    PendularWallModel<CloudType>::adddictionaryConstructorToTable<SS<CloudType> > \
        add##SS##CloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PendularWallModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
