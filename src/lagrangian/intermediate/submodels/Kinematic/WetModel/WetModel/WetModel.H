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
    Foam::WetModel

Description
    Templated wet model class.

SourceFiles
    WetModel.C
    WetModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef WetModel_H
#define WetModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "SubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class WetModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class WetModel
:
    public SubModelBase<CloudType>
{
    // Protected data

        //- Convenience typedef for parcel type
        typedef typename CloudType::parcelType parcelType;


public:

    //- Runtime type information
    TypeName("wetModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        WetModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        WetModel(CloudType& owner);

        //- Construct from components
        WetModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        WetModel(WetModel<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<WetModel<CloudType> > clone()
        {
            return autoPtr<WetModel<CloudType> >
            (
                new WetModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~WetModel();


    //- Selector
    static autoPtr<WetModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        //- Indicates whether model determines wall bonding or not,
        //  used to determine what value to use for wallImpactDistance
        virtual bool controlsWallInteraction() const;

        // Bonding function
        virtual void bond();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeWetModel(CloudType)                                         \
                                                                              \
    typedef CloudType::wetCloudType wetCloudType;                 \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        WetModel<wetCloudType>,                                   \
        0                                                                     \
    );                                                                        \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        WetModel<wetCloudType>,                                   \
        dictionary                                                            \
    );


#define makeWetModelType(SS, CloudType)                                 \
                                                                              \
    typedef CloudType::wetCloudType wetCloudType;                 \
    defineNamedTemplateTypeNameAndDebug(SS<wetCloudType>, 0);           \
                                                                              \
    WetModel<wetCloudType>::                                      \
        adddictionaryConstructorToTable<SS<wetCloudType> >              \
            add##SS##CloudType##wetCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WetModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
