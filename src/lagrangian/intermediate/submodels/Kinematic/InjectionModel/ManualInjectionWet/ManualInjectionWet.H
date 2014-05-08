/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::ManualInjectionWet

Description
    Manual injection
    - User specifies
      - Total mass to inject
      - Parcel positions in file \c positionsFile
      - Initial parcel velocity
    - Parcel diameters obtained by distribution model
    - All parcels introduced at SOI

SourceFiles
    ManualInjectionWet.C

\*---------------------------------------------------------------------------*/

#ifndef ManualInjectionWet_H
#define ManualInjectionWet_H

#include "InjectionModel.H"
#include "distributionModel.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ManualInjectionWet Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ManualInjectionWet
:
    public InjectionModel<CloudType>
{
    // Private data

        //- Name of file containing positions data
        const word positionsFile_;

        //- Field of parcel positions
        vectorIOField positions_;

        //- Field of parcel diameters
        scalarList diameters_;

        //- List of cell labels corresponding to injector positions
        labelList injectorCells_;

        //- List of tetFace labels corresponding to injector positions
        labelList injectorTetFaces_;

        //- List of tetPt labels corresponding to injector positions
        labelList injectorTetPts_;

        //- Initial parcel velocity
        const vector U0_;

        //- Parcel size distribution model
        const autoPtr<distributionModels::distributionModel> sizeDistribution_;

        //- Flag to suppress errors if particle injection site is out-of-bounds
        Switch ignoreOutOfBounds_;

        //- Initial liquid volume per each particle
        const scalar iniVliq_;

        //- Initial liquid positions for each particle
        const List<vector> iniLiquidPositions_;

        //- min vector of wet position
        const vector minWet_;

        //- max vector of wet position
        const vector maxWet_;


public:

    //- Runtime type information
    TypeName("manualInjectionWet");


    // Constructors

        //- Construct from dictionary
        ManualInjectionWet
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        ManualInjectionWet(const ManualInjectionWet<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType> > clone() const
        {
            return autoPtr<InjectionModel<CloudType> >
            (
                new ManualInjectionWet<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ManualInjectionWet();


    // Member Functions

        //- Set injector locations when mesh is updated
        virtual void updateMesh();

        //- Return the end-of-injection time
        scalar timeEnd() const;

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);


        // Injection geometry

            //- Set the injection position and owner cell, tetFace and tetPt
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner,
                label& tetFaceI,
                label& tetPtI
            );

            //- Set the parcel properties
            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;

            //- Return flag to identify whether or not injection of parcelI is
            //  permitted
            virtual bool validInjection(const label parcelI);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ManualInjectionWet.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
