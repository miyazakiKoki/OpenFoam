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
    Foam::WetCloud

Description
    Adds coolisions to kinematic clouds

SourceFiles
    WetCloudI.H
    WetCloud.C

\*---------------------------------------------------------------------------*/

#ifndef WetCloud_H
#define WetCloud_H

#include "particle.H"
#include "Cloud.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class WetModel;

/*---------------------------------------------------------------------------*\
                       Class WetCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class WetCloud
:
    public CloudType
{
public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef WetCloud<CloudType> wetCloudType;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        WetCloud(const WetCloud&);

        //- Disallow default bitwise assignment
        void operator=(const WetCloud&);

protected:
    // Protected data

        // References to the cloud sub-models

            //- Wet model
            autoPtr<WetModel<WetCloud<CloudType> > > wetModel_;

            // Initialisation

                //- Set cloud sub-models
                void setModels();

            // Cloud evolution functions

                //- Move-wet particles
                template<class TrackData>
                void moveWetCollide(TrackData& td, const scalar deltaT);

public:

    // Constructors

        //- Construct given carrier gas fields
        WetCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g,
            bool readFields = true
        );

        //- Copy constructor with new name
        WetCloud
        (
            WetCloud<CloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        WetCloud
        (
            const fvMesh& mesh,
            const word& name,
            const WetCloud<CloudType>& c
        );


    //- Destructor
    virtual ~WetCloud();

    // Member Functions

        // Access

            // Sub-models

                //- Return const access to the wet model
                inline const WetModel<WetCloud<CloudType> >&
                    wetModel() const;

                //- Return reference to the wet model
                inline WetModel<WetCloud<CloudType> >&
                    wetModel();

        // Evolution

            //- Evolve the cloud
            void evolve();

            //- Particle motion
            template<class TrackData>
            void motion(TrackData& td);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "WetCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WetCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
