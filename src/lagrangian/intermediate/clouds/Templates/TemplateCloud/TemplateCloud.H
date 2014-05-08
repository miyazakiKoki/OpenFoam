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
    Foam::TemplateCloud

Description
    Adds coolisions to kinematic clouds

SourceFiles
    TemplateCloudI.H
    TemplateCloud.C

\*---------------------------------------------------------------------------*/

#ifndef TemplateCloud_H
#define TemplateCloud_H

#include "particle.H"
#include "Cloud.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class TemplateCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class TemplateCloud
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
        typedef TemplateCloud<CloudType> templateCloudType;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        TemplateCloud(const TemplateCloud&);

        //- Disallow default bitwise assignment
        void operator=(const TemplateCloud&);


public:

    // Constructors

        //- Construct given carrier gas fields
        TemplateCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g,
            bool readFields = true
        );

        //- Copy constructor with new name
        TemplateCloud
        (
            TemplateCloud<CloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        TemplateCloud
        (
            const fvMesh& mesh,
            const word& name,
            const TemplateCloud<CloudType>& c
        );


    //- Destructor
    virtual ~TemplateCloud();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "TemplateCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
