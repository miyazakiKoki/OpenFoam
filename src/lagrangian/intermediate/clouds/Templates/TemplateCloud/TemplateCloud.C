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

\*---------------------------------------------------------------------------*/

#include "TemplateCloud.H"
#include "CollisionModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TemplateCloud<CloudType>::TemplateCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g,
    bool readFields
)
:
    CloudType(cloudName, rho, U, mu, g, false)
{
    if (this->solution().steadyState())
    {
        FatalErrorIn
        (
            "Foam::TemplateCloud<CloudType>::TemplateCloud"
            "("
                "const word&, "
                "const volScalarField&, "
                "const volVectorField&, "
                "const volScalarField&, "
                "const dimensionedVector&, "
                "bool"
            ")"
        )   << "Collision modelling not currently available for steady state "
            << "calculations" << exit(FatalError);
    }

    if (this->solution().active())
    {
        if (readFields)
        {
            parcelType::readFields(*this);
        }
    }

}


template<class CloudType>
Foam::TemplateCloud<CloudType>::TemplateCloud
(
    TemplateCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name)
{}


template<class CloudType>
Foam::TemplateCloud<CloudType>::TemplateCloud
(
    const fvMesh& mesh,
    const word& name,
    const TemplateCloud<CloudType>& c
)
:
    CloudType(mesh, name, c)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::TemplateCloud<CloudType>::~TemplateCloud()
{}
