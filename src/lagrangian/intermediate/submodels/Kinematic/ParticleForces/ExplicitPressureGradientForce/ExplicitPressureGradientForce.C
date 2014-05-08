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

#include "ExplicitPressureGradientForce.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ExplicitPressureGradientForce<CloudType>::ExplicitPressureGradientForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, false)
{
    const fvMesh& ms = this->mesh();

    const volScalarField& p
        = ms.lookupObject<volScalarField>("p");

    if(p.dimensions() == dimPressure)
    {
        isDividedByRho_ = false;
    }
    else
    {
        isDividedByRho_ = true;
    }

}


template<class CloudType>
Foam::ExplicitPressureGradientForce<CloudType>::ExplicitPressureGradientForce(const ExplicitPressureGradientForce& gf)
:
    ParticleForce<CloudType>(gf)
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ExplicitPressureGradientForce<CloudType>::~ExplicitPressureGradientForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::ExplicitPressureGradientForce<CloudType>::calcNonCoupled
(
    const typename CloudType::parcelType& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
    forceSuSp value(vector::zero, 0.0);

    const fvMesh& ms = this->mesh();

    const volVectorField& gradP
        = ms.lookupObject<volVectorField>("gradP");

    const vector& gp = gradP.internalField()[p.cell()];

    scalar rhoFactor = 1;
    if(isDividedByRho_)
    {
        rhoFactor = p.rhoc();
    }

    value.Su() = -mass*rhoFactor*gp/p.rho();

    return value;
}


// ************************************************************************* //
