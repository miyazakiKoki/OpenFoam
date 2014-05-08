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

\*---------------------------------------------------------------------------*/

#include "CellAverageParticleVelocity.H"

// * * * * * * * * * * * * * Protectd Member Functions * * * * * * * * * * * //

template<class CloudType>
void Foam::CellAverageParticleVelocity<CloudType>::write()
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CellAverageParticleVelocity<CloudType>::CellAverageParticleVelocity
(
    const dictionary& dict,
    CloudType& owner
)
:
    CloudFunctionObject<CloudType>(owner),
    UpPtr_(NULL)
{}


template<class CloudType>
Foam::CellAverageParticleVelocity<CloudType>::CellAverageParticleVelocity
(
    const CellAverageParticleVelocity<CloudType>& vf
)
:
    CloudFunctionObject<CloudType>(vf),
    UpPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CellAverageParticleVelocity<CloudType>::~CellAverageParticleVelocity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CellAverageParticleVelocity<CloudType>::preEvolve()
{
    if (UpPtr_.valid())
    {
        UpPtr_->internalField() = vector::zero;
    }
    else
    {
        const fvMesh& mesh = this->owner().mesh();

        UpPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    this->owner().name() + "Up",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedVector("zeroVector", dimVelocity, vector::zero)
            )
        );
    }

    if (pVolPtr_.valid())
    {
        pVolPtr_->internalField() = scalar(0);
    }
    else
    {
        const fvMesh& mesh = this->owner().mesh();

        pVolPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "pVol",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zeroVolume", dimVolume, scalar(0))
            )
        );
    }
}


template<class CloudType>
void Foam::CellAverageParticleVelocity<CloudType>::postEvolve()
{
    volVectorField& U = UpPtr_();

    volScalarField& pVol = pVolPtr_();

    pVol.internalField() += scalar(VSMALL);

    const fvMesh& mesh = this->owner().mesh();

    U.internalField() /= mesh.time().deltaTValue()*pVol.internalField();

    CloudFunctionObject<CloudType>::postEvolve();
}


template<class CloudType>
void Foam::CellAverageParticleVelocity<CloudType>::postMove
(
    const parcelType& p,
    const label cellI,
    const scalar dt,
    const point& position0,
    bool&
)
{
    volVectorField& U = UpPtr_();

    volScalarField& pVol = pVolPtr_();

    scalar pv = p.mass()/p.rho();

    U[cellI] += dt*p.nParticle()*p.U()*pv;

    pVol[cellI] += pv;
}


// ************************************************************************* //
