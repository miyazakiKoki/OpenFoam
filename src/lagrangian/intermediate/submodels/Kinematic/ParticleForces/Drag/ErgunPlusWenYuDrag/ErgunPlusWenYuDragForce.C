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

#include "ErgunPlusWenYuDragForce.H"
#include "GeometricFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::ErgunPlusWenYuDragForce<CloudType>::Cd(const scalar Re) const
{
    if (Re > 1000.0)
    {
        return 0.424;
    }
    else
    {
        return 24.0*(1.0 + 1.0/6.0*pow(Re, 2.0/3.0))/(Re+ROOTVSMALL);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ErgunPlusWenYuDragForce<CloudType>::ErgunPlusWenYuDragForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, false)
{}


template<class CloudType>
Foam::ErgunPlusWenYuDragForce<CloudType>::ErgunPlusWenYuDragForce
(
    const ErgunPlusWenYuDragForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ErgunPlusWenYuDragForce<CloudType>::~ErgunPlusWenYuDragForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::ErgunPlusWenYuDragForce<CloudType>::calcCoupled
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
    const volScalarField& vf =
        ms.lookupObject<volScalarField>("voidFraction");
    const scalar& vfp = vf.internalField()[p.cell()];

    value.Sp() = beta(p,Re,muc,vfp)*mass/(p.rho()*(1-vfp+ROOTVSMALL));

    return value;
}


template<class CloudType>
Foam::scalar Foam::ErgunPlusWenYuDragForce<CloudType>::beta
(
    const typename CloudType::parcelType& p,
    const scalar Re,
    const scalar muc,
    const scalar vfp
) const
{
    if (vfp <= 0.8)
    {
        return ErgunBeta(p,Re*vfp,muc,vfp);
    }
    else
    {
        return WenYuBeta(p,Re*vfp,muc,vfp);
    }
}


template<class CloudType>
Foam::scalar Foam::ErgunPlusWenYuDragForce<CloudType>::ErgunBeta
(
    const typename CloudType::parcelType& p,
    const scalar Rep,
    const scalar muc,
    const scalar vfp
) const
{
    return muc*(1-vfp)*(150*(1-vfp)+1.75*Rep)/(sqr(p.d())*vfp+ROOTVSMALL);
}


template<class CloudType>
Foam::scalar Foam::ErgunPlusWenYuDragForce<CloudType>::WenYuBeta
(
    const typename CloudType::parcelType& p,
    const scalar Rep,
    const scalar muc,
    const scalar vfp
) const
{
    return 0.75*Cd(Rep)*muc*(1-vfp)*(pow(vfp,-2.7)+ROOTVSMALL)*Rep/sqr(p.d());
}
// ************************************************************************* //
