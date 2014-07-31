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

#include "LambertWallNearestPoint.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::LambertWallNearestPoint<CloudType>::evaluatePendularWall
(
    typename CloudType::parcelType& p,
    const point& site,
    const WallSiteData<vector>& data,
    scalar pREff
) const
{
    const scalar& st = this->surfaceTension();
    const scalar& ca = this->contactAngle();
    const scalar& lf = this->liqFrac();
    const scalar& vis = this->viscosity();
    const scalar& ms = this->minSep();

    scalar Vtot = lf*(p.Vliq());

    vector r_PW = p.position() - site;

    vector U_PW = p.U() - data.wallData();

    scalar r_PW_mag = mag(r_PW);

    scalar normalOverlapMag = pREff - r_PW_mag;

    scalar S = -normalOverlapMag;

    vector rHat_PW = r_PW/(r_PW_mag + VSMALL);


    //The closest wet point forms liquid bridge with wall

    // calculate the minimum length between wet point and wall
    List<scalar> lengthParticleWall(p.liquidPositionVectors().size());

    forAll(lengthParticleWall,i)
    {    
        lengthParticleWall[i] = mag( site - p.liquidPositionVectors()[i] );
    }
    // find the minimum length from lengthParticleWall[i]
    // l is the number of array which the value of minLength is stored
    scalar minLength;
    minLength = lengthParticleWall[0];
    scalar l=0;
                
    forAll(lengthParticleWall,i)
    {
        if(lengthParticleWall[i] < minLength)
        {
            minLength = lengthParticleWall[i];
            l = i;
        }
    }

    //If the closest wet point is dry, then the liquid volume will be 0
    //so the force fN_PW and ft_PW will be 0
    //else, the volume of the liquid bridge will be p.partVliq[l]

    scalar Srup = (1+0.5*ca)*pow(p.partVliq()[l], 1./3.);

    if(p.partVliq()[l]!=0)
    {
        // Normal force
        scalar capMag =
            4*mathematical::pi*pREff*st*cos(ca)/
            (1+max(S, 0)*sqrt(mathematical::pi*pREff/p.partVliq()[l]));


        scalar Svis = max(pREff*ms, S);

        scalar etaN = 6*mathematical::pi*vis*pREff*pREff/Svis;

        vector fN_PW = (-capMag - etaN*(U_PW & rHat_PW)) * rHat_PW;

        p.f() += fN_PW;

        vector UT_PW = U_PW - (U_PW & rHat_PW)*rHat_PW
                      - ((pREff*p.omega()) ^ rHat_PW);

        scalar etaT =
            6*mathematical::pi*vis*pREff*(8./15.*log(pREff/Svis) + 0.9588);

        vector fT_PW = -etaT * UT_PW;

        p.f() += fT_PW;

        p.torque() += (pREff*-rHat_PW) ^ fT_PW;
        //Info << " particle is in contact with vector "<<p.liquidPositionVectors()[l] << endl;
        //Info << " The volume of this liquid bridge is "<< p.partVliq()[l] << endl;
    }





    //-----------------------------------------------------------------//
  //-------------------original------------------------------
    // Normal force
  /*  scalar capMag =
        4*mathematical::pi*pREff*st*cos(ca)/
        (1+max(S, 0)*sqrt(mathematical::pi*pREff/Vtot));

    Info << "the value of capMag is " << capMag << endl;
    Info << " the value of overlapMag S is " << S << endl;
    Info << " the volume of Vtot is " << Vtot << endl;

    scalar Svis = max(pREff*ms, S);

    scalar etaN = 6*mathematical::pi*vis*pREff*pREff/Svis;

    vector fN_PW = (-capMag - etaN*(U_PW & rHat_PW)) * rHat_PW;

    p.f() += fN_PW;

    vector UT_PW = U_PW - (U_PW & rHat_PW)*rHat_PW
                  - ((pREff*p.omega()) ^ rHat_PW);

    scalar etaT =
        6*mathematical::pi*vis*pREff*(8./15.*log(pREff/Svis) + 0.9588);

    vector fT_PW = -etaT * UT_PW;

    p.f() += fT_PW;

    p.torque() += (pREff*-rHat_PW) ^ fT_PW;*/
    //----------------------------------------------------------------

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LambertWallNearestPoint<CloudType>::LambertWallNearestPoint
(
    const dictionary& dict,
    CloudType& cloud,
    const scalar& surfaceTension,
    const scalar& contactAngle,
    const scalar& liqFrac,
    const scalar& viscosity,
    const scalar& minSep
)
:
    PendularWallModel<CloudType>
    (
        dict,
        cloud,
        typeName,
        surfaceTension,
        contactAngle,
        liqFrac,
        viscosity,
        minSep
    ),
    useEquivalentSize_(Switch(this->coeffDict().lookup("useEquivalentSize")))
{
    if (useEquivalentSize_)
    {
        volumeFactor_ = readScalar(this->coeffDict().lookup("volumeFactor"));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LambertWallNearestPoint<CloudType>::~LambertWallNearestPoint()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::LambertWallNearestPoint<CloudType>::pREff
(
    const typename CloudType::parcelType& p
) const
{
    if (useEquivalentSize_)
    {
        return p.d()/2*cbrt(p.nParticle()*volumeFactor_);
    }
    else
    {
        return p.d()/2;
    }
}

template<class CloudType>
void Foam::LambertWallNearestPoint<CloudType>::evaluatePendularWall
(
    typename CloudType::parcelType& p,
    const List<point>& flatSitePoints,
    const List<WallSiteData<vector> >& flatSiteData,
    const List<point>& sharpSitePoints,
    const List<WallSiteData<vector> >& sharpSiteData
) const
{
    scalar pREff = this->pREff(p);

    forAll(flatSitePoints, siteI)
    {
        evaluatePendularWall
        (
            p,
            flatSitePoints[siteI],
            flatSiteData[siteI],
            pREff
        );
    }

    forAll(sharpSitePoints, siteI)
    {
        // Treating sharp sites like flat sites

        evaluatePendularWall
        (
            p,
            sharpSitePoints[siteI],
            sharpSiteData[siteI],
            pREff
        );
    }

}



// ************************************************************************* //
