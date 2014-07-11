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

#include "LambertWall.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::LambertWall<CloudType>::evaluatePendularWall
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


    // Normal force
    scalar capMag =
        4*mathematical::pi*pREff*st*cos(ca)/
        (1+max(S, 0)*sqrt(mathematical::pi*pREff/Vtot));

    //Info << "the value of capMag is " << capMag << endl;
    //Info << " the value of overlapMag S is " << S << endl;
    //Info << " the volume of Vtot is " << Vtot << endl;

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

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LambertWall<CloudType>::LambertWall
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
Foam::LambertWall<CloudType>::~LambertWall()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::LambertWall<CloudType>::pREff
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
void Foam::LambertWall<CloudType>::evaluatePendularWall
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
