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

#include "LambertDispersionNearestPoint.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::LambertDispersionNearestPoint<CloudType>::LambertDispersionNearestPoint
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
    PendularModel<CloudType>
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
Foam::LambertDispersionNearestPoint<CloudType>::~LambertDispersionNearestPoint()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::LambertDispersionNearestPoint<CloudType>::evaluatePendular
(
    typename CloudType::parcelType& pA,
    typename CloudType::parcelType& pB
) const
{
    const scalar& st = this->surfaceTension();
    const scalar& ca = this->contactAngle();
    const scalar& lf = this->liqFrac();
    const scalar& vis = this->viscosity();
    const scalar& ms = this->minSep();

    scalar Vtot = lf*(pA.Vliq() + pB.Vliq());

    if(Vtot > VSMALL)
    {
        //------------------my part-----------------------


        scalar pi = Foam::constant::mathematical::pi;

        vector r_AB = (pA.position() - pB.position());

        scalar dAEff = pA.d();

        scalar dBEff = pB.d();

        scalar r_AB_mag = mag(r_AB);

        scalar normalOverlapMag = 0.5*(dAEff + dBEff) - r_AB_mag;

        scalar S = -normalOverlapMag;


        // calculating the middle point between particle A and B
        vector middlePointAB;
        middlePointAB = ( pA.position() + pB.position() )/2;
        // calculate the minimum length between middle point and liquidPositions of particle A
        List<scalar> lengthMidPointLiqPosA(pA.liquidPositions().size());

        forAll(lengthMidPointLiqPosA,i)
        {
            lengthMidPointLiqPosA[i] = mag( middlePointAB - pA.liquidPositions()[i] );
        }

        // find the minimum length from lengthMidPointLiqPosA[i]
        // lA is the number of array which the value of minLengthA is stored
        scalar minLengthA;
        minLengthA = lengthMidPointLiqPosA[0];
        scalar lA=0;

        forAll(lengthMidPointLiqPosA,i)
        {
            if(lengthMidPointLiqPosA[i] < minLengthA)
            {
                minLengthA = lengthMidPointLiqPosA[i];
                lA = i;
            }
        }


        // calculate the minimum length between middle point and liquidPositions of particle B
        List<scalar> lengthMidPointLiqPosB(pB.liquidPositions().size());

        forAll(lengthMidPointLiqPosB,i)
        {
            lengthMidPointLiqPosB[i] = mag( middlePointAB - pB.liquidPositions()[i] );
        }

        // find the minimum length from lengthMidPointLiqPosB[i]
        // lB is the number of array which the value of minLengthB is stored
        scalar minLengthB;
        minLengthB = lengthMidPointLiqPosB[0];
        scalar lB=0;

        forAll(lengthMidPointLiqPosB,i)
        {
            if(lengthMidPointLiqPosB[i] < minLengthB)
            {
                minLengthB = lengthMidPointLiqPosB[i];
                lB = i;
            }
        }

        scalar VliqBrid;
        VliqBrid = pA.partVliq()[lA]+pB.partVliq()[lB];


        scalar Srup = (1+0.5*ca)*pow(VliqBrid, 1./3.);



        if (S < Srup)
        {

            vector rHat_AB = r_AB/(r_AB_mag + VSMALL);

            // Effective radius
            scalar R = 0.5*dAEff*dBEff/(dAEff + dBEff);
            
            
            // Normal force
            scalar capMag =
                4*mathematical::pi
                *R*st*cos(ca);
            


            if(VliqBrid > 0)
            {
            
                if(S > 0)
                {
                   capMag /= 1 + 1 /
                       (sqrt(1+VliqBrid/(mathematical::pi*R*S*S))-1);
                }



                // Relative translational velocity
                vector U_AB = pA.U() - pB.U();

                scalar Svis = max(R*ms, S);

                scalar etaN = 6*mathematical::pi*vis*R*R/Svis;

                vector fN_AB = (-capMag - etaN*(U_AB & rHat_AB)) * rHat_AB;

                pA.f() += fN_AB;
                pB.f() += -fN_AB;

                vector UT_AB = U_AB - (U_AB & rHat_AB)*rHat_AB
                         - ((0.5*(dAEff*pA.omega() + dBEff*pB.omega()))
                           ^ rHat_AB);

                scalar etaT =
                    6*mathematical::pi*vis*R*(8./15.*log(R/Svis) + 0.9588);

                vector fT_AB = -etaT * UT_AB;

                pA.f() += fT_AB;
                pB.f() += -fT_AB;

                pA.torque() += (dAEff/2*-rHat_AB) ^ fT_AB;
                pB.torque() += (dBEff/2*rHat_AB) ^ -fT_AB;



                //---------------------------liquid Distribution----------------------------------------
                //distributing liquid to each particle position
                //the amount of liquid distributed is VliqBrid/2
                pA.partVliq()[lA] = VliqBrid/2;
                pB.partVliq()[lB] = VliqBrid/2;


            }
            
            




        //------------------------------------------------//
        //------------------original washino sennsei model--------------------
            //Pendular bridge formed
/*
            vector rHat_AB = r_AB/(r_AB_mag + VSMALL);

            // Effective radius
            scalar R = 0.5*dAEff*dBEff/(dAEff + dBEff);

            // Normal force
            scalar capMag =
                4*mathematical::pi
                *R*st*cos(ca);

            if(S > 0)
            {
                capMag /= 1 + 1 /
                    (sqrt(1+Vtot/(mathematical::pi*R*S*S))-1);
            }

            // Relative translational velocity
            vector U_AB = pA.U() - pB.U();

            scalar Svis = max(R*ms, S);

            scalar etaN = 6*mathematical::pi*vis*R*R/Svis;

            vector fN_AB = (-capMag - etaN*(U_AB & rHat_AB)) * rHat_AB;

            pA.f() += fN_AB;
            pB.f() += -fN_AB;

            vector UT_AB = U_AB - (U_AB & rHat_AB)*rHat_AB;

            scalar etaT =
                6*mathematical::pi*vis*R*(8./15.*log(R/Svis) + 0.9588);

            vector fT_AB = -etaT * UT_AB;

            pA.f() += fT_AB;
            pB.f() += -fT_AB;

            pA.torque() += (dAEff/2*-rHat_AB) ^ fT_AB;
            pB.torque() += (dBEff/2*rHat_AB) ^ -fT_AB;
*/
        }
    }
}


// ************************************************************************* //
