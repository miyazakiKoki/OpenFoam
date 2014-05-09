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

#include "Lambert.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::Lambert<CloudType>::Lambert
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
Foam::Lambert<CloudType>::~Lambert()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::Lambert<CloudType>::evaluatePendular
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
    {/*
        //------------------my part-----------------------
        //find minimum distance of particle position
        label n = pA.liquidPositions().size()*pB.liquidPositions().size();
        scalarList liqPosDistance(n);
        forAll(pA.liquidPositions(),i)
        {
             forAll(pB.liquidPositions(),j)
             {
                 liqPosDistance[i+j]=mag(pA.liquidPositions()[i]-pB.liquidPositions()[j]);
             }
        }

        scalar liqPosMinDistance = 0;
        forAll(liqPosDistance,c)
        {
            if(liqPosDistance[c] < liqPosMinDistance)
            {
                liqPosMinDistance = liqPosDistance[c];
            }
        }
*/
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        //teimenn calculation//
        //height of droplet on the surface of particle
        scalar pi = Foam::constant::mathematical::pi;

        //radius of cap pA and pB
        List<scalar> pACapR(pA.liquidPositionVectors().size());

        forAll(pACapR,i)
        {
            pACapR[i]=pow(6*pA.partVliq()[i]/(pi*tan(ca/2)*(3+pow(tan(ca/2),2))),0.33333);
        }
        List<scalar> pBCapR(pB.liquidPositionVectors().size());

        forAll(pBCapR,i)
        {
            pBCapR[i]=pow(6*pB.partVliq()[i]/(pi*tan(ca/2)*(3+pow(tan(ca/2),2))),0.33333);
        }
        //calculating solid angle  cos(thetaC) pA and pB
        List<scalar> cosThetaCpA(pA.liquidPositionVectors().size());
        List<scalar> cosThetaCpB(pB.liquidPositionVectors().size());
        //calculating cos thetaCp
        forAll(cosThetaCpA,i)
        {
            cosThetaCpA[i]=sqrt((pA.d())*(pA.d())-pACapR[i]*pACapR[i])/(pA.d());
        }

        forAll(cosThetaCpB,i)
        {
            cosThetaCpB[i]=sqrt((pB.d())*(pB.d())-pBCapR[i]*pBCapR[i])/(pB.d() );
        }







        //---------------------------------------------------




        vector r_AB = (pA.position() - pB.position());

        scalar dAEff = pA.d();

        scalar dBEff = pB.d();

        scalar r_AB_mag = mag(r_AB);

        scalar normalOverlapMag = 0.5*(dAEff + dBEff) - r_AB_mag;

        scalar S = -normalOverlapMag;

        scalar Srup = (1+0.5*ca)*pow(Vtot, 1./3.);



        if (S < Srup)
        {

            //-----------------------my part------------------//
            //calculating angle between liquid position and r_AB
            List<scalar> cosThetaA(pA.liquidPositionVectors().size());
            List<scalar> cosThetaB(pB.liquidPositionVectors().size());

            forAll(cosThetaA,i)
            {
                cosThetaA[i]= (-r_AB/(r_AB_mag+VSMALL))&(pA.liquidPositionVectors()[i]/(mag(pA.liquidPositionVectors()[i])+VSMALL));
            }
            
            forAll(cosThetaB,i)
            {
                cosThetaB[i]= (r_AB/(r_AB_mag+VSMALL))&(pB.liquidPositionVectors()[i]/(mag(pB.liquidPositionVectors()[i])+VSMALL));
            }

            
            //----------- counter array for distributing liquid volume----//
            List<scalar> counterA(pA.liquidPositionVectors().size());
            List<scalar> counterB(pB.liquidPositionVectors().size());

            forAll(counterA,i)
            {
                counterA[i]=0;
            }
            forAll(counterB,i)
            {
                counterB[i]=0;
            }
            //determing which part on the surface of particles form liquid bridge.
            //calculating the volume of liquid used to form liquid bridge on each particle
            //give 1 to counter array when forming liquid bridge at i's position of array
            scalar VliqBrid;
            scalar VliqBridA;
            scalar VliqBridB;            

            VliqBrid = 0;
            VliqBridA = 0;
            VliqBridB = 0;
            forAll(cosThetaCpA,i)
            {
                if(cosThetaCpA[i]<cosThetaA[i])
                {
                    VliqBridA += pA.partVliq()[i];
                    counterA[i] = 1;
                }
            }
            forAll(cosThetaCpB,i)
            {
                if(cosThetaCpB[i]<cosThetaB[i])
                {
                    VliqBridB += pB.partVliq()[i];
                    counterB[i] = 1;
                }
            }
            // to check the value of counter
/*
            Info<<"counterA[0]="<<counterA[0]<<endl;
            Info<<"counterA[1]="<<counterA[1]<<endl;
            Info<<"counterA[2]="<<counterA[2]<<endl;
            Info<<"counterA[3]="<<counterA[3]<<endl;
            Info<<"counterA[4]="<<counterA[4]<<endl;
            Info<<"counterA[5]="<<counterA[5]<<endl;

            Info<<"counterB[0]="<<counterB[0]<<endl;
            Info<<"counterB[1]="<<counterB[1]<<endl;
            Info<<"counterB[2]="<<counterB[2]<<endl;
            Info<<"counterB[3]="<<counterB[3]<<endl;
            Info<<"counterB[4]="<<counterB[4]<<endl;
            Info<<"counterB[5]="<<counterB[5]<<endl;
*/


            //calculating total points needed to form liquid bridge
            scalar totalPoints;
            totalPoints = 0;
            forAll(counterA,i)
            {
                totalPoints += (counterA[i] + counterB[i]);
            }

//            Info << "totalPoints" << totalPoints << endl;
        

            //calculating the total liquid volume needed to form liquid bridge
            VliqBrid = VliqBridA + VliqBridB;

            vector rHat_AB = r_AB/(r_AB_mag + VSMALL);

            // Effective radius
            scalar R = 0.5*dAEff*dBEff/(dAEff + dBEff);
            
            
            // Normal force
            scalar capMag =
                4*mathematical::pi
                *R*st*cos(ca);
            
            scalar a = 1+VliqBrid/(mathematical::pi*R*S*S);
            scalar b = 1+Vtot/(mathematical::pi*R*S*S);
//            Info<<"inside my sqrt="<<a<<endl;
//            Info<<"inside original sqrt="<<b<<endl;
            if(VliqBrid > 0)
            {
            
                if(S > 0)
                {
                   capMag /= 1 + 1 /
                       (sqrt(1+VliqBrid/(mathematical::pi*R*S*S))-1);
                }
                /*
                if(S > 0)
                {
                    capMag /= 1 + 1 /
                        (sqrt(1+Vtot/(mathematical::pi*R*S*S))-1);
                }
                */
/*
                Info<<"aho------------------------------------------"<<endl;
                Info<<"VliqBridA"<<VliqBridA<<endl;
                Info<<"VliqBridB"<<VliqBridB<<endl;
                Info<<"VliqBrid"<<VliqBrid<<endl;
*/
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

                //distributing liquid to each particle position
                //liquid volume is distributed equally to each position
                //the amount of liquid distributed is VliqBrid/totalPoints





/*
                //-----------------------------------this is not right because particleA or B is competely dry
                List<scalar> dummy(6);
                forAll(dummy,i)
                {
                    dummy[i]=0;
                }
                

                if((pA.partVliq()==dummy)&&(pB.partVliq()!=dummy))
                {
                    Info<<"aho1--------------------------"<<endl;
                    
                    // calculating the middle point between particle A and B
                    vector middlePointAB;
                    middlePointAB = ( pA.position() + pB.position() )/2;

                    // defining new list to distribute liquid.
                    // particle A is comple
                    List<scalar> counterAdry(pA.liquidPositionVectors().size());
                    List<scalar> lengthMidPointLiqPos(pA.liquidPositionVectors().size());
                    
                    forAll(counterAdry,i)
                    {
                        counterAdry[i]=0;
                    }
                    
                    forAll(lengthMidPointLiqPos,i)
                    {    
                        lengthMidPointLiqPos[i] = mag( middlePointAB - pA.liquidPositionVectors()[i] );
                    }
                    scalar minLength;
                    minLength = lengthMidPointLiqPos[0];
                    forAll(lengthMidPointLiqPos,i)
                    {
                        if(lengthMidPointLiqPos[i] < minLength)
                        {
                            minLength = lengthMidPointLiqPos[i];
                        }
                    }
                    
                    forAll(lengthMidPointLiqPos,i)
                    {
                        if(lengthMidPointLiqPos[i]==minLength)
                        {
                            counterAdry[i] = 1;
                        }else
                        {
                            counterAdry[i] = 0;
                        }
                    }
                    

                    forAll(counterAdry,i)
                    {
                        if(counterAdry[i]==1)
                        {
                            pA.partVliq()[i]=VliqBrid/(totalPoints+1);
                        }
                    }

                    forAll(counterB,i)
                    {
                        if(counterB[i]==1)
                        {
                            pB.partVliq()[i]=VliqBrid/(totalPoints+1);
                        }
                    }
                                        
                    
                
                }
                else if((pB.partVliq()==dummy)&&(pA.partVliq()!=dummy))
                {
                    Info<<"aho2--------------------------"<<endl;
                }
                else if((pA.partVliq()!=dummy)&&(pB.partVliq()!=dummy))
                {
                    Info<<"aho3--------------------------"<<endl;

                    forAll(counterA,i)
                    {
                        if(counterA[i]==1)
                        {
                            pA.partVliq()[i]=VliqBrid/totalPoints;
                        }
                    }

                    forAll(counterB,i)
                    {
                        if(counterB[i]==1)
                        {
                            pB.partVliq()[i]=VliqBrid/totalPoints;
                        }
                    }

                }*/
                 //-------------------------------------------------//


                //-----probably the correct liquid transportation method----
                // calculating the middle point between particle A and B
                vector middlePointAB;
                middlePointAB = ( pA.position() + pB.position() )/2;
                // calculate the minimum length between middle point and liquidPositionVectors of particle A
                List<scalar> lengthMidPointLiqPosA(pA.liquidPositionVectors().size());

                forAll(lengthMidPointLiqPosA,i)
                {    
                    lengthMidPointLiqPosA[i] = mag( middlePointAB - pA.liquidPositionVectors()[i] );
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


                // calculate the minimum length between middle point and liquidPositionVectors of particle B
                List<scalar> lengthMidPointLiqPosB(pB.liquidPositionVectors().size());

                forAll(lengthMidPointLiqPosB,i)
                {    
                    lengthMidPointLiqPosB[i] = mag( middlePointAB - pB.liquidPositionVectors()[i] );
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

                // liquid transportation routine
                // 3 patterns.

                if ( (pA.partVliq()[lA] == 0)&&(pB.partVliq()[lB] != 0) )
                {
                    //Info<<"aho1--------------------------"<<endl;
                    
                    pA.partVliq()[lA] = VliqBrid/(totalPoints+1);

                    forAll(counterB,i)
                    {
                        if(counterB[i]==1)
                        {
                            pB.partVliq()[i]=VliqBrid/(totalPoints+1);
                        }
                    }

                }
                else if ( (pA.partVliq()[lA] != 0)&&(pB.partVliq()[lB] == 0) )
                {
                    //Info<<"aho2--------------------------"<<endl;

                    pB.partVliq()[lB] = VliqBrid/(totalPoints+1);

                    forAll(counterA,i)
                    {
                        if(counterA[i]==1)
                        {
                            pA.partVliq()[i]=VliqBrid/(totalPoints+1);
                        }
                    }
                }
                else if ( (pA.partVliq()[lA] != 0)&&(pB.partVliq()[lB] != 0) )
                {
                    //Info<<"aho3--------------------------"<<endl;
                    forAll(counterA,i)
                    {
                        if(counterA[i]==1)
                        {
                            pA.partVliq()[i]=VliqBrid/totalPoints;
                        }
                    }

                    forAll(counterB,i)
                    {
                        if(counterB[i]==1)
                        {
                            pB.partVliq()[i]=VliqBrid/totalPoints;
                        }
                    }


                }
  


                

                

/*
                forAll(counterA,i)
                {
                    if(counterA[i]==1)
                        {
                            pA.partVliq()[i]=VliqBrid/totalPoints;
                        }
                }

                forAll(counterB,i)
                {
                    if(counterB[i]==1)
                        {
                            pB.partVliq()[i]=VliqBrid/totalPoints;
                        }
                }
                */


            }else
            {
                /* to exclude VligBrid = 0, which causes error */
             /*nothing is done here. Same condition as none pendular bridge, by giving 0 to force */
                vector U_AB = pA.U() - pB.U();
                vector UT_AB = U_AB - (U_AB & rHat_AB)*rHat_AB;
                vector fT_AB = 0 * UT_AB;
                pA.f() += fT_AB;
                pB.f() += -fT_AB;

                pA.torque() += (dAEff/2*-rHat_AB) ^ fT_AB;
                pB.torque() += (dBEff/2*rHat_AB) ^ -fT_AB;
            }
            
            




        //------------------------------------------------//

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
