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

#include "WetRegion.H"

// * * * * * * * * * * * * * Protectd Member Functions * * * * * * * * * * * //

template<class CloudType>
void Foam::WetRegion<CloudType>::write()
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::WetRegion<CloudType>::WetRegion
(
    const dictionary& dict,
    CloudType& owner
)
:
    CloudFunctionObject<CloudType>(dict, owner, typeName),
    UpPtr_(NULL),
    bb_
    (
        this->coeffDict().lookup("minBb"),
        this->coeffDict().lookup("maxBb")
    ),
    liquidAddRate_
    (
        readScalar(this->coeffDict().lookup("LiquidAdditionRate"))
    ),
    SOI_
    (
        readScalar(this->coeffDict().lookup("SOI"))
    ),
    EOI_
    (
        readScalar(this->coeffDict().lookup("EOI"))
    )
{
    if (EOI_ < VSMALL)
    {
        EOI_ = GREAT;
    }
}


template<class CloudType>
Foam::WetRegion<CloudType>::WetRegion
(
    const WetRegion<CloudType>& vf
)
:
    CloudFunctionObject<CloudType>(vf),
    UpPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::WetRegion<CloudType>::~WetRegion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::WetRegion<CloudType>::preEvolve()
{
}


template<class CloudType>
void Foam::WetRegion<CloudType>::postEvolve()
{
    
    scalar t(this->owner().mesh().time().value());

    if (t > SOI_ && t < EOI_)
    {
        label Np = 0;

        scalar VliqTotal = 0;

        forAllIter(typename CloudType, this->owner(), iter)
        {
            typename CloudType::parcelType& p = iter();

            if (bb_.contains(p.position()))
            {
                Np++;
            }

            VliqTotal += p.Vliq();
        }
        //Info << "VliqTotal is " << VliqTotal<<endl;
        reduce(Np, sumOp<label>());
        reduce(VliqTotal, sumOp<scalar>());

        scalar Vtarget = (t - SOI_)*liquidAddRate_;

        scalar addedVolume = Vtarget - VliqTotal;
        //Info << "addedVolume is " << addedVolume << endl;

        if (addedVolume < 0)
        {
            Info << "Warning: liquid addition went negative.\n" << endl;
        }

        if (Np == 0)
        {
            Info << "Warning: No particle is in the wetting region.\n" << endl;
        }
        
        //Info << "number of Np is " <<Np << endl; 

        forAllIter(typename CloudType, this->owner(), iter)
        {
            typename CloudType::parcelType& p = iter();

            if (bb_.contains(p.position()))
            {
                //p.Vliq() += addedVolume / Np;
                
                forAll(p.partVliq(),j)
                {
                    p.partVliq()[j] += addedVolume / (Np*p.partVliq().size());                  
                }
                
                p.Vliq() = 0;
                
                forAll(p.partVliq(),i)
                {
                    p.Vliq() += p.partVliq()[i];
                }
                
                
            }
        }

    }
    

    CloudFunctionObject<CloudType>::postEvolve();
}


template<class CloudType>
void Foam::WetRegion<CloudType>::postMove
(
    const parcelType& p,
    const label cellI,
    const scalar dt,
    const point& position0,
    bool&
)
{
}


// ************************************************************************* //
