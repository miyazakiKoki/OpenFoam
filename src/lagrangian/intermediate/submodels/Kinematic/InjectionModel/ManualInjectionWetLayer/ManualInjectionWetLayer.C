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

#include "ManualInjectionWetLayer.H"
#include "mathematicalConstants.H"
#include "PackedBoolList.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ManualInjectionWetLayer<CloudType>::ManualInjectionWetLayer
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    positionsFile_(this->coeffDict().lookup("positionsFile")),
    positions_
    (
        IOobject
        (
            positionsFile_,
            owner.db().time().constant(),
            owner.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    diameters_(positions_.size()),
    injectorCells_(positions_.size(), -1),
    injectorTetFaces_(positions_.size(), -1),
    injectorTetPts_(positions_.size(), -1),
    U0_(this->coeffDict().lookup("U0")),
    sizeDistribution_
    (
        distributionModels::distributionModel::New
        (
            this->coeffDict().subDict("sizeDistribution"),
            owner.rndGen()
        )
    ),
    ignoreOutOfBounds_
    (
        this->coeffDict().lookupOrDefault("ignoreOutOfBounds", false)
    ),
    iniVliq_
    (
        readScalar
        (
            this->coeffDict().lookup("initialLiquidVolumePerParticle")
        )
    ),
    iniLiquidPositions_
    (
        (
            this->coeffDict().lookup("initialLiquidPositionVector")
        )
    ),
    minWet_
    (
        (
            this->coeffDict().lookup("minWetVector")
        )
    ),
    maxWet_
    (
        (
            this->coeffDict().lookup("maxWetVector")
        )
    ),
    layerAxis_
    (
        (
            this->coeffDict().lookup("layerAxis")
            )
    ),
    layerNum_
    (
        readLabel
        (
            this->coeffDict().lookup("layerNum")
        )
    )
{
    updateMesh();

    // Construct parcel diameters
    forAll(diameters_, i)
    {
        diameters_[i] = sizeDistribution_->sample();
    }

    // Determine volume of particles to inject
    this->volumeTotal_ = sum(pow3(diameters_))*pi/6.0;
}


template<class CloudType>
Foam::ManualInjectionWetLayer<CloudType>::ManualInjectionWetLayer
(
    const ManualInjectionWetLayer<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    positionsFile_(im.positionsFile_),
    positions_(im.positions_),
    diameters_(im.diameters_),
    injectorCells_(im.injectorCells_),
    injectorTetFaces_(im.injectorTetFaces_),
    injectorTetPts_(im.injectorTetPts_),
    U0_(im.U0_),
    sizeDistribution_(im.sizeDistribution_().clone().ptr()),
    ignoreOutOfBounds_(im.ignoreOutOfBounds_),
    iniVliq_(im.iniVliq_),
    iniLiquidPositions_(im.iniLiquidPositions_),
    minWet_(im.minWet_),
    maxWet_(im.maxWet_),
    layerAxis_(im.layerAxis_),
    layerNum_(im.layerNum_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ManualInjectionWetLayer<CloudType>::~ManualInjectionWetLayer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ManualInjectionWetLayer<CloudType>::updateMesh()
{
    label nRejected = 0;

    PackedBoolList keep(positions_.size(), true);

    forAll(positions_, pI)
    {
        if
        (
            !this->findCellAtPosition
            (
                injectorCells_[pI],
                injectorTetFaces_[pI],
                injectorTetPts_[pI],
                positions_[pI],
                !ignoreOutOfBounds_
            )
        )
        {
            keep[pI] = false;
            nRejected++;
        }
    }

    if (nRejected > 0)
    {
        inplaceSubset(keep, positions_);
        inplaceSubset(keep, diameters_);
        inplaceSubset(keep, injectorCells_);
        inplaceSubset(keep, injectorTetFaces_);
        inplaceSubset(keep, injectorTetPts_);

        Info<< "    " << nRejected
            << " particles ignored, out of bounds" << endl;
    }
}


template<class CloudType>
Foam::scalar Foam::ManualInjectionWetLayer<CloudType>::timeEnd() const
{
    // Not used
    return this->SOI_;
}


template<class CloudType>
Foam::label Foam::ManualInjectionWetLayer<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((0.0 >= time0) && (0.0 < time1))
    {
        return positions_.size();
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::ManualInjectionWetLayer<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    // All parcels introduced at SOI
    if ((0.0 >= time0) && (0.0 < time1))
    {
        return this->volumeTotal_;
    }
    else
    {
        return 0.0;
    }
}


template<class CloudType>
void Foam::ManualInjectionWetLayer<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label,
    const scalar,
    vector& position,
    label& cellOwner,
    label& tetFaceI,
    label& tetPtI
)
{
    position = positions_[parcelI];
    cellOwner = injectorCells_[parcelI];
    tetFaceI = injectorTetFaces_[parcelI];
    tetPtI = injectorTetPts_[parcelI];
}


template<class CloudType>
void Foam::ManualInjectionWetLayer<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{

    // set particle velocity
    parcel.U() = U0_;

    parcel.liquidPositionVectors().setSize(iniLiquidPositions_.size()); 

    // set particle diameter
    parcel.d() = diameters_[parcelI];
    
    List<scalar> a(parcel.liquidPositionVectors().size());
    forAll(a,i)
    {
        a[i]=iniVliq_/parcel.liquidPositionVectors().size();
    }
    List<scalar> b(parcel.liquidPositionVectors().size());
    forAll(b,i)
    {
        b[i]=0;
    }

    // set initial liquid volume
    //parcel.partVliq() = iniVliq_/parcel.liquidPositionVectors().size();
    //parcel.partVliq() = a;

    // set initial liquid position vector

    parcel.liquidPositionVectors() = iniLiquidPositions_;

    //-------------------specify wet area------------------
    // give initial liquid to particles that are in the area between wetMin and wetMax
    /*if((minWet_.x()<=parcel.position().x()<=maxWet_.x())&&(minWet_.y()<=parcel.position().y()<=maxWet_.y())&&(minWet_.z()<=parcel.position().z()<=maxWet_.z()))
    {
        parcel.partVliq() = a;
    }
    else{
        parcel.partVliq() = b;
    }*/


        vector x (1,0,0);
        vector y (0,1,0);
        vector z (0,0,1);

        //When using this class, maxWet and minWet of the layerAxis must be set exactly at the point which particles exist.
        //If, for example the minWet.x and minWet., where set enough or too large, the only the partial of particles would be layer

        if(layerAxis_ == x)
        {
        	label i;
        	scalar layerLength = (maxWet_.x() - minWet_.x() )/layerNum_;
        	for (i=1; i<=layerNum_;i++ )
        	{
        	    if(((minWet_.x()+layerLength*(i-1))<=parcel.position().x())&&((parcel.position().x()<=minWet_.x()+layerLength*i))&&(minWet_.y()<=parcel.position().y())&&(parcel.position().y()<=maxWet_.y())&&(minWet_.z()<=parcel.position().z())&&(parcel.position().z()<=maxWet_.z()))
        	    {
            	    parcel.partVliq() = a/i;
        	    }
        	}
        	if (parcel.partVliq() < a/(layerNum_+1))
        	{
        		parcel.partVliq() = b;
        	}
        }else if(layerAxis_ == y)
        {
        	label i;
            scalar layerLength = (maxWet_.y() - minWet_.y() )/layerNum_;
            for (i=1; i<=layerNum_;i++ )
            {
        	    if((minWet_.x()<=parcel.position().x())&&(parcel.position().x()<=maxWet_.x())&&(minWet_.y()+layerLength*(i-1)<=parcel.position().y())&&(parcel.position().y()<=minWet_.y()+layerLength*i)&&(minWet_.z()<=parcel.position().z())&&(parcel.position().z()<=maxWet_.z()))
        	    {
        	        parcel.partVliq() = a/i;
        	    }
        	}
        	if (parcel.partVliq() < a/(layerNum_+1))
        	{
        	    parcel.partVliq() = b;
        	}
        }else if(layerAxis_ == z)
        {
        	label i;
            scalar layerLength = (maxWet_.z() - minWet_.z() )/layerNum_;
            for (i=1; i<=layerNum_;i++ )
            {
        	    if((minWet_.x()<=parcel.position().x())&&(parcel.position().x()<=maxWet_.x())&&(minWet_.y()<=parcel.position().y())&&(parcel.position().y()<=minWet_.y())&&(minWet_.z()+layerLength*(i-1)<=parcel.position().z())&&(parcel.position().z()<=maxWet_.z()+layerLength*i))
        	    {
        	        parcel.partVliq() = a/i;
        	    }
        	}
        	if (parcel.partVliq() < a/(layerNum_+1))
        	{
        	    parcel.partVliq() = b;
        	}
        }

}


template<class CloudType>
bool Foam::ManualInjectionWetLayer<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::ManualInjectionWetLayer<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
