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

#include "WetCloud.H"
#include "WetModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::WetCloud<CloudType>::setModels()
{
    wetModel_.reset
    (
        WetModel<WetCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}

template<class CloudType>
template<class TrackData>
void  Foam::WetCloud<CloudType>::moveWetCollide
(
    TrackData& td,
    const scalar deltaT
)
{
    td.part() = TrackData::tpVelocityHalfStep;
    CloudType::move(td,  deltaT);

    td.part() = TrackData::tpLinearTrack;
    CloudType::move(td,  deltaT);

    // td.part() = TrackData::tpRotationalTrack;
    // CloudType::move(td);

    this->updateCellOccupancy();

    this->collision().collide();

    this->wetModel().bond();

    td.part() = TrackData::tpVelocityHalfStep;
    CloudType::move(td,  deltaT);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::WetCloud<CloudType>::WetCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& mu,
    const dimensionedVector& g,
    bool readFields
)
:
    CloudType(cloudName, rho, U, mu, g, false),
    wetModel_(NULL)
{
    Info << "boke-------------" << endl;
    if (this->solution().steadyState())
    {
        FatalErrorIn
        (
            "Foam::WetCloud<CloudType>::WetCloud"
            "("
                "const word&, "
                "const volScalarField&, "
                "const volVectorField&, "
                "const volScalarField&, "
                "const dimensionedVector&, "
                "bool"
            ")"
        )   << "Wet collision modelling not currently available for steady state "
            << "calculations" << exit(FatalError);
    }

    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this);
        }
    }

    Info<<"aho-------------------------------"<<endl;

}


template<class CloudType>
Foam::WetCloud<CloudType>::WetCloud
(
    WetCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    wetModel_(c.wetModel_->clone())
{}


template<class CloudType>
Foam::WetCloud<CloudType>::WetCloud
(
    const fvMesh& mesh,
    const word& name,
    const WetCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    wetModel_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::WetCloud<CloudType>::~WetCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::WetCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::template
            TrackingData<WetCloud<CloudType> > td(*this);

        this->solve(td);
    }
}


template<class CloudType>
template<class TrackData>
void  Foam::WetCloud<CloudType>::motion(TrackData& td)
{
    // Sympletic leapfrog integration of particle forces:
    // + apply half deltaV with stored force
    // + move positions with new velocity
    // + calculate forces in new position
    // + apply half deltaV with new force

    label nSubCycles = this->collision().nSubCycles();

    if (nSubCycles > 1)
    {
        Info<< "    " << nSubCycles << " move-wet-collide subCycles" << endl;

        subCycleTime moveWetCollideSubCycle
        (
            const_cast<Time&>(this->db().time()),
            nSubCycles
        );

        while(!(++moveWetCollideSubCycle).end())
        {
            moveWetCollide(td, this->db().time().deltaTValue());
        }

        moveWetCollideSubCycle.endSubCycle();
    }
    else
    {
        moveWetCollide(td, this->db().time().deltaTValue());
    }
}


// ************************************************************************* //
