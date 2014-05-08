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

#include "VoidFractionSubDivision.H"

// * * * * * * * * * * * * * Protectd Member Functions * * * * * * * * * * * //

template<class CloudType>
void Foam::VoidFractionSubDivision<CloudType>::write()
{
    if (thetaPtr_.valid())
    {
        thetaPtr_->write();
    }
    else
    {
        FatalErrorIn("void Foam::VoidFractionSubDivision<CloudType>::write()")
            << "thetaPtr not valid" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VoidFractionSubDivision<CloudType>::VoidFractionSubDivision
(
    const dictionary& dict,
    CloudType& owner
)
:
    CloudFunctionObject<CloudType>(owner),
    thetaPtr_(NULL),
    refinementLevel_
    (
        // Keep this as it is until new update so that coeffs can be specified within CloudFunctionObject scope
        readLabel(this->owner().particleProperties().subDict("voidFractionSubDivisionCoeffs").lookup("refinementLevel"))
    ),
    annulusThicknessFactor_
    (
        // Keep this as it is until new update so that coeffs can be specified within CloudFunctionObject scope
        readScalar(this->owner().particleProperties().subDict("voidFractionSubDivisionCoeffs").lookup("annulusThicknessFactor"))
    )
{}


template<class CloudType>
Foam::VoidFractionSubDivision<CloudType>::VoidFractionSubDivision
(
    const VoidFractionSubDivision<CloudType>& vf
)
:
    CloudFunctionObject<CloudType>(vf),
    thetaPtr_(NULL),
    refinementLevel_
    (
        // Keep this as it is until new update so that coeffs can be specified within CloudFunctionObject scope
        readLabel(this->owner().particleProperties().subDict("voidFractionSubDivisionCoeffs").lookup("refinementLevel"))
    ),
    annulusThicknessFactor_
    (
        // Keep this as it is until new update so that coeffs can be specified within CloudFunctionObject scope
        readScalar(this->owner().particleProperties().subDict("voidFractionSubDivisionCoeffs").lookup("annulusThicknessFactor"))
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::VoidFractionSubDivision<CloudType>::~VoidFractionSubDivision()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::VoidFractionSubDivision<CloudType>::preEvolve()
{
    if (thetaPtr_.valid())
    {
        thetaPtr_->internalField() = 0.0;
    }
    else
    {
        const fvMesh& mesh = this->owner().mesh();

        thetaPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "Theta",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
    }
}


template<class CloudType>
void Foam::VoidFractionSubDivision<CloudType>::postEvolve()
{
    volScalarField& theta = thetaPtr_();

    const fvMesh& mesh = this->owner().mesh();

    theta.internalField() /= mesh.time().deltaTValue()*mesh.V();

    CloudFunctionObject<CloudType>::postEvolve();
}


template<class CloudType>
void Foam::VoidFractionSubDivision<CloudType>::postMove
(
    const parcelType& p,
    const label cellI,
    const scalar dt,
    bool&
)
{
    labelHashSet lhash;

    recursivelyCalcVoidFraction
    (
        p,
        cellI,
        dt,
        lhash
    );

/*
    const scalar &d = p.d();
    const vector &position = p.position();

    scalar delRef = scalar(1) / refinementLevel_;

    scalar eleVol = pow(delRef * d, 3);

    forAll(elementPosition_,i)
    {
        vector elePos
        (
            position + d * elementPosition_[i]
        );

        // findCell is slow. Should be replaced with more efficient algorithm in the future.
        label eleCellI = this->owner().mesh().findCell(elePos);
        if(eleCellI >= 0)
        {
            theta[eleCellI] += dt*eleVol;
        }
    }
*/
    //theta[cellI] += dt*p.nParticle()*p.volume();
}


template<class CloudType>
void Foam::VoidFractionSubDivision<CloudType>::recursivelyCalcVoidFraction
(
    const parcelType& p,
    const label cellI,
    const scalar &dt,
    labelHashSet &lhash
)
{
    if(!lhash.found(cellI))
    {
        lhash.insert(cellI);

        const fvMesh &mesh(this->owner().mesh());

        volScalarField& theta = thetaPtr_();

        const scalar &d = p.d();

        const vector &position = p.position();

        const vector &cellPosition = mesh.C()[cellI];

        const scalar &cellVol = mesh.V()[cellI];

        scalar dist = mag(cellPosition - position);

        bool flagRecursive(false);

        if(dist < 0.5*d*(1 - annulusThicknessFactor_))
        {
            theta[cellI] += cellVol * dt;

            flagRecursive = true;
        }
        else if(dist < 0.5*d*(1 + annulusThicknessFactor_))
        {
            vector maxVertex = vector::min;

            vector minVertex = vector::max;

            const labelList& verticesID = mesh.cellPoints()[cellI];

            forAll(verticesID, i)
            {
                vector vertexPoint(mesh.points()[verticesID[i]]);

                maxVertex = max(maxVertex, vertexPoint);

                minVertex = min(minVertex, vertexPoint);
            }

            vector delRef = (maxVertex - minVertex) / refinementLevel_;

            scalar refVol = delRef.x() * delRef.y() * delRef.z();

            for(label ix = 0; ix < refinementLevel_; ix++)
            {
                for(label iy = 0; iy < refinementLevel_; iy++)
                {
                    for(label iz = 0; iz < refinementLevel_; iz++)
                    {
                        vector refinePosition
                        (
                            minVertex.x() + (ix+0.5) * delRef.x(),
                            minVertex.y() + (iy+0.5) * delRef.y(),
                            minVertex.z() + (iz+0.5) * delRef.z()
                        );

                        scalar distRef = mag(refinePosition - position);
                        if(distRef <= 0.5*d)
                        {
                            if(mesh.pointInCell(refinePosition, cellI))
                            {
                                theta[cellI] += refVol * dt;
                                flagRecursive = true;
                            }
                        }
                    }
                }
            }

        }

        if(flagRecursive)
        {
            const labelList &neighCells(mesh.cellCells()[cellI]);

            forAll(neighCells, i)
            {
                const label neighI = neighCells[i];

                recursivelyCalcVoidFraction
                (
                    p,
                    neighI,
                    dt,
                    lhash
                );
            }
        }
    }

}
// ************************************************************************* //
