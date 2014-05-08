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

#include "PendularWallModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::autoPtr<Foam::PendularWallModel<CloudType> >
Foam::PendularWallModel<CloudType>::New
(
    const dictionary& dict,
    CloudType& owner,
    const scalar& surfaceTension,
    const scalar& contactAngle,
    const scalar& liqFrac,
    const scalar& viscosity,
    const scalar& minSep
)
{
    word PendularWallModelType(dict.lookup("pendularWallModel"));

    Info<< "Selecting pendular wall model " << PendularWallModelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(PendularWallModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "PendularWallModel<CloudType>::New"
            "("
                "const dictionary&, "
                "CloudType&, "
                "const scalar&, "
                "const scalar&, "
                "const scalar&, "
                "const scalar&, "
                "const scalar&"
            ")"
        )   << "Unknown pendular wall model type type " << PendularWallModelType
            << ", constructor not in hash table" << nl << nl
            << "    Valid wall model types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << exit(FatalError);
    }

    return autoPtr<PendularWallModel<CloudType> >
    (
        cstrIter()
        (
            dict,
            owner,
            surfaceTension,
            contactAngle,
            liqFrac,
            viscosity,
            minSep
        )
    );
}


// ************************************************************************* //
