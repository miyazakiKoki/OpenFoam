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

Class
    Foam::VoidFractionDivided

Description
    Creates particle void fraction field on carrier phase

SourceFiles
    VoidFractionDivided.C

\*---------------------------------------------------------------------------*/

#ifndef VoidFractionDivided_H
#define VoidFractionDivided_H

#include "CloudFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class VoidFractionDivided Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class VoidFractionDivided
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;


        //- Void fraction field
        autoPtr<volScalarField> thetaPtr_;


        //- Refinement level means how finely satellite particles are distributed
        label refinementLevel_;

        //- Number of the satellite particles
        label satParNo_;

        //- Position of the satellite particles
        List<vector> satPositions_;

        //- TreeBoundBox for sub-domain allocated to processor
        treeBoundBox procBb_;

        //- Maximum interaction distance
        scalar maxDistance_;

        //- TreeBoundBox for extended sub-domain allocated to processor
        treeBoundBox extendedProcBb_;

        //- List of all extendedProcBbs
        treeBoundBoxList allExtendedProcBbs_;

        //- List of treeBoundBoxes for all the cells
        treeBoundBoxList cellBbs_;

        //- List of extended treeBoundBoxes for all the cells
        treeBoundBoxList extendedCellBbs_;

        //- List of cells where a particle inside the cell may interact
        labelListList allInteractionList_;

        labelList neighbouringProcs_;

        List<DynamicList<vector> > positionTransferList_;

        List<DynamicList<labelList> > cellLabelTransferList_;

        //- List of neighbouring cell labels in other processors
        labelListList interactionCellLabelList_;

       //- List of neighouring cell processors in other processors
       labelListList interactionCellProcList_;

protected:

    // Protected Member Functions

        //- Write post-processing info
        virtual void write();


public:

    //- Runtime type information
    TypeName("voidFractionDivided");


    // Constructors

        //- Construct from dictionary
        VoidFractionDivided(const dictionary& dict, CloudType& owner);

        //- Construct copy
        VoidFractionDivided(const VoidFractionDivided<CloudType>& vf);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType> > clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType> >
            (
                new VoidFractionDivided<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~VoidFractionDivided();


    // Member Functions

        // Evaluation

            //- Pre-evolve hook
            virtual void preEvolve();

            //- Post-evolve hook
            virtual void postEvolve();

            //- Post-move hook
            virtual void postMove
            (
                const parcelType& p,
                const label cellI,
                const scalar dt,
                const point& position0,
                bool& keepParticle
            );

            //- Build cellBbs for interaction
            void buildInteractionList();

            //- Check if point is inside cell
            bool isInside(const point&, label) const;


        // Access

            // Return const access to the refinement level
            label refinementLevel() const;

            // Return const access to the number of the satellite particles
            label satParNo() const;

            // Return const access to the position of the satellite particles
            List<vector> satPositions() const;

            //Return const access to the maximum interaction distance
            scalar maxDistance() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VoidFractionDividedI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "VoidFractionDivided.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
