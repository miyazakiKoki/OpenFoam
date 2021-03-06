/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "WetParcel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::WetParcel<ParcelType>::WetParcel
(
    const WetParcel<ParcelType>& p
)
:
    ParcelType(p),
    Vliq_(p.Vliq_),
    partVliq_(p.partVliq_),
    liquidPositionVectors_(p.liquidPositionVectors_),
    liquidPositions_(p.liquidPositions_),
    contactList_(p.contactList_),
    previousContactList_(p.previousContactList_)
{}


template<class ParcelType>
Foam::WetParcel<ParcelType>::WetParcel
(
    const WetParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    Vliq_(p.Vliq_),
    partVliq_(p.partVliq_),
    liquidPositionVectors_(p.liquidPositionVectors_),
    liquidPositions_(p.liquidPositions_),
    contactList_(p.contactList_),
    previousContactList_(p.previousContactList_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
bool Foam::WetParcel<ParcelType>::move
(
    TrackData& td,
    const scalar trackTime
)
{
    typename TrackData::cloudType::parcelType& p =
        static_cast<typename TrackData::cloudType::parcelType&>(*this);

    td.switchProcessor = false;
    td.keepParticle = true;

    const polyMesh& mesh = td.cloud().pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();
    const scalarField& V = mesh.cellVolumes();

    switch (td.part())
    {
        case TrackData::tpVelocityHalfStep:
        {
            // First and last leapfrog velocity adjust part, required
            // before and after tracking and force calculation

            p.U() += 0.5*trackTime*p.f()/p.mass();

            p.angularMomentum() += 0.5*trackTime*p.torque();

            break;
        }

        case TrackData::tpLinearTrack:
        {
            scalar tEnd = (1.0 - p.stepFraction())*trackTime;
            const scalar dtMax = tEnd;

            while (td.keepParticle && !td.switchProcessor && tEnd > ROOTVSMALL)
            {
                // Apply correction to position for reduced-D cases
                meshTools::constrainToMeshCentre(mesh, p.position());

                const point start(p.position());

                // Set the Lagrangian time-step
                scalar dt = min(dtMax, tEnd);

                // Remember which cell the parcel is in since this
                // will change if a face is hit
                const label cellI = p.cell();

                const scalar magU = mag(p.U());
                if (p.active() && magU > ROOTVSMALL)
                {
                    const scalar d = dt*magU;
                    const scalar maxCo = td.cloud().solution().maxCo();
                    const scalar dCorr = min(d, maxCo*cbrt(V[cellI]));
                    dt *=
                        dCorr/d
                       *p.trackToFace(p.position() + dCorr*p.U()/magU, td);
                }

                tEnd -= dt;
                p.stepFraction() = 1.0 - tEnd/trackTime;

                // Avoid problems with extremely small timesteps
                if (dt > ROOTVSMALL)
                {
                    // Update cell based properties
                    p.setCellValues(td, dt, cellI);

                    if (td.cloud().solution().cellValueSourceCorrection())
                    {
                        p.cellValueSourceCorrection(td, dt, cellI);
                    }

                    p.calc(td, dt, cellI);

/*-----------------------*/
                    p.liquidPositions().setSize(p.liquidPositionVectors().size()); 
                    List<vector> liquidPositionVelocitys(p.liquidPositionVectors().size());
                    forAll(liquidPositionVectors_,i)
                    {

                        if(p.omega()==vector::zero)
                        {
                            p.liquidPositions()[i]=p.liquidPositionVectors()[i]*p.d()/2+p.position();
                            p.liquidPositionVectors()[i]=p.liquidPositionVectors()[i]/mag(p.liquidPositionVectors()[i]);
                        }else
                        {
 /*          liquidPositionVelocitys[i]=(p.d()/2)*(p.omega()^p.liquidPositionVectors()[i]);
           p.liquidPositions()[i] = liquidPositionVelocitys[i]*dt+p.liquidPositionVectors()[i]*p.d()/2+p.position();
           p.liquidPositionVectors()[i] = p.liquidPositions()[i]-p.position();
           p.liquidPositionVectors()[i]=p.liquidPositionVectors()[i]/mag(p.liquidPositionVectors()[i]);
*/
                            vector angle = dt*p.omega();
                            p.liquidPositionVectors()[i]=p.liquidPositionVectors()[i]/mag(p.liquidPositionVectors()[i]);
                            vector liquidPositionVectorsDummy = p.liquidPositionVectors()[i];     

                            p.liquidPositionVectors()[i].x() = ( cos(angle.x())*cos(angle.y())*cos(angle.z())-sin(angle.x())*sin(angle.z()) )*liquidPositionVectorsDummy.x()*p.d()/2 + ( -cos(angle.x())*cos(angle.y())*sin(angle.z())-sin(angle.x())*cos(angle.z()) )*liquidPositionVectorsDummy.y()*p.d()/2 + cos(angle.x())*sin(angle.y())*liquidPositionVectorsDummy.z()*p.d()/2;
           
                            p.liquidPositionVectors()[i].y() = ( sin(angle.x())*cos(angle.y())*cos(angle.z())+cos(angle.x())*sin(angle.z()) )*liquidPositionVectorsDummy.x()*p.d()/2 + ( -sin(angle.x())*cos(angle.y())*sin(angle.z())+cos(angle.x())*cos(angle.z()) )*liquidPositionVectorsDummy.y()*p.d()/2 + sin(angle.x())*sin(angle.y())*liquidPositionVectorsDummy.z()*p.d()/2;

                            p.liquidPositionVectors()[i].z() = -sin(angle.y())*cos(angle.z())*liquidPositionVectorsDummy.x()*p.d()/2 + sin(angle.y())*sin(angle.z())*liquidPositionVectorsDummy.y()*p.d()/2 + cos(angle.y())*liquidPositionVectorsDummy.z()*p.d()/2;

                            p.liquidPositions()[i] = p.liquidPositionVectors()[i] + p.position();
           
                            p.liquidPositionVectors()[i]=p.liquidPositionVectors()[i]/mag(p.liquidPositionVectors()[i]);
            
                       }
                   }

                   //copying contactList to previousContactList
                   p.previousContactList().clear();
                   //Info<<"size of contactList of particle p before copying to previousContactList of particle" << p.origId() << "is"<< p.contactList().size()<<endl;
                   p.previousContactList().setSize(p.contactList().size());
                   p.previousContactList() = p.contactList();
                   //forAll(p.previousContactList_,i)
                   //{
                   //	Info<<"previousContactList of this particle is "<< p.previousContactList()<<endl;
                   // }
                   p.contactList().clear();
                   //Info << "contactList of particle " << p.origId() << "is cleared"<<endl;

/*-----------------------*/
                }

                if (p.onBoundary() && td.keepParticle)
                {
                    if (isA<processorPolyPatch>(pbMesh[p.patch(p.face())]))
                    {
                        td.switchProcessor = true;
                    }
                }

                p.age() += dt;

                td.cloud().functions().postMove
                (
                    p,
                    cellI,
                    dt,
                    start,
                    td.keepParticle
                );
            }

            break;
        }

        case TrackData::tpRotationalTrack:
        {
            notImplemented("TrackData::tpRotationalTrack");

            break;
        }

        default:
        {
            FatalErrorIn
            (
                "WetParcel<ParcelType>::move(TrackData&, const scalar)"
            )   << td.part() << " is an invalid part of the tracking method."
                << abort(FatalError);
        }
    }

    p.Vliq()=0;
    forAll(partVliq_,i)
    {
        p.Vliq()+=p.partVliq()[i];
    }





    return td.keepParticle;
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "WetParcelIO.C"

// ************************************************************************* //
