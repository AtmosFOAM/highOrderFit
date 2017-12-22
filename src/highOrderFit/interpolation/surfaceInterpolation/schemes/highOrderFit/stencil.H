/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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
    Foam::highOrderFit::stencil

Description
    A stencil of \ref Foam::highOrderFit::cell "cells" used to approximate
    the high order surface interpolation.

SourceFiles
    stencil.C

\*---------------------------------------------------------------------------*/

#ifndef stencil_H
#define stencil_H

#include "cell.H"
#include "label.H"
#include "order.H"
#include "scalarList.H"
#include "targetFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace highOrderFit
{

/*---------------------------------------------------------------------------*\
                         Class stencil Declaration
\*---------------------------------------------------------------------------*/

class stencil
:
    public List<cell>
{

    // Private data

        targetFace targetFace_;
    

    // Private Member Functions

        //- Disallow default bitwise copy construct
        stencil(const stencil&);


public:

    // Constructors

        //- Construct an empty stencil with no cells.
        //  Necessary only for initialising collections of stencils.
        stencil();

        //- Construct a stencil comprising the given cells, transformed into the local coordinate system.
        //  The target face centre is used as the origin of the local
        //  coordinate system, and the stencil is rotated such that the
        //  primary direction of the local coordinate system is parallel to the
        //  target face area vector.
        //  @param targetFace the stencil's target face that determines the
        //                    origin and orientation of the local coordinate
        //                    system
        stencil
        (
            const targetFace& targetFace,
            const List<cell>& cells
        );


    //- Destructor
    ~stencil();


    // Member Functions

        void transform(const vector x, const vector from, const vector to);

        const targetFace& target() const
        {
            return targetFace_;
        }
    
        point origin() const
        {
            return targetFace_.Cf();
        }

        vector primaryDirection() const
        {
            return targetFace_.unitNormal();
        }

        //- Calculate a moment for the stencil's target face, divided by the face area.
        //  \sa Foam::highOrderFit::targetFace::momentAverage
        scalar targetFaceMomentAverage(const order&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
