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
    // Private Member Functions

        //- Disallow default bitwise copy construct
        stencil(const stencil&);

        //Disallow default bitwise assignment
        //void operator=(const stencil&);

        void translateVerticesWithOrigin(const point o);

        void rotateVerticesWithPrimaryDirection(const vector x);


public:

    // Constructors

        //- Construct an empty stencil with no cells.
        //  Necessary only for initialising collections of stencils.
        stencil();

        //- Construct a stencil comprising the given cells, transformed into the local coordinate system.
        //  @param targetCf the stencil's target face centre to be used as
        //                  the origin of the local coordinate system
        //  @param Sf       the target face area vector.  The stencil is rotated
        //                  such that primary direction of the local coordinate
        //                  system is parallel to <tt>Sf</tt>.
        stencil
        (
            const point targetCf,
            const vector Sf,
            const List<cell>& cells
        );


    //- Destructor
    ~stencil();


    // Member Functions

        scalar targetFaceMoment(const order&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
