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
    Foam::highOrderFit::tet

Description
    A tetrahedron with one point at the coordinate origin.

SourceFiles
    tet.C

\*---------------------------------------------------------------------------*/

#ifndef tet_H
#define tet_H

#include "Vector.H"
#include "order.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace highOrderFit
{

/*---------------------------------------------------------------------------*\
                         Class tet Declaration
\*---------------------------------------------------------------------------*/

class tet
    : public Vector<point>
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        tet(const tet&);

        // Disallow default bitwise assignment
        //void operator=(const tet&);
        
        scalar moment
        (
            const order& order,
            const tensor& A,
            const label dimensions
        ) const;

        scalar area() const;

public:

    // Constructors

        //- Construct an uninitialised instance.
        //  Necessary only for initialising collections of tets.
        tet();

        //- Construct from three points, with the fourth point fixed at the coordinate origin.
        //  The point order establishes a right-hand orientation for the face
        //  <tt>(a, b, c)</tt>.
        tet(const point a, const point b, const point c);


    //- Destructor
    ~tet();


    // Member Functions

        //- Calculate the volume moment of the given order.
        //  Volumes may be negative depending on the orientation of the face
        //  <tt>(a, b, c)</tt>.
        //  \sa Foam::highOrderFit::cell::moment
        scalar volumeMoment(const order& o) const;

        //- Calculate the surface moment of the given order.
        //  \sa Foam::highOrderFit::face::moment
        scalar surfaceMoment(const order& o) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
