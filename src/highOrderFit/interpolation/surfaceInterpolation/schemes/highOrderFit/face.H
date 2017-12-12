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
    Foam::highOrderFit::face

Description
    A right-hand oriented list of face vertex points.

SourceFiles
    face.C

\*---------------------------------------------------------------------------*/

#ifndef highOrderFit_face_H
#define highOrderFit_face_H

#include "List.H"
#include "point.H"
#include "tet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace highOrderFit
{

/*---------------------------------------------------------------------------*\
                         Class face Declaration
\*---------------------------------------------------------------------------*/

class face
:
    public List<point>
{

    // Private Member Functions

        // Disallow default bitwise copy construct
        //face(const face&);

        // Disallow default bitwise assignment
        //void operator=(const face&);


public:

    // Constructors

        //- Construct an empty instance with no vertex points.
        //  Necessary only for initialising collections of faces.
        face();

        //- Construct from the given initialiser list.
        face(std::initializer_list<point> lst);

    //- Destructor
    ~face();


    // Member Functions

        //- Decompose into a collection of tetrahedra.
        //  Each tetrahedron has a point at the origin and the opposite face
        //  having the same orientation as the original polygonal face.
        //  @param [out] tets on entry, an empty list; on exit, the tets will
        //               be initialised
        void decompose(List<tet>& tets) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
