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

#include "fvMesh.H"
#include "List.H"
#include "order.H"
#include "point.H"
#include "quaternion.H"
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

public:

    // Constructors

        //- Construct an instance for a given face.
        //  @param mesh a mesh instance
        //  @param facei an index into <tt>mesh.faces()</tt>
        face(const fvMesh& mesh, const label facei);

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
        
        //- Translate the face by vector \c x.
        void translate(const vector x);

        //- Rotate the face to align \c from with \c to.
        //  The transformation is constrained by the quaternion rotation, so
        //  faces cannot be rotated 180&deg; or more.
        void rotate(const vector from, const vector to);

        //- Rotate the face by the given quaternion.
        void rotate(const quaternion& q);

        //- Calculate the surface moment of the given order.
        scalar moment(const order&) const;

        //- Flip this face by reversing the order of its points.
        void flip();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
