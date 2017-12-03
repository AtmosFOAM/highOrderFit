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

\*---------------------------------------------------------------------------*/

#include "cellVertices.H"
#include "face.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highOrderFit::cellVertices::cellVertices()
{}

Foam::highOrderFit::cellVertices::cellVertices
(
    const Foam::primitiveMesh& mesh,
    const Foam::label celli
)
:
List<List<point>>(mesh.cells()[celli].size())
{
    const labelList& faces = mesh.cells()[celli];
    forAll(faces, i)
    {
        List<point>& faceVertices = (*this)[i];
        const face& face = mesh.faces()[faces[i]];

        faceVertices.setSize(face.size());
        forAll(face, pointi)
        {
            faceVertices[pointi] = mesh.points()[face[pointi]];
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highOrderFit::cellVertices::~cellVertices()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
