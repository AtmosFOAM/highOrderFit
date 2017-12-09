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

#include "cell.H"
#include "face.H"
#include "labelList.H"
#include "quaternion.H"
#include "transform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highOrderFit::cell::cell
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


Foam::highOrderFit::cell::cell
(
    const Foam::List<Foam::List<Foam::point>>& points
)
:
List<List<point>>(points)
{}


Foam::highOrderFit::cell::cell()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highOrderFit::cell::~cell()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::highOrderFit::cell::translate(const Foam::vector x)
{
    forAll((*this), facei)
    {
        forAll((*this)[facei], pointi)
        {
            (*this)[facei][pointi] += x;
        }
    }
}


void Foam::highOrderFit::cell::rotate
(
    const Foam::vector from,
    const Foam::vector to
)
{
    const quaternion q(rotationTensor(from, to));

    forAll((*this), facei)
    {
        forAll((*this)[facei], pointi)
        {
            (*this)[facei][pointi] = q.transform((*this)[facei][pointi]);
        }
    }
}


Foam::scalar Foam::highOrderFit::cell::moment
(
    const Foam::highOrderFit::order& o
) const
{
    if (o == order(0, 0, 0))
    {
        return 1.0;
    }
    else if (o == order(1, 0, 0))
    {
        return average().x();
    }
    else
    {
        const scalar x = average().x();
        if (x <= -2.5 + SMALL)
        {
            return 19.0/3.0;
        }
        else if (x <= -1.5 + SMALL)
        {
            return 7.0/3.0;
        }
        else
        {
            return 1.0/3.0;
        }
    }
}


Foam::point Foam::highOrderFit::cell::average() const
{
    point average = point(0, 0, 0);
    label points = 0;

    forAll((*this), facei)
    {
        forAll((*this)[facei], pointi)
        {
            average += (*this)[facei][pointi];
            points++;
        }
    }
    average /= points;

    return average;
}

// ************************************************************************* //
