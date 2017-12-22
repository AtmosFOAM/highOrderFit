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
#include "tet.H"
#include "transform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highOrderFit::cell::cell
(
    const Foam::fvMesh& mesh,
    const Foam::label celli
)
:
List<face>(mesh.cells()[celli].size())
{
    const labelList& faces = mesh.cells()[celli];
    forAll(faces, i)
    {
        const label facei = faces[i];
        (*this)[i] = face(mesh, facei);

        if (facei < mesh.nInternalFaces() && mesh.neighbour()[facei] == celli)
        {
            (*this)[i].flip();
        }
    }
}


Foam::highOrderFit::cell::cell
(
    const Foam::List<Foam::highOrderFit::face>& faces
)
:
Foam::List<Foam::highOrderFit::face>(faces)
{}


Foam::highOrderFit::cell::cell(Istream& is)
{
    // Check state of Istream
    is.check("Foam::highOrderFit::cell::cell(Foam::Istream&)");

    operator>>(is, *this);
}


Foam::highOrderFit::cell::cell()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highOrderFit::cell::~cell()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::highOrderFit::cell::transform
(
    const Foam::vector x,
    const Foam::vector from,
    const Foam::vector to
)
{
    const quaternion q(rotationTensor(from, to));

    forAll((*this), facei)
    {
        (*this)[facei].transform(x, from, to);
    }
}


Foam::scalar Foam::highOrderFit::cell::moment
(
    const Foam::highOrderFit::order& o
) const
{
    scalar moment = 0.0;

    forAll((*this), facei)
    {
        List<tet> faceTets;
        (*this)[facei].decompose(faceTets);

        forAll(faceTets, teti)
        {
            moment += faceTets[teti].volumeMoment(o);
        }
    }

    if (o == order(0, 0, 0) && moment < SMALL)
    {
        FatalErrorInFunction
            << "Cell has zero volume "
            << (*this) << endl
            << exit(FatalError);
    }

    return moment;
}


// ************************************************************************* //
