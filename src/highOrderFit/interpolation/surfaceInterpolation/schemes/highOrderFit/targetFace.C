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

#include "targetFace.H"
#include "surfaceMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::highOrderFit::targetFace::transform()
{
    face::translate(-Cf_);
    face::rotate(unitNormal_, vector(1, 0, 0));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highOrderFit::targetFace::targetFace()
:
Cf_(0, 0, 0),
unitNormal_(0, 0, 0),
area_(0)
{}

Foam::highOrderFit::targetFace::targetFace
(
    const Foam::fvMesh& mesh,
    const Foam::label facei
)
:
face(mesh, facei),
Cf_(mesh.Cf()[facei]),
unitNormal_(mesh.Sf()[facei]/mag(mesh.Sf()[facei])),
area_(mag(mesh.Sf()[facei]))
{
    transform();
}

Foam::highOrderFit::targetFace::targetFace
(
    std::initializer_list<point> lst,
    const Foam::point Cf,
    const Foam::vector Sf
)
:
face(lst),
Cf_(Cf),
unitNormal_(Sf/mag(Sf)),
area_(mag(Sf))
{
    transform();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highOrderFit::targetFace::~targetFace()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::highOrderFit::targetFace::transform
(
    Foam::List<Foam::highOrderFit::cell>& stencil
) const
{
    forAll(stencil, i)
    {
        stencil[i].transform(-Cf_, unitNormal_, vector(1, 0, 0));
    }
}


Foam::scalar Foam::highOrderFit::targetFace::momentAverage
(
    const Foam::highOrderFit::order& o
) const
{
    return face::moment(o) / area_;
}

// ************************************************************************* //
