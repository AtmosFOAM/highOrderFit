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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highOrderFit::targetFace::targetFace
(
    const Foam::fvMesh& mesh,
    const Foam::label facei
)
:
Cf_(mesh.Cf()[facei]),
Sf_(mesh.Sf()[facei])
{}

Foam::highOrderFit::targetFace::targetFace
(
    const Foam::point Cf,
    const Foam::vector Sf,
    const Foam::highOrderFit::face f
)
:
Cf_(Cf),
Sf_(Sf)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highOrderFit::targetFace::~targetFace()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::highOrderFit::targetFace::translate
(
    Foam::highOrderFit::stencil& stencil
) const
{
    forAll(stencil, i)
    {
        stencil[i].translate(-Cf_);
    }
}


void Foam::highOrderFit::targetFace::rotate
(
    Foam::highOrderFit::stencil& stencil
) const
{
    forAll(stencil, i)
    {
        stencil[i].rotate(Sf_/mag(Sf_), vector(1, 0, 0));
    }
}


// ************************************************************************* //
