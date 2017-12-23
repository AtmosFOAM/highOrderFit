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

#include "weightsFieldPair.H"

#include "stencilField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highOrderFit::weightsFieldPair::weightsFieldPair
(
    const Foam::fvMesh& mesh,
    const Foam::extendedUpwindCellToFaceStencil& stencils,
    const Foam::highOrderFit::transformer& transformer,
    const Foam::List<Foam::highOrderFit::order>& moments,
    const Foam::highOrderFit::multipliers& multipliers
)
:
Foam::MeshObject
<
    Foam::fvMesh,
    Foam::MoveableMeshObject,
    Foam::highOrderFit::weightsFieldPair
>
(
    mesh
),
owner_
(
    mesh,
    highOrderFit::stencilField
    (
        stencils.ownStencil(),
        stencils.ownMap(),
        mesh,
        transformer
    ),
    moments,
    multipliers
),
neighbour_
(
    mesh,
    highOrderFit::stencilField
    (
        stencils.neiStencil(),
        stencils.neiMap(),
        mesh,
        transformer
    ),
    moments,
    multipliers
)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highOrderFit::weightsFieldPair::~weightsFieldPair()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::highOrderFit::weightsFieldPair::movePoints()
{
    return true;
}

// ************************************************************************* //
