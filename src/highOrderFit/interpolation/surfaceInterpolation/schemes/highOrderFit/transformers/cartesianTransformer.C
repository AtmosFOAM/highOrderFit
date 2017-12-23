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

#include "cartesianTransformer.H"
#include "addToRunTimeSelectionTable.H"
#include "axesRotation.H"
#include "cartesianCS.H"
#include "transform.H"

namespace Foam
{
namespace highOrderFit
{

defineTypeNameAndDebug(cartesianTransformer, 0);
addToRunTimeSelectionTable
(
    transformer,
    cartesianTransformer,
    word
);

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highOrderFit::cartesianTransformer::cartesianTransformer()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highOrderFit::cartesianTransformer::~cartesianTransformer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::highOrderFit::cartesianTransformer::transform
(
    Foam::highOrderFit::stencil& stencil
) const
{
    vector from = stencil.primaryDirection();
    const vector to(1, 0, 0);

    const scalar s = from & to;
    const scalar magSqrN3 = magSqr(from ^ to);

    if (magSqrN3 <= SMALL && s < 0)
    {
        from = vector(1, 0, 0);
    }

    const axesRotation* rotation = new axesRotation(rotationTensor(to, from));
    const cartesianCS coordinates
    (
        "stencilCoordinates",
        stencil.origin(),
        *rotation
    );
    stencil.transform(coordinates);
}

// ************************************************************************* //
