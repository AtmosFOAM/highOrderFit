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

#include "sphericalTransformer.H"
#include "addToRunTimeSelectionTable.H"
#include "axesRotation.H"
#include "cartesianCS.H"
#include "transform.H"

namespace Foam
{
namespace highOrderFit
{

defineTypeNameAndDebug(sphericalTransformer, 0);
addToRunTimeSelectionTable
(
    transformer,
    sphericalTransformer,
    word
);

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highOrderFit::sphericalTransformer::sphericalTransformer()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highOrderFit::sphericalTransformer::~sphericalTransformer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::highOrderFit::sphericalTransformer::transform
(
    Foam::highOrderFit::stencil& stencil
) const
{
    const vector i = stencil.primaryDirection()/mag(stencil.primaryDirection());
    const vector k = stencil.origin()/mag(stencil.origin());
    const vector j = k ^ i;

    const tensor R(i, j, k);
    const axesRotation* rotation = new axesRotation(R.T());

    const cartesianCS coordinates
    (
        "stencilCoordinates",
        stencil.origin(),
        *rotation
    );
    stencil.transform(coordinates);
}

// ************************************************************************* //
