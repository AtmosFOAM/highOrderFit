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

#include "stencil.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highOrderFit::stencil::stencil()
{
}

Foam::highOrderFit::stencil::stencil
(
    const Foam::highOrderFit::targetFace& targetFace,
    const Foam::List<Foam::highOrderFit::cell>& cells
)
:
Foam::List<Foam::highOrderFit::cell>(cells),
targetFace_(targetFace)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highOrderFit::stencil::~stencil()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::highOrderFit::stencil::transform
(
    const Foam::vector x,
    const Foam::vector from,
    const Foam::vector to
)
{
    targetFace_.transform(x, from, to);

    forAll((*this), celli)
    {
        (*this)[celli].transform(x, from, to);
    }
}

Foam::scalar Foam::highOrderFit::stencil::targetFaceMomentAverage
(
    const Foam::highOrderFit::order& o
) const
{
    return targetFace_.momentAverage(o);
}

// ************************************************************************* //
