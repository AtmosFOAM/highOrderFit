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
{
    targetFace_.translate(*this);
    targetFace_.rotate(*this);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highOrderFit::stencil::~stencil()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::highOrderFit::stencil::targetFaceMoment
(
    const Foam::highOrderFit::order& o
) const
{
    return targetFace_.moment(o);
}

// ************************************************************************* //
