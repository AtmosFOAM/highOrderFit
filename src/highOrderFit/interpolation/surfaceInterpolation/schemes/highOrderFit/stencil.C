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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::highOrderFit::stencil::translateVerticesWithOrigin
(
    const Foam::point o
)
{
    forAll((*this), i)
    {
        (*this)[i].translate(-o);
    }
}


void Foam::highOrderFit::stencil::rotateVerticesWithPrimaryDirection
(
    const Foam::vector x
)
{
    forAll((*this), i)
    {
        (*this)[i].rotate(x/mag(x), vector(1, 0, 0));
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highOrderFit::stencil::stencil()
:
size_(0)
{}


Foam::highOrderFit::stencil::stencil
(
    const Foam::point targetCf,
    const Foam::vector Sf,
    const Foam::List<Foam::highOrderFit::cellVertices>& cellVertices
)
:
Foam::List<Foam::highOrderFit::cellVertices>(cellVertices),
size_(cellVertices.size()),
zeroMoment_(cellVertices.size(), 1.0)
{
    translateVerticesWithOrigin(targetCf);
    rotateVerticesWithPrimaryDirection(Sf);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highOrderFit::stencil::~stencil()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::highOrderFit::stencil::size() const
{
    return size_;
}

const Foam::scalarList& Foam::highOrderFit::stencil::moment
(
    const order& order
) const
{
    return zeroMoment_;
}

// ************************************************************************* //
