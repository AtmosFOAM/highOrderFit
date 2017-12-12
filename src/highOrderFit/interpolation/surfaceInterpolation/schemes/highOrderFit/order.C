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

#include "order.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highOrderFit::order::order
(
    const Foam::label x,
    const Foam::label y,
    const Foam::label z
)
:
    Foam::labelVector(x, y, z)
{}


Foam::highOrderFit::order::order(Istream& is)
{
    // Check state of Istream
    is.check("Foam::highOrderFit::order::order(Foam::Istream&)");

    operator>>(is, *this);
}


Foam::highOrderFit::order::order()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highOrderFit::order::~order()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::highOrderFit::order::factorialRatio
(
    const label dimensions
) const
{
    return
        scalar(factorial(x()) * factorial(y()) * factorial(z()))
        /
        scalar(factorial(x() + y() + z() + dimensions));
}


void Foam::highOrderFit::order::calculateExponentTensors
(
    Foam::List<Foam::labelTensor>& K
) const
{
    K.setSize(1);
    K[0] = labelTensor
    (
        labelVector(0, 0, 0),
        labelVector(0, 0, 0),
        labelVector(0, 0, 0)
    );
}


// ************************************************************************* //
