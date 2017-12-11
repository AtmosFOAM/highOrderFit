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
#include "labelList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highOrderFit::order::order
(
    const Foam::label x,
    const Foam::label y,
    const Foam::label z
)
:
x_(x),
y_(y),
z_(z)
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


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::highOrderFit::operator>>
(
    Foam::Istream& is,
    Foam::highOrderFit::order& o
)
{
    labelList l;
    is >> l;

    if (l.size() != 3)
    {
        FatalErrorInFunction
            << "order must have three components but found " << l.size()
            << abort(FatalError);
    }

    o.x_ = l[0];
    o.y_ = l[1];
    o.z_ = l[2];

    // Check state of Istream
    is.check
    (
        "Foam::Ostream& Foam::highOrderFit::operator>>(Foam::Istream&, "
        "Foam::highOrderFit::order&)"
    );

    return is;
}

Foam::Ostream& Foam::highOrderFit::operator<<
(
    Foam::Ostream& os,
    const Foam::highOrderFit::order& o
)
{
    const labelList l({o.x_, o.y_, o.z_});
    os << l;

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::highOrderFit::operator<<(Foam::Ostream&, "
        "const Foam::highOrderFit::order&)"
    );

    return os;
}


// ************************************************************************* //
