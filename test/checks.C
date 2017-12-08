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

#include "checks.H"

void Test::checkEqual
(
    const Foam::scalarList& actual,
    const Foam::scalarList& expected,
    Approx approx
)
{
    forAll(actual, i)
    {
        CHECK( actual[i] == approx(expected[i]) );
    }
}

void Test::checkEqual
(
    const Foam::scalarList& actual,
    const Foam::scalar expected,
    Approx approx
)
{
    forAll(actual, i)
    {
        CHECK( actual[i] == approx(expected) );
    }
}

void Test::checkEqual
(
    const Foam::vector actual,
    const Foam::vector expected,
    Approx approx
)
{
    CHECK( actual.x() == approx(expected.x()) );
    CHECK( actual.y() == approx(expected.y()) );
    CHECK( actual.z() == approx(expected.z()) );
}

void Test::checkEqual
(
    const Foam::scalarRectangularMatrix& actual,
    const Foam::scalarRectangularMatrix& expected,
    Approx approx
)
{
    REQUIRE( actual.m() == expected.m() );
    REQUIRE( actual.n() == expected.n() );

    for (Foam::label row = 0; row < actual.m(); row++)
    {
        for (Foam::label col = 0; col < actual.n(); col++)
        {
            CHECK( actual(row, col) == approx(expected(row, col)) );
        }
    }
}

// ************************************************************************* //
