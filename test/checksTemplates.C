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
#include "label.H"

template<class Type>
void Test::checkEqual
(
    const Foam::List<Type>& actual,
    const Type expected
)
{
    forAll(actual, i)
    {
        CHECK( actual[i] == expected );
    }
}

template<class Type>
void Test::checkEqual
(
    const Foam::List<Type>& actual,
    const Foam::List<Type>& expected
)
{
    forAll(actual, i)
    {
        CHECK( actual[i] == expected[i] );
    }
}

template<class Type>
void Test::checkEqual
(
    const Foam::Tensor<Type>& actual,
    const Foam::Tensor<Type>& expected
)
{
    CHECK( actual.xx() == expected.xx() );
}

template<class Type>
Foam::label Test::countMatches
(
    const Foam::List<Type>& list,
    const Type item
)
{
    Foam::label matches = 0;

    forAll(list, i)
    {
        if (list[i] == item) matches++;
    }
     
    return matches;
}

template<class Type>
Foam::label Test::countMatches
(
    const Foam::Vector<Type>& vector,
    const Type item
)
{
    Foam::label matches = 0;

    if (vector.x() == item) matches++;
    if (vector.y() == item) matches++;
    if (vector.z() == item) matches++;
     
    return matches;
}

// ************************************************************************* //
