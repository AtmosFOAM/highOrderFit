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

#include "catch.hpp"
#include "checks.H"
#include "exponentTensor.H"

using namespace Foam;

namespace Test
{

TEST_CASE("exponentTensor_calculates_product_of_exponentials")
{
    const tensor A
    (
        10, 11, 12,
        13, 14, 15,
        16, 17, 18
    );

    const highOrderFit::exponentTensor K
    (
        1, 2, 3,
        4, 5, 6,
        7, 8, 9
    );

    CHECK( K.productOfExponentials(A) ==
            approx(13588579074408578757989397706630967514844102656e7) );
}

TEST_CASE("exponentTensor_calculates_factorial_ratio")
{
    const highOrderFit::exponentTensor K
    (
        0, 1, 0,
        8, 7, 6,
        3, 4, 5
    );

    CHECK( K.factorialRatio() == approx(3.018708e8) );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test


// ************************************************************************* //
