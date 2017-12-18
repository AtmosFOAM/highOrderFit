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
#include "HashTable.H"
#include "tet.H"

using namespace Foam;

namespace Test
{

TEST_CASE("tet_calculates_zeroth_volume_moment")
{
    const point a(1, 0 ,0), b(0, 1, 0), c(0, 0, 1);
    highOrderFit::tet t(a, b, c);

    CHECK( t.volumeMoment(highOrderFit::order(0, 0, 0)) == approx(1.0/6.0) );
}

TEST_CASE("tet_calculates_1_0_0_volume_moment")
{
    const point a(2, 3, 1), b(1, 1, 4), c(4, 4, 2);
    highOrderFit::tet t(a, b, c);

    tensor A(a, b, c);
    A = A.T();

    CHECK( t.volumeMoment(highOrderFit::order(1, 0, 0)) ==
            approx(det(A)/24.0 * (a.x() + b.x() + c.x())) );
}

TEST_CASE("tet_calculates_0_2_0_volume_moment")
{
    const point a(2, 3, 1), b(1, 1, 4), c(4, 4, 2);
    highOrderFit::tet t(a, b, c);

    tensor A(a, b, c);
    A = A.T();

    CHECK
    (
        t.volumeMoment(highOrderFit::order(0, 2, 0))
        ==
        approx
        (
            det(A)/60.0 * 
            (
                sqr(a.y()) + sqr(b.y()) + sqr(c.y())
              + a.y()*b.y() + a.y()*c.y() + b.y()*c.y()
            )
        )
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

// ************************************************************************* //
