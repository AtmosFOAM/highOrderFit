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
#include "face.H"
#include "faceCache.H"
#include "tet.H"

using namespace Foam;

namespace Test
{

TEST_CASE("faceCache_round_trips_new_face")
{
    const highOrderFit::face face({point(0, 0, 0)});
    const highOrderFit::tet tet
    (
        point(1, 2, 3),
        point(4, 5, 6),
        point(7, 8, 9)
    );

    List<highOrderFit::tet> tetsIn;
    tetsIn.append(tet);

    highOrderFit::faceCache cache;

    const highOrderFit::faceCache::decomposer func =
        [&tet](List<highOrderFit::tet>& target)
    {
        target.append(tet);
    };

    const List<highOrderFit::tet>& tetsOut =
        cache.calculateIfNecessary(face, func);

    REQUIRE( tetsOut.size() == 1 );
    checkEqual(tetsOut[0].x(), point(1, 2, 3));
    checkEqual(tetsOut[0].y(), point(4, 5, 6));
    checkEqual(tetsOut[0].z(), point(7, 8, 9));
}

TEST_CASE("faceCache_caches_identical_face")
{
    const highOrderFit::face face({point(0, 0, 0)});
    const highOrderFit::tet tet
    (
        point(1, 2, 3),
        point(4, 5, 6),
        point(7, 8, 9)
    );

    List<highOrderFit::tet> tetsIn;
    tetsIn.append(tet);

    highOrderFit::faceCache cache;
    label funcCalls = 0;

    const highOrderFit::faceCache::decomposer func =
        [&tet, &funcCalls](List<highOrderFit::tet>& target)
    {
        target.append(tet);
        funcCalls++;
    };

    cache.calculateIfNecessary(face, func);
    cache.calculateIfNecessary(face, func);

    CHECK( funcCalls == 1 );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

// ************************************************************************* //
