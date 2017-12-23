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
#include "axesRotation.H"
#include "cartesianCS.H"
#include "checks.H"
#include "face.H"

using namespace Foam;

namespace Test
{

TEST_CASE("face_triangle_face_decomposes_into_three_tets_with_common_centre")
{
    const highOrderFit::face face
    (
        {
            point(0, 0, 0),
            point(1, 0, 0),
            point(0, 1, 0)
        }
    );

    List<highOrderFit::tet> tets;
    face.decompose(tets);

    const point expectedCentre(1.0/3.0, 1.0/3.0, 0);

    CHECK( tets.size() == 3 );
    forAll(tets, teti)
    {
        CHECK( countMatches(tets[teti], expectedCentre) == 1 );
    }
}

TEST_CASE("face_triangle_face_decomposes_into_tets_having_original_orientation")
{
    const highOrderFit::face face
    (
        {
            point(0, 0, 0),
            point(1, 0, 0),
            point(0, 1, 0)
        }
    );
      
    List<highOrderFit::tet> tets;
    face.decompose(tets);

    const vector expectedOrientation(0, 0, 1);

    forAll(tets, teti)
    {
        const highOrderFit::tet& t = tets[teti];
        REQUIRE( t.size() == 3 );
           
        vector actualOrientation = (t[1] - t[0]) ^ (t[2] - t[1]);
        actualOrientation /= mag(actualOrientation);

        checkEqual(actualOrientation, expectedOrientation);
    }
}

TEST_CASE("face_calculates_1_0_0_moment_for_unit_square_centred_at_0.5_0_0")
{
    highOrderFit::face face
    (
        {
            point(-0.5, -0.5, 0),
            point( 0.5, -0.5, 0),
            point( 0.5,  0.5, 0),
            point(-0.5,  0.5, 0)
        }
    );
    cartesianCS coordinates
    (
        "faceTransform",
        point(-0.5, 0, 0),
        *(new axesRotation(I))
    );
    face.transform(coordinates);

    CHECK( face.moment(highOrderFit::order(1, 0, 0)) == approx(0.5) );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

// ************************************************************************* //
