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
#include "cell.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "labelList.H"
#include "mesh.H"
#include "OStringStream.H"
#include "testCase.H"

using namespace Foam;

namespace Test
{

TEST_CASE("cell_has_six_faces_for_cuboid_cell")
{
    const Test::testCase c("cartesian4x3Mesh");
    const Test::mesh testMesh(c.mesh());

    const label celli = testMesh.indexOfCellWithCentreAt(point(0.5, 0.5, 0));
    const highOrderFit::cell cell(c.mesh(), celli);

    CHECK( cell.size() == 6 );
}

TEST_CASE("cell_has_four_cell_for_each_face_of_cuboid_cell")
{
    const Test::testCase c("cartesian4x3Mesh");
    const Test::mesh testMesh(c.mesh());

    const label celli = testMesh.indexOfCellWithCentreAt(point(0.5, 0.5, 0));
    const highOrderFit::cell cell(c.mesh(), celli);

    forAll(cell, i)
    {
        CHECK( cell[i].size() == 4 );
    }
}

TEST_CASE("cell_has_cell_of_face")
{
    const Test::testCase c("cartesian4x3Mesh");
    const Test::mesh testMesh(c.mesh());

    const label celli = testMesh.indexOfCellWithCentreAt(point(0.5, 0.5, 0));
    const highOrderFit::cell cell(c.mesh(), celli);

    const List<point>& firstFace = cell[0];
    Test::checkEqual(firstFace[0], point(1, 0, -0.5));
    Test::checkEqual(firstFace[1], point(1, 1, -0.5));
    Test::checkEqual(firstFace[2], point(1, 1,  0.5));
    Test::checkEqual(firstFace[3], point(1, 0,  0.5));
}

TEST_CASE("cell_round_trips_through_IOstreams")
{
    List<List<point>> points(1);
    points[0].setSize(1);
    points[0][0] = point(2, 3, 4);

    const highOrderFit::cell cellOut(points);

    OStringStream o;
    o << cellOut << endl;

    highOrderFit::cell cellIn;

    IStringStream i(o.str());
    i >> cellIn;

    REQUIRE( cellIn.size() == 1 );
    REQUIRE( cellIn[0].size() == 1 );
    checkEqual(cellIn[0][0], point(2, 3, 4));
}

TEST_CASE("cell_calculates_zeroth_volume_moment_for_unit_cube")
{
    highOrderFit::cell cell;
    IFstream i("resources/unitCube");
    i >> cell;

    CHECK( cell.moment(highOrderFit::order(0, 0, 0)) == approx(1.0) );
}

TEST_CASE("cell_calculates_x_volume_moment_for_unit_cube_centred_at_origin",
          "[!mayfail]")
{
    highOrderFit::cell cell;
    IFstream i("resources/unitCube");
    i >> cell;

    CHECK( cell.moment(highOrderFit::order(1, 0, 0)) == approx(0.0) );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

// ************************************************************************* //
