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
#include "face.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "labelList.H"
#include "mesh.H"
#include "OStringStream.H"
#include "testCase.H"

using namespace Foam;

namespace Test
{

void translate(highOrderFit::cell& c, const vector x)
{
    c.transform(x, vector(1, 0, 0), vector(1, 0, 0));
}

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

TEST_CASE("cell_has_face_points")
{
    const Test::testCase c("cartesian4x3Mesh");
    const Test::mesh testMesh(c.mesh());

    const label celli = testMesh.indexOfCellWithCentreAt(point(0.5, 0.5, 0));
    const highOrderFit::cell cell(c.mesh(), celli);

    const List<point>& firstFace = cell[0];
    checkEqual(firstFace[0], point(1, 0, -0.5));
    checkEqual(firstFace[1], point(1, 1, -0.5));
    checkEqual(firstFace[2], point(1, 1,  0.5));
    checkEqual(firstFace[3], point(1, 0,  0.5));
}

TEST_CASE("cell_has_centre_from_mesh")
{
    const Test::testCase c("cartesian4x3Mesh");
    const Test::mesh testMesh(c.mesh());

    const label celli = testMesh.indexOfCellWithCentreAt(point(0.5, 0.5, 0));
    const highOrderFit::cell cell(c.mesh(), celli);

    checkEqual(cell.centre(), point(0.5, 0.5, 0)); 
}

TEST_CASE("cell_has_centre_from_explicit_constructor")
{
    const point expectedCentre(2, 3, 4);
    const List<highOrderFit::face> faces;
    const highOrderFit::cell cell(faces, expectedCentre);

    checkEqual(cell.centre(), expectedCentre);
}

TEST_CASE("cell_round_trips_through_IOstreams")
{
    List<highOrderFit::face> faces(1);
    faces[0].setSize(1);
    faces[0][0] = point(2, 3, 4);
    const point centreNotSerialised;

    const highOrderFit::cell cellOut(faces, centreNotSerialised);

    OStringStream o;
    o << cellOut << endl;

    IStringStream i(o.str());
    const highOrderFit::cell cellIn(i);

    REQUIRE( cellIn.size() == 1 );
    REQUIRE( cellIn[0].size() == 1 );
    checkEqual(cellIn[0][0], point(2, 3, 4));
}

TEST_CASE("cell_calculates_zeroth_volume_moment_for_unit_cube_at_origin")
{
    IFstream is("resources/unitCube");
    const highOrderFit::cell cell(is);

    CHECK( cell.moment(highOrderFit::order(0, 0, 0)) == approx(1.0) );
}

TEST_CASE("cell_preserves_right_hand_orientation_when_neighbour_of_face")
{
    const Test::testCase c("cartesian4x3Mesh");
    const Test::mesh testMesh(c.mesh());

    const label celli = testMesh.indexOfCellWithCentreAt(point(2.5, 0.5, 0));
    const highOrderFit::cell cell(c.mesh(), celli);
    CHECK( cell.moment(highOrderFit::order(0, 0, 0)) == approx(1.0) );
}

TEST_CASE("cell_calculates_zeroth_volume_moment_for_2x2x2_cube")
{
    IFstream is("resources/2x2x2Cube");
    const highOrderFit::cell cell(is);

    CHECK( cell.moment(highOrderFit::order(0, 0, 0)) == approx(8.0) );
}

TEST_CASE("cell_calculates_x_volume_moment_for_unit_cube_centred_at_origin")
{
    IFstream is("resources/unitCube");
    const highOrderFit::cell cell(is);

    CHECK( cell.moment(highOrderFit::order(1, 0, 0)) == approx(0.0) );
}

TEST_CASE("cell_calculates_x_volume_moment_for_unit_cube_centred_at_0.5,0,0")
{
    IFstream is("resources/unitCube");
    highOrderFit::cell cell(is);
    translate(cell, vector(0.5, 0, 0));

    CHECK( cell.moment(highOrderFit::order(1, 0, 0)) == approx(0.5) );
}

TEST_CASE("cell_calculates_x^2_volume_moment_for_unit_cube_centred_at_0.5,0,0")
{
    IFstream is("resources/unitCube");
    highOrderFit::cell cell(is);
    translate(cell, vector(0.5, 0, 0));

    CHECK( cell.moment(highOrderFit::order(2, 0, 0)) == approx(1.0/3.0) );
}

TEST_CASE("cell_calculates_x^2_volume_moment_for_unit_cube_centred_at_-1.5,0,0")
{
    IFstream is("resources/unitCube");
    highOrderFit::cell cell(is);
    translate(cell, vector(-1.5, 0, 0));

    CHECK( cell.moment(highOrderFit::order(2, 0, 0)) == approx(7.0/3.0) );
}

TEST_CASE("cell_calculates_x^2_volume_moment_for_unit_cube_centred_at_1.5,0,0")
{
    IFstream is("resources/unitCube");
    highOrderFit::cell cell(is);
    translate(cell, vector(1.5, 0, 0));

    CHECK( cell.moment(highOrderFit::order(2, 0, 0)) == approx(7.0/3.0) );
}

TEST_CASE("cell_calculates_x^2_volume_moment_for_unit_cube_centred_at_-2.5,0,0")
{
    IFstream is("resources/unitCube");
    highOrderFit::cell cell(is);
    translate(cell, vector(-2.5, 0, 0));

    CHECK( cell.moment(highOrderFit::order(2, 0, 0)) == approx(19.0/3.0) );
}

TEST_CASE("cell_calculates_x^2_volume_moment_for_unit_cube_centred_at_2.5,0,0")
{
    IFstream is("resources/unitCube");
    highOrderFit::cell cell(is);
    translate(cell, vector(2.5, 0, 0));

    CHECK( cell.moment(highOrderFit::order(2, 0, 0)) == approx(19.0/3.0) );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

// ************************************************************************* //
