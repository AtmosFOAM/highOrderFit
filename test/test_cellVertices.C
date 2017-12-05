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
#include "cellVertices.H"
#include "labelList.H"
#include "mesh.H"
#include "testCase.H"

using namespace Foam;

namespace Test
{

TEST_CASE("cellVertices_has_six_faces_for_cuboid_cell")
{
    const Test::testCase c("cartesian4x3Mesh");
    const Test::mesh testMesh(c.mesh());

    const label celli = testMesh.indexOfCellWithCentreAt(point(0.5, 0.5, 0));
    const highOrderFit::cellVertices vertices(c.mesh(), celli);

    CHECK( vertices.size() == 6 );
}

TEST_CASE("cellVertices_has_four_vertices_for_each_face_of_cuboid_cell")
{
    const Test::testCase c("cartesian4x3Mesh");
    const Test::mesh testMesh(c.mesh());

    const label celli = testMesh.indexOfCellWithCentreAt(point(0.5, 0.5, 0));
    const highOrderFit::cellVertices vertices(c.mesh(), celli);

    forAll(vertices, i)
    {
        CHECK( vertices[i].size() == 4 );
    }
}

TEST_CASE("cellVertices_has_vertices_of_face")
{
    const Test::testCase c("cartesian4x3Mesh");
    const Test::mesh testMesh(c.mesh());

    const label celli = testMesh.indexOfCellWithCentreAt(point(0.5, 0.5, 0));
    const highOrderFit::cellVertices vertices(c.mesh(), celli);

    const List<point>& firstFace = vertices[0];
    Test::checkEqual(firstFace[0], point(1, 0, -0.5));
    Test::checkEqual(firstFace[1], point(1, 1, -0.5));
    Test::checkEqual(firstFace[2], point(1, 1,  0.5));
    Test::checkEqual(firstFace[3], point(1, 0,  0.5));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

// ************************************************************************* //
