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
#include "interpolation.H"
#include "labelList.H"
#include "mesh.H"
#include "stencil.H"
#include "stencilField.H"
#include "testCase.H"

using namespace Foam;

namespace Test
{

static Approx pointApprox = Approx::custom().margin(1e-12);

TEST_CASE("stencilField_creates_stencil_for_each_in_stencilCellsList")
{
    const Test::interpolation highOrderFit("cartesian4x3Mesh");

    highOrderFit::stencilField stencilField
    (
        highOrderFit.stencils().ownStencil(),
        highOrderFit.stencils().ownMap(),
        highOrderFit.mesh()
    );

    CHECK( stencilField.size() == 55 );
}

TEST_CASE("stencilField_populates_stencil_cells")
{
    const Test::interpolation highOrderFit("cartesian4x3Mesh");
    const Test::mesh testMesh(highOrderFit.mesh());

    highOrderFit::stencilField stencilField
    (
        highOrderFit.stencils().ownStencil(),
        highOrderFit.stencils().ownMap(),
        highOrderFit.mesh()
    );

    const label facei = testMesh.indexOfFaceWithCentreAt(point(3, 1.5, 0));
    const highOrderFit::stencil& stencil = stencilField[facei];

    REQUIRE( stencil.size() == 12 );
    REQUIRE( stencil[0].size() == 6); // 6 faces
    CHECK( stencil[0][0].size() == 4); // 4 points
}

TEST_CASE("stencilField_translate_stencil_with_targetCf_as_origin")
{
    Test::interpolation highOrderFit("cartesian4x3Mesh");
    const Test::mesh testMesh(highOrderFit.mesh());

    highOrderFit::stencilField stencilField
    (
        highOrderFit.stencils().ownStencil(),
        highOrderFit.stencils().ownMap(),
        highOrderFit.mesh()
    );

    const label facei = testMesh.indexOfFaceWithCentreAt(point(3, 1.5, 0));
    const highOrderFit::stencil& stencil = stencilField[facei];

    REQUIRE( stencil.size() == 12 );
    const List<highOrderFit::face>& upwindCell = stencil[0];

    forAll(upwindCell, i)
    {
        forAll(upwindCell[i], pointi)
        {
            CHECK( upwindCell[i][pointi].x() <= 0 );
        }
    }
}

TEST_CASE("stencilField_rotates_stencil_for_horizontal_face")
{
    Test::interpolation highOrderFit("cartesian3x4Mesh");
    const Test::mesh testMesh(highOrderFit.mesh());

    highOrderFit::stencilField stencilField
    (
        highOrderFit.stencils().ownStencil(),
        highOrderFit.stencils().ownMap(),
        highOrderFit.mesh()
    );

    const label facei = testMesh.indexOfFaceWithCentreAt(point(1.5, 3, 0));
    const highOrderFit::stencil& stencil = stencilField[facei];

    REQUIRE( stencil.size() == 12 );
    const List<highOrderFit::face>& upwindCell = stencil[0];

    forAll(upwindCell, i)
    {
        forAll(upwindCell[i], pointi)
        {
            CHECK( upwindCell[i][pointi].x() >= pointApprox(-1.0) );
            CHECK( upwindCell[i][pointi].x() <= pointApprox(0.0) );
            CHECK( upwindCell[i][pointi].y() >= pointApprox(-0.5) );
            CHECK( upwindCell[i][pointi].y() <= pointApprox(0.5) );
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

// ************************************************************************* //
