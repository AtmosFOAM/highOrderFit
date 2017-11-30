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
#include "interpolation.H"
#include "mesh.H"

#include "stencilField.H"
#include "weightsField.H"

using namespace Foam;

namespace Test
{

static void initialiseStencilCellsList(labelListList& stencilCellsList)
{
    forAll(stencilCellsList, i)
    {
        stencilCellsList[i].setSize(2);
    }
}

TEST_CASE("weightsField_has_same_size_as_stencil_field")
{
    Test::interpolation highOrderFit("cartesian4x3Mesh");
    labelListList stencilCellsList(55);
    Test::initialiseStencilCellsList(stencilCellsList);

    highOrderFit::stencilField stencils(stencilCellsList, highOrderFit.mesh());
    highOrderFit::weightsField weights(stencils);

    CHECK( weights.size() == 55 );
}

TEST_CASE("weightsField_has_weight_for_each_cell_in_stencil")
{
    Test::interpolation highOrderFit("cartesian4x3Mesh");
    const Test::mesh testMesh(highOrderFit.mesh());
    const label facei = testMesh.indexOfFaceWithCentreAt(point(3, 1.5, 0));
    
    labelListList stencilCellsList(55);
    Test::initialiseStencilCellsList(stencilCellsList);
    stencilCellsList[facei].setSize(12);

    highOrderFit::stencilField stencils(stencilCellsList, highOrderFit.mesh());
    highOrderFit::weightsField weights(stencils);

    CHECK( weights[facei].size() == 12 );
}

TEST_CASE("weightsField_calculates_correction_on_upwind")
{
    Test::interpolation highOrderFit("cartesian4x3Mesh");
    const Test::mesh testMesh(highOrderFit.mesh());
    const label facei = testMesh.indexOfFaceWithCentreAt(point(3, 1.5, 0));
    
    labelListList stencilCellsList(55);
    Test::initialiseStencilCellsList(stencilCellsList);
    stencilCellsList[facei].setSize(12);

    highOrderFit::stencilField stencils(stencilCellsList, highOrderFit.mesh());
    highOrderFit::weightsField weights(stencils);

    const label upwind = 0;
    CHECK( weights[facei][upwind] == Test::approx(-11.0/12.0) );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

// ************************************************************************* //