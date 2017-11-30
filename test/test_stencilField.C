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
#include "labelList.H"
#include "stencilField.H"
#include "testCase.H"

using namespace Foam;

namespace Test
{

TEST_CASE("stencilField_creates_stencil_for_each_in_stencilCellsList")
{
    const Test::testCase c("cartesian4x3Mesh");
    const labelListList stencilCellsList(5);

    highOrderFit::stencilField stencils(stencilCellsList, c.mesh());

    CHECK( stencils.size() == 5 );
}

TEST_CASE("stencilField_creates_each_stencil_with_size_of_stencilCells")
{
    const Test::testCase c("cartesian4x3Mesh");
    labelListList stencilCellsList(1);
    stencilCellsList[0].setSize(3);

    highOrderFit::stencilField stencils(stencilCellsList, c.mesh());

    CHECK( stencils[0].size() == 3 );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

// ************************************************************************* //
