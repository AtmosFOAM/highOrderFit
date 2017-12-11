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
#include "List.H"
#include "mesh.H"
#include "order.H"
#include "stencil.H"
#include "testCase.H"

using namespace Foam;

namespace Test
{

TEST_CASE("stencil_is_translated_such_that_targetCf_is_coordinate_origin")
{
    const point targetCf(2, 1, 0);
    const vector Sf(1, 0, 0);

    List<highOrderFit::cell> cells(1);
    List<List<point>> points(1, List<point>({point(0, 0, 0)}));
    cells[0] = highOrderFit::cell(points, points[0][0]);

    const highOrderFit::stencil stencil(targetCf, Sf, cells);

    checkEqual(stencil[0][0][0], point(-2, -1, 0));
    checkEqual(stencil[0].centre(), point(-2, -1, 0));
}

TEST_CASE("stencil_is_rotated_such_that_primary_direction_is_downwind")
{
    const point targetCf(0, 0, 0);
    const vector Sf(0, 3, 0);
    List<highOrderFit::cell> cells(2);
    List<List<point>> upwindPoints(1, List<point>({point(0, -1, 0)}));
    List<List<point>> downwindPoints(1, List<point>({point(0, 2, 0)}));
    cells[0] = highOrderFit::cell(upwindPoints, upwindPoints[0][0]);
    cells[1] = highOrderFit::cell(downwindPoints, downwindPoints[0][0]);

    const highOrderFit::stencil stencil(targetCf, Sf, cells);

    checkEqual(stencil[0][0][0], point(-1, 0, 0));
    checkEqual(stencil[1][0][0], point(2, 0, 0));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

// ************************************************************************* //
