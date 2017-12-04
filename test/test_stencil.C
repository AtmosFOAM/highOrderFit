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
#include "List.H"
#include "mesh.H"
#include "stencil.H"
#include "testCase.H"

using namespace Foam;

namespace Test
{

TEST_CASE("stencil_calculates_zeroth_cell_moments")
{
    const point targetCf;
    const vector Sf;
    const List<highOrderFit::cellVertices> vertices;
    const highOrderFit::stencil stencil(targetCf, Sf, vertices);

    const scalarList& zerothCellMoments =
        stencil.moment(/*highOrderFit::order(0, 0, 0)*/);

    checkEqual(zerothCellMoments, 1.0);
}

TEST_CASE("stencil_is_translated_such_that_targetCf_is_coordinate_origin")
{
    const point targetCf(2, 1, 0);
    const vector Sf;

    List<highOrderFit::cellVertices> vertices(1);
    List<List<point>> points(1);
    points[0] = List<point>({point(0, 0, 0)});
    vertices[0] = highOrderFit::cellVertices(points);

    const highOrderFit::stencil stencil(targetCf, Sf, vertices);

    checkEqual(stencil.vertices()[0][0][0], point(-2, -1, 0));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

// ************************************************************************* //
