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
#include "targetFace.H"
#include "testCase.H"

using namespace Foam;

namespace Test
{

TEST_CASE("stencil_is_translated_such_that_targetCf_is_coordinate_origin")
{
    const highOrderFit::targetFace targetFace
    (
        point(2, 1, 0),
        vector(1, 0, 0),
        highOrderFit::face({})
    );

    List<highOrderFit::cell> cells(1);
    List<highOrderFit::face> faces(1, highOrderFit::face({point(0, 0, 0)}));
    cells[0] = highOrderFit::cell(faces, faces[0][0]);

    const highOrderFit::stencil stencil(targetFace, cells);

    checkEqual(stencil[0][0][0], point(-2, -1, 0));
    checkEqual(stencil[0].centre(), point(-2, -1, 0));
}

TEST_CASE("stencil_is_rotated_such_that_primary_direction_is_downwind")
{
    const highOrderFit::targetFace targetFace
    (
        point(0, 0, 0),
        vector(0, 3, 0),
        highOrderFit::face({})
    );

    List<highOrderFit::cell> cells(2);
    List<highOrderFit::face> upwindFaces
    (
        1,
        highOrderFit::face({point(0, -1, 0)})
    );
    List<highOrderFit::face> downwindFaces
    (
        1,
        highOrderFit::face({point(0, 2, 0)})
    );
    cells[0] = highOrderFit::cell(upwindFaces, upwindFaces[0][0]);
    cells[1] = highOrderFit::cell(downwindFaces, downwindFaces[0][0]);

    const highOrderFit::stencil stencil(targetFace, cells);

    checkEqual(stencil[0][0][0], point(-1, 0, 0));
    checkEqual(stencil[0].centre(), point(-1, 0, 0));
    
    checkEqual(stencil[1][0][0], point(2, 0, 0));
    checkEqual(stencil[1].centre(), point(2, 0, 0));
}

TEST_CASE("stencil_calculates_zeroth_moment_of_unit_square_target_face")
{
    const highOrderFit::targetFace targetFace
    (
        point(0, 0, 0),
        vector(1, 0, 0),
        highOrderFit::face
        (
            {
                point(-0.5, -0.5, 0),
                point( 0.5, -0.5, 0),
                point( 0.5,  0.5, 0),
                point(-0.5,  0.5, 0)
            }
        )
    );

    const List<highOrderFit::cell> cells;

    const highOrderFit::stencil stencil(targetFace, cells);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

// ************************************************************************* //
