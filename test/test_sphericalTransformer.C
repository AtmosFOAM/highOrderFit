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
#include "sphericalTransformer.H"
#include "checks.H"
#include "stencil.H"
#include "targetFace.H"

using namespace Foam;

namespace Test
{

TEST_CASE("sphericalTransformer_translates_stencil_so_targetCf_is_origin")
{
    const highOrderFit::sphericalTransformer transformer;
    const highOrderFit::targetFace targetFace
    (
        {point(10, 11, 12)},
        point(0, 0, 2),
        vector(1, 0, 0)
    );

    const List<highOrderFit::cell> cells
    (
        {
            highOrderFit::cell({highOrderFit::face({point(7, 8, 9)})})
        }
    );

    highOrderFit::stencil stencil(targetFace, cells);

    transformer.transform(stencil);

    checkEqual(stencil[0][0][0], point(7, 8, 7));
    checkEqual(stencil.target()[0], point(10, 11, 10));
}

TEST_CASE("sphericalTransformer_rotates_so_that_primary_direction_downwind")
{
    const highOrderFit::sphericalTransformer transformer;
    const highOrderFit::targetFace targetFace
    (
        {point(1, 1, 0)},
        point(1, 0, 0),
        vector(0, 1, 0)
    );

    const List<highOrderFit::cell> cells
    (
        {
            highOrderFit::cell({highOrderFit::face({point(1, -1, 0)})})
        }
    );

    highOrderFit::stencil stencil(targetFace, cells);

    transformer.transform(stencil);

    checkEqual(stencil[0][0][0], point(-1, 0, 0));
    checkEqual(stencil.target()[0], point(1, 0, 0));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

// ************************************************************************* //
