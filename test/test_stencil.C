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
    List<highOrderFit::cellVertices> stencilCellVertices;
    const highOrderFit::stencil stencil(targetCf, Sf, stencilCellVertices);

    const scalarList& zerothCellMoments =
        stencil.moment(/*highOrderFit::order(0, 0, 0)*/);

    checkEqual(zerothCellMoments, 1.0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

// ************************************************************************* //
