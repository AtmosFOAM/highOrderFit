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
#include "cellVertices.H"
#include "checks.H"
#include "inverseDistanceMultipliers.H"
#include "stencil.H"
#include "uniformMultipliers.H"
#include "weights.H"

using namespace Foam;

namespace Test
{

TEST_CASE("weights_with_uniform_multipliers_average_all_cells_in_stencil")
{
    scalarList w(12);
    const point targetFace;
    const vector Sf;
    const List<highOrderFit::cellVertices> stencilCellVertices(12);
    highOrderFit::stencil stencil(targetFace, Sf, stencilCellVertices);
    highOrderFit::uniformMultipliers multipliers(12);
    highOrderFit::weights weights;

    weights.calculate(w, stencil, multipliers);

    checkEqual( w, 1.0/12.0 );
}

TEST_CASE("weights_with_inverse_distance_multipliers_fit_central_cells_closely")
{
    scalarList w(5);
    const point targetFace;
    const vector Sf;
    const List<highOrderFit::cellVertices> stencilCellVertices(5);
    highOrderFit::stencil stencil(targetFace, Sf, stencilCellVertices);
    highOrderFit::inverseDistanceMultipliers multipliers(5);
    highOrderFit::weights weights;

    weights.calculate(w, stencil, multipliers);

    const label upwind = 0, downwind = 1;
    CHECK( w[upwind] == approx(0.4999992847) ); 
    CHECK( w[downwind] == approx(0.4999992847) ); 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

// ************************************************************************* //
