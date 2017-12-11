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
#include "cell.H"
#include "checks.H"
#include "inverseDistanceMultipliers.H"
#include "IStringStream.H"
#include "order.H"
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
    const List<highOrderFit::cell> cells(12);
    const highOrderFit::stencil stencil(targetFace, Sf, cells);
    const scalarList multipliers(12, 1.0);
    const highOrderFit::weights weights({highOrderFit::order(0, 0, 0)});

    weights.calculate(w, stencil, multipliers);

    checkEqual( w, 1.0/12.0 );
}

TEST_CASE("weights_with_inverse_distance_multipliers_fit_central_cells_closely")
{
    scalarList w(5);
    const point targetFace;
    const vector Sf;
    const List<highOrderFit::cell> cells(5);
    const highOrderFit::stencil stencil(targetFace, Sf, cells);
    const highOrderFit::inverseDistanceMultipliers multipliers(5);
    const highOrderFit::weights weights({highOrderFit::order(0, 0, 0)});

    weights.calculate(w, stencil, multipliers);

    const label upwind = 0, downwind = 1;
    CHECK( w[upwind] == approx(0.4999992847) ); 
    CHECK( w[downwind] == approx(0.4999992847) ); 
}

TEST_CASE("weights_populates_matrix_with_zeroth_and_x_volume_moments")
{
    const List<highOrderFit::order> moments(
    {
        highOrderFit::order(0, 0, 0),
        highOrderFit::order(1, 0, 0)
    });
    const point targetFace;
    const vector Sf;

    List<List<point>> cell0(1);
    cell0[0] = List<point>({point(-1.5, 0, 0)});
    List<List<point>> cell1(1);
    cell1[0] = List<point>({point(-0.5, 0, 0)});
    List<List<point>> cell2(1);
    cell2[0] = List<point>({point(0.5, 0, 0)});

    List<highOrderFit::cell> cells(3);
    cells[0] = highOrderFit::cell(cell0, cell0[0][0]);
    cells[1] = highOrderFit::cell(cell1, cell1[0][0]);
    cells[2] = highOrderFit::cell(cell2, cell2[0][0]);

    const highOrderFit::stencil stencil(targetFace, Sf, cells);
    const scalarList multipliers(3, 1.0);
    const highOrderFit::weights weights(moments);

    autoPtr<scalarRectangularMatrix> actual =
        weights.createMatrix(stencil, multipliers);

    IStringStream in("3 2((1 -1.5)(1 -0.5)(1 0.5))");
    scalarRectangularMatrix expected(in);
    checkEqual(actual(), expected);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

// ************************************************************************* //
