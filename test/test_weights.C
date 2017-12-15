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
#include "mesh.H"
#include "inverseDistanceMultipliers.H"
#include "IStringStream.H"
#include "order.H"
#include "stencil.H"
#include "testCase.H"
#include "uniformMultipliers.H"
#include "weights.H"

using namespace Foam;

namespace Test
{

TEST_CASE("weights_with_uniform_multipliers_average_all_cells_in_stencil")
{
    const Test::testCase c("cartesian4x3Mesh");
    const Test::mesh testMesh(c.mesh());
    const label facei = testMesh.indexOfFaceWithCentreAt(point(3, 1.5, 0));
    
    const label uCelli = testMesh.indexOfCellWithCentreAt(point(2.5, 1.5, 0));
    const label dCelli = testMesh.indexOfCellWithCentreAt(point(3.5, 1.5, 0));

    List<highOrderFit::cell> cells(12);
    label i = 0;
    cells[i++] = highOrderFit::cell(c.mesh(), uCelli);
    cells[i++] = highOrderFit::cell(c.mesh(), dCelli);
    for (label celli = 0; celli < 12; celli++)
    {
        if (celli != uCelli && celli != dCelli)
        {
            cells[i++] = highOrderFit::cell(c.mesh(), celli);
        }
    }

    const highOrderFit::targetFace targetFace(c.mesh(), facei);
    const highOrderFit::stencil stencil(targetFace, cells);
    
    scalarList w(12);
    const scalarList multipliers(12, 1.0);
    const highOrderFit::weights weights({highOrderFit::order(0, 0, 0)});

    weights.calculate(w, stencil, multipliers);

    checkEqual( w, 1.0/12.0 );
}

TEST_CASE("weights_with_inverse_distance_multipliers_fit_central_cells_closely")
{
    const Test::testCase c("cartesian4x3Mesh");
    const Test::mesh testMesh(c.mesh());
    const label facei = testMesh.indexOfFaceWithCentreAt(point(3, 1.5, 0));
    
    const label uCelli = testMesh.indexOfCellWithCentreAt(point(2.5, 1.5, 0));
    const label dCelli = testMesh.indexOfCellWithCentreAt(point(3.5, 1.5, 0));
    const label uuCelli = testMesh.indexOfCellWithCentreAt(point(1.5, 1.5, 0));

    List<highOrderFit::cell> cells(3);
    cells[0] = highOrderFit::cell(c.mesh(), uCelli);
    cells[1] = highOrderFit::cell(c.mesh(), dCelli);
    cells[2] = highOrderFit::cell(c.mesh(), uuCelli);

    const highOrderFit::targetFace targetFace(c.mesh(), facei);
    const highOrderFit::stencil stencil(targetFace, cells);

    scalarList w(3);
    const highOrderFit::inverseDistanceMultipliers multipliers;
    scalarList m(3);
    multipliers.calculate(stencil, m);
    const highOrderFit::weights weights({highOrderFit::order(0, 0, 0)});

    weights.calculate(w, stencil, m);

    const label upwind = 0, downwind = 1;
    CHECK( w[upwind] == approx(0.4999992847) ); 
    CHECK( w[downwind] == approx(0.4999992847) ); 
}

TEST_CASE("weights_populates_matrix_with_zeroth_and_x_volume_moments")
{
    const Test::testCase c("cartesian4x3Mesh");
    const Test::mesh testMesh(c.mesh());
    const label facei = testMesh.indexOfFaceWithCentreAt(point(3, 1.5, 0));
    
    const label uCelli = testMesh.indexOfCellWithCentreAt(point(2.5, 1.5, 0));
    const label dCelli = testMesh.indexOfCellWithCentreAt(point(3.5, 1.5, 0));
    const label uuCelli = testMesh.indexOfCellWithCentreAt(point(1.5, 1.5, 0));

    const List<highOrderFit::order> moments(
    {
        highOrderFit::order(0, 0, 0),
        highOrderFit::order(1, 0, 0)
    });

    List<highOrderFit::cell> cells(3);
    cells[0] = highOrderFit::cell(c.mesh(), uCelli);
    cells[1] = highOrderFit::cell(c.mesh(), dCelli);
    cells[2] = highOrderFit::cell(c.mesh(), uuCelli);

    const highOrderFit::targetFace targetFace(c.mesh(), facei);
    const highOrderFit::stencil stencil(targetFace, cells);

    const scalarList multipliers(3, 1.0);
    const highOrderFit::weights weights(moments);

    autoPtr<scalarRectangularMatrix> actual =
        weights.createMatrix(stencil, multipliers);

    IStringStream in("3 2((1 -0.5)(1 0.5)(1 -1.5))");
    scalarRectangularMatrix expected(in);
    checkEqual(actual(), expected);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

// ************************************************************************* //
