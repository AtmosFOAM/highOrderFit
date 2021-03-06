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
#include "IStringStream.H"
#include "order.H"
#include "OStringStream.H"

using namespace Foam;

namespace Test
{

TEST_CASE("order_calculates_factorial_ratio")
{
    const highOrderFit::order order(2, 3, 4);
    CHECK( order.factorialRatio(3) == 288.0/479001600.0 );
}

TEST_CASE("order_calculates_exponent_tensor_for_zeroth_moment")
{
    const highOrderFit::order order(0, 0, 0);
    List<highOrderFit::exponentTensor> K;
    order.calculateExponentTensors(K);

    IStringStream in("(0 0 0 0 0 0 0 0 0)");
    const highOrderFit::exponentTensor expectedExponentTensor(in);

    REQUIRE( K.size() == 1);
    checkEqual( K[0], expectedExponentTensor );
}

TEST_CASE("order_calculates_exponent_tensor_for_1_0_0_moment")
{
    const highOrderFit::order order(1, 0, 0);
    List<highOrderFit::exponentTensor> K;
    order.calculateExponentTensors(K);

    IStringStream in(R"END(((0 0 1 0 0 0 0 0 0)
                            (0 1 0 0 0 0 0 0 0)
                            (1 0 0 0 0 0 0 0 0)))END");
    const List<highOrderFit::exponentTensor> expectedExponentTensors(in);

    checkEqual
    (
        static_cast<List<labelTensor>>(K),
        static_cast<List<labelTensor>>(expectedExponentTensors)
    );
}

TEST_CASE("order_round_trips_through_IOstreams")
{
    const highOrderFit::order orderOut(2, 3, 4);

    OStringStream o;
    o << orderOut << endl;

    IStringStream i(o.str());
    const highOrderFit::order orderIn(i);

    CHECK( orderIn == orderOut );
}

TEST_CASE("order_round_trips_list_through_IOstreams")
{
    List<highOrderFit::order> orderOut(1);
    orderOut[0] = highOrderFit::order(2, 3, 4);

    OStringStream o;
    o << orderOut << endl;

    IStringStream i(o.str());
    List<highOrderFit::order> orderIn;

    checkEqual(orderIn, orderOut);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

// ************************************************************************* //
