#include "catch.hpp"
#include "checks.H"
#include "interpolation.H"

#include "IOobject.H"
#include "tmp.H"

using namespace Foam;

namespace Test
{

TEST_CASE("highOrderFit_interpolates_constant_scalar_field")
{
    Test::interpolation highOrderFit("cartesian4x3Mesh");
    highOrderFit.T() = dimensionedScalar("T", dimless, scalar(1));
    const surfaceScalarField expectedTf
    (
        IOobject
        (
            "expectedTf",
            highOrderFit.runTime().timeName(),
            highOrderFit.mesh()
        ),
        highOrderFit.mesh(),
        dimensionedScalar("expectedTf", dimless, scalar(1))
    );

    const tmp<surfaceScalarField> Tf = highOrderFit.interpolateT();

    Test::checkEqual(Tf, expectedTf);
}

}
