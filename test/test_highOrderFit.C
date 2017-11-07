#include "catch.hpp"
#include "checks.H"
#include "mesh.H"
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

TEST_CASE("highOrderFit_exactly_reconstructs_linear_in_x_for_vertical_face",
          "[!mayfail]")
{
    Test::interpolation highOrderFit("cartesian4x3Mesh");
    const Test::mesh testMesh(highOrderFit.mesh());
    forAll(highOrderFit.T(), cellI)
    {
        highOrderFit.T()[cellI] = 3*highOrderFit.mesh().C()[cellI].x() + 4;
    }

    const tmp<surfaceScalarField> Tf = highOrderFit.interpolateT();

    const label faceI = testMesh.indexOfFaceWithCentreAt(point(3, 1.5, 0));
    CHECK(Tf()[faceI] == Test::approx(13.0));
}

}
