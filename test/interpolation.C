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

#include "interpolation.H"
#include "stencilBoundaryExclusion.H"
#include "tmp.H"
#include "upwindCPCCellToFaceStencilObject.H"

using namespace Foam;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const tmp<Foam::surfaceScalarField>
Test::interpolation::initialiseFaceFlux()
{
    const surfaceVectorField Uf
    (
        IOobject
        (
            "Uf",
            c_.runTime().constant(),
            c_.mesh(),
            IOobject::MUST_READ
        ),
        c_.mesh()
    );

    const tmp<surfaceScalarField> tFaceFlux
    (
        new surfaceScalarField(Uf & c_.mesh().Sf())
    );

    return tFaceFlux;
}

const tmp<surfaceInterpolationScheme<scalar>>
Test::interpolation::initialiseScheme
(
    const word& schemeName
)
{
    IStringStream schemeSpecification(schemeName + " ((0 0 0)) " + 
            "uniformMultipliers");
    return surfaceInterpolationScheme<scalar>::New
        (
            c_.mesh(),
            faceFlux_,
            schemeSpecification
        );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Test::interpolation::interpolation
(
    const Foam::fileName& caseName,
    const word& schemeName
)
:
c_(caseName),
faceFlux_(initialiseFaceFlux()),
scheme_(initialiseScheme(schemeName)),
stencils_
(
    Foam::upwindCPCCellToFaceStencilObject::New
    (
        mesh(),
        false, // not pureUpwind
        0.5,   // minOpposedness
        false, // not linearCorr
        Foam::autoPtr<Foam::stencilBoundaryPolicy>
        (
            new Foam::stencilBoundaryExclusion()
        )
    )
),
T_
(
    IOobject
    (
        "T",
        c_.runTime().timeName(),
        c_.mesh()
    ),
    c_.mesh(),
    dimensionedScalar("T", dimless, scalar(0))
)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::tmp<Foam::surfaceScalarField>
Test::interpolation::interpolateT() const
{
    return scheme_().interpolate(T_);
}

void Test::interpolation::negateFaceFlux()
{
    faceFlux_.ref() == -faceFlux_();
}

void Test::interpolation::setTlinearInX(scalar m, scalar c)
{
    forAll(T_, cellI)
    {
        T_[cellI] = m*c_.mesh().C()[cellI].x() + c;
    }
}

// ************************************************************************* //
