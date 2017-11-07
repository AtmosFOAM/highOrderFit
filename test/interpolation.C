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
#include "tmp.H"

using namespace Foam;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const tmp<Foam::surfaceScalarField>
Test::interpolation::initialisePhi()
{
    const surfaceVectorField Uf
    (
        IOobject
        (
            "Uf",
            runTime_.constant(),
            mesh_,
            IOobject::MUST_READ
        ),
        mesh_
    );

    const tmp<surfaceScalarField> tPhi
    (
        new surfaceScalarField(Uf & mesh_.Sf())
    );

    return tPhi;
}

const tmp<surfaceInterpolationScheme<scalar>>
Test::interpolation::initialiseScheme
(
    const word& schemeName
)
{
    IStringStream schemeNameStream(schemeName);
    return surfaceInterpolationScheme<scalar>::New
        (
            mesh_,
            phi_,
            schemeNameStream
        );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Test::interpolation::interpolation
(
    const fileName& caseName,
    const word& schemeName
)
:
runTime_(
    Time::controlDictName,
    "resources",
    caseName
),
mesh_
(
    IOobject
    (
        fvMesh::defaultRegion,
        runTime_.constant(),
        runTime_,
        IOobject::MUST_READ
    )
),
phi_(initialisePhi()),
scheme_(initialiseScheme(schemeName)),
T_
(
    IOobject
    (
        "T",
        runTime_.timeName(),
        mesh_
    ),
    mesh_,
    dimensionedScalar("T", dimless, scalar(0))
)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::Time& Test::interpolation::runTime() const
{
    return runTime_;
}

const Foam::fvMesh& Test::interpolation::mesh() const
{
    return mesh_;
}

Foam::volScalarField& Test::interpolation::T()
{
    return T_;
}

const Foam::tmp<Foam::surfaceScalarField>
Test::interpolation::interpolateT() const
{
    return scheme_().interpolate(T_);
}

// ************************************************************************* //
