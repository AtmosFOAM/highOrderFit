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

Class
    Test::interpolation

Description
    Convenient test setup for testing the highOrderFit surfaceInterpolationScheme

SourceFiles
    interpolation.C

\*---------------------------------------------------------------------------*/

#ifndef interpolation_H
#define interpolation_H

#include "scalar.H"
#include "surfaceInterpolationScheme.H"
#include "fileName.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Test
{

/*---------------------------------------------------------------------------*\
                         Class interpolation Declaration
\*---------------------------------------------------------------------------*/

class interpolation
{
    // Private data

        const Foam::Time runTime_;
        const Foam::fvMesh mesh_;
        const Foam::tmp<Foam::surfaceScalarField> faceFlux_;
        const Foam::tmp
        <
            Foam::surfaceInterpolationScheme<Foam::scalar>
        >
        scheme_;
        Foam::volScalarField T_;

    // Private Member Functions

        const Foam::tmp<Foam::surfaceScalarField> initialiseFaceFlux();

        const Foam::tmp
        <
            Foam::surfaceInterpolationScheme<Foam::scalar>
        >
        initialiseScheme(const Foam::word& schemeName);

        //- Disallow default bitwise copy construct
        interpolation(const interpolation&);

        //- Disallow default bitwise assignment
        void operator=(const interpolation&);

public:

    // Constructors

        interpolation
        (
            const Foam::fileName& caseName,
            const Foam::word& schemeName = "highOrderFit"
        );

    // Member Functions

        const Foam::Time& runTime() const;
        const Foam::fvMesh& mesh() const;
        Foam::volScalarField& T();
        const Foam::tmp<Foam::surfaceScalarField> interpolateT() const;

        //- reverse the velocity field.  This is useful for testing owner
        //- and neighbour sides for the upwind-biased stencil.
        void negateFaceFlux();

        void setTlinearInX(Foam::scalar m=3, Foam::scalar c=4);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Test

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
