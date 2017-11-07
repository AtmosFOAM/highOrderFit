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
    Foam::testMesh

Description
    Utility methods for tests that use fvMesh instances

SourceFiles
    testMesh.C

\*---------------------------------------------------------------------------*/

#ifndef testMesh_H
#define testMesh_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Test
{

/*---------------------------------------------------------------------------*\
                         Class testMesh Declaration
\*---------------------------------------------------------------------------*/

class mesh
{
    // Private data

        const Foam::fvMesh& mesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        mesh(const mesh&);

        //- Disallow default bitwise assignment
        void operator=(const mesh&);


public:

    // Constructors

        //- Construct from components
        mesh(const Foam::fvMesh& mesh);


    // Member Functions

        // Access
        Foam::label indexOfFaceWithCentreAt
        (
            const Foam::point& Cf,
            const Foam::scalar epsilon = Foam::SMALL
        ) const;
        // Check

        // Edit

        // Write

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
