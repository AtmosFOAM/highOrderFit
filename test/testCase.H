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
    Foam::testCase

Description

SourceFiles
    testCase.C

\*---------------------------------------------------------------------------*/

#ifndef testCase_H
#define testCase_H

#include "fileName.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Test
{

/*---------------------------------------------------------------------------*\
                       Class testCase Declaration
\*---------------------------------------------------------------------------*/

class testCase
{
    // Private data

        const Foam::Time runTime_;
        const Foam::fvMesh mesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        testCase(const testCase&);

        //- Disallow default bitwise assignment
        void operator=(const testCase&);


public:

    // Constructors

        testCase(const Foam::fileName& name);


    //- Destructor
    ~testCase();


    // Member Functions

        const Foam::Time& runTime() const
        {
            return runTime_;
        };

        const Foam::fvMesh& mesh() const
        {
            return mesh_;
        };


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
