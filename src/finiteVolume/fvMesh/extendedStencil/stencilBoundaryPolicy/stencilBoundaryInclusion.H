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
    Foam::stencilBoundaryInclusion

Description

SourceFiles
    stencilBoundaryInclusion.C

\*---------------------------------------------------------------------------*/

#ifndef stencilBoundaryInclusion_H
#define stencilBoundaryInclusion_H

#include "polyMesh.H"
#include "stencilBoundaryPolicy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class stencilBoundaryInclusion Declaration
\*---------------------------------------------------------------------------*/

class stencilBoundaryInclusion
:
    public stencilBoundaryPolicy
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        stencilBoundaryInclusion(const stencilBoundaryInclusion&);

        //- Disallow default bitwise assignment
        void operator=(const stencilBoundaryInclusion&);


public:

    // Constructors

        stencilBoundaryInclusion();


    //- Destructor
    ~stencilBoundaryInclusion();


    // Member Functions

        void applyTo
        (
            const polyMesh& mesh,
            boolList& includedBoundaryFaces
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
