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
    Foam::stencilBoundaryExclusion

Description

SourceFiles
    stencilBoundaryExclusion.C

\*---------------------------------------------------------------------------*/

#ifndef stencilBoundaryExclusion_H
#define stencilBoundaryExclusion_H

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class stencilBoundaryExclusion Declaration
\*---------------------------------------------------------------------------*/

class stencilBoundaryExclusion
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        stencilBoundaryExclusion(const stencilBoundaryExclusion&);

        //- Disallow default bitwise assignment
        void operator=(const stencilBoundaryExclusion&);


public:

    // Constructors

        //- Construct null
        stencilBoundaryExclusion();


    //- Destructor
    ~stencilBoundaryExclusion();


    // Member Functions

        void applyTo
        (
            const polyMesh& mesh,
            boolList& includeBoundaryFaces
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
