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
    Foam::highOrderFitStencilField

Description
    Stencils of cells for cellToFace surface interpolation indexed by face.

SourceFiles
    highOrderFitStencilField.C

\*---------------------------------------------------------------------------*/

#ifndef highOrderFitStencilField_H
#define highOrderFitStencilField_H

#include "fvMesh.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class highOrderFitStencilField Declaration
\*---------------------------------------------------------------------------*/

class highOrderFitStencilField
:
    public labelListList
{
    // Private data

        const fvMesh& mesh_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        highOrderFitStencilField(const highOrderFitStencilField&);

        //- Disallow default bitwise assignment
        void operator=(const highOrderFitStencilField&);


public:

    // Constructors

        //- Construct from a specified list of stencils of cell labels.
        highOrderFitStencilField
        (
            const fvMesh& mesh,
            const labelListList& stencilCellsList
        );


    //- Destructor
    ~highOrderFitStencilField();


    // Member Functions

        const fvMesh& mesh() const
        {
            return mesh_;
        }

        // Check
        // Edit
        // Write


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
