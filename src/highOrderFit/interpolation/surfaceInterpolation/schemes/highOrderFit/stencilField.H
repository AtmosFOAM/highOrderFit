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
    Foam::stencilField

Description
    Stencils of cells for cellToFace surface interpolation indexed by face.

SourceFiles
    stencilField.C

\*---------------------------------------------------------------------------*/

#ifndef stencilField_H
#define stencilField_H

#include "fvMesh.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class stencilField Declaration
\*---------------------------------------------------------------------------*/

class stencilField
:
    public labelListList
{
    // Private data

        const fvMesh& mesh_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        stencilField(const stencilField&);

        //- Disallow default bitwise assignment
        void operator=(const stencilField&);


public:

    // Constructors

        //- Construct from a specified list of stencils of cell labels.
        stencilField
        (
            const labelListList& stencilCellsList,
            const fvMesh& mesh
        );


    //- Destructor
    ~stencilField();


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
