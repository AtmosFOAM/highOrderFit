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
    Foam::highOrderFit::cellVertices

Description
    Points of the cell vertices per face.

SourceFiles
    cellVertices.C

\*---------------------------------------------------------------------------*/

#ifndef cellVertices_H
#define cellVertices_H

#include "primitiveMesh.H"
#include "List.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace highOrderFit
{

/*---------------------------------------------------------------------------*\
                         Class cellVertices Declaration
\*---------------------------------------------------------------------------*/

class cellVertices
:
    public List<List<point>>
{

    // Private Member Functions

        // Disallow default bitwise copy construct
        //cellVertices(const cellVertices&);

        // Disallow default bitwise assignment
        //void operator=(const cellVertices&);


public:

    // Constructors

        //- Construct an instance for a given cell.
        //  @param mesh a mesh instance
        //  @param celli an index into <code>mesh.cells()</code>
        cellVertices(const primitiveMesh& mesh, const label celli);

        //- Construct an empty instance with no vertices.
        //  Necessary only for initialising collections of cellVertices.
        cellVertices();


    //- Destructor
    ~cellVertices();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
