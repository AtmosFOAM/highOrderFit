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
    Foam::highOrderFit::cell

Description
    A cell belonging to a stencil with vertex points stored per face.

SourceFiles
    cellI.H
    cell.C

\*---------------------------------------------------------------------------*/

#ifndef highOrderFit_cell_H
#define highOrderFit_cell_H

#include "order.H"
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
                         Class cell Declaration
\*---------------------------------------------------------------------------*/

class cell
:
    public List<List<point>>
{

    // Private Member Functions

        // Disallow default bitwise copy construct
        //cell(const cell&);

        // Disallow default bitwise assignment
        //void operator=(const cell&);

public:

    // Constructors

        //- Construct an instance for a given cell.
        //  @param mesh a mesh instance
        //  @param celli an index into <tt>mesh.cells()</tt>
        cell(const primitiveMesh& mesh, const label celli);

        //- Construct directly from point data.
        explicit cell(const List<List<point>>& points);

        //- Construct an empty instance with no vertex points.
        //  Necessary only for initialising collections of cells.
        cell();


    //- Destructor
    ~cell();


    // Member Functions

        //- Translate the cell by vector \c x.
        void translate(const vector x);

        //- Rotate the cell to align \c from with \c to.
        //  The transformation is constrained by the quaternion rotation, so
        //  cells cannot be rotated 180&deg; or more.
        void rotate(const vector from, const vector to);

        scalar moment(const order& o) const;

        point average() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
