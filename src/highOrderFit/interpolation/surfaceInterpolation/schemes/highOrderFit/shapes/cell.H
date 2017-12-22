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
    A cell belonging to a stencil, comprising a collection of faces.

SourceFiles
    cell.C

\*---------------------------------------------------------------------------*/

#ifndef highOrderFit_cell_H
#define highOrderFit_cell_H

#include "face.H"
#include "fvMesh.H"
#include "order.H"
#include "IOstreams.H"
#include "List.H"
#include "point.H"
#include "volFields.H"

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
    public List<face>
{
    // Private data

        point C_;

public:

    // Constructors

        //- Construct an instance for a given cell.
        //  @param mesh a mesh instance
        //  @param celli an index into <tt>mesh.cells()</tt>
        cell(const fvMesh& mesh, const label celli);

        //- Construct directly from a collection of faces and specified centre.
        cell(const List<face>& faces, const point centre);

        //- Construct by reading from an input stream.
        //  The cell centre is assumed to be <tt>(0, 0, 0)</tt>.
        cell(Istream&);

        //- Construct an empty instance with no faces.
        //  Necessary only for initialising collections of cells.
        cell();


    //- Destructor
    ~cell();


    // Member Functions

        //- Transform the cell into local coordinates.
        void transform(const vector x, const vector from, const vector to);

        //- The volume moment of the given order.
        //  Volume moments are calculated using the method of
        //  <a href="https://doi.org/10.1016/S0031-3203(03)00127-4">
        //  Tuzikov et al. 2003</a>.
        scalar moment(const order& o) const;

        //- The cell centre as specified during construction.
        point centre() const
        {
            return C_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
