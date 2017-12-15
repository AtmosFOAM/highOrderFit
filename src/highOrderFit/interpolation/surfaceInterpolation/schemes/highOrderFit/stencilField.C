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

\*---------------------------------------------------------------------------*/

#include "stencilField.H"
#include "extendedCellToFaceStencil.H"
#include "targetFace.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::highOrderFit::stencilField::stencilField
(
    const Foam::labelListList& stencilCellsList,
    const Foam::mapDistribute& map,
    const Foam::fvMesh& mesh
)
:
    List<stencil>(stencilCellsList.size()),
    mesh_(mesh)
{
    List<cell> myCells(map.constructSize());
    forAll(mesh_.cells(), celli)
    {
        myCells[celli] = cell(mesh_, celli);
    }

    map.distribute(myCells);

    forAll(stencilCellsList, facei)
    {
        const labelList& stencilCells = stencilCellsList[facei];
        List<cell> cells(stencilCells.size());

        forAll(stencilCells, i)
        {
            cells[i] = myCells[stencilCells[i]];
        }

        const targetFace targetFace(mesh_, facei);
        (*this)[facei] = stencil(targetFace, cells);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highOrderFit::stencilField::~stencilField()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
