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

#include "stencilBoundaryInclusion.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stencilBoundaryInclusion::stencilBoundaryInclusion()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::stencilBoundaryInclusion::~stencilBoundaryInclusion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::stencilBoundaryInclusion::applyTo
(
    const Foam::polyMesh& mesh,
    Foam::boolList& includedBoundaryFaces
) const
{
	const polyBoundaryMesh& patches = mesh.boundaryMesh();

    includedBoundaryFaces.setSize
    (
        mesh.nFaces() - mesh.nInternalFaces(),
        true
    );

	forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled() || isA<emptyPolyPatch>(pp))
        {
            label bFacei = pp.start()-mesh.nInternalFaces();
            forAll(pp, i)
            {                                                                   
                includedBoundaryFaces[bFacei++] = false;
            }
        }
    }

}

// ************************************************************************* //
