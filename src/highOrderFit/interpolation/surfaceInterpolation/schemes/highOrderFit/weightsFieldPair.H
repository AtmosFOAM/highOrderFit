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
    Foam::highOrderFit::weightsFieldPair

Description
    Contains \ref weightsField "weightsFields" for owner and neighbour sides.

    A single weightsFieldsPair instance can be associated with an fvMesh
    using the MeshObject run-time selector.

SourceFiles
    weightsFieldPair.C

\*---------------------------------------------------------------------------*/

#ifndef weightsFieldPair_H
#define weightsFieldPair_H

#include "extendedUpwindCellToFaceStencil.H"
#include "fvMesh.H"
#include "MeshObject.H"
#include "order.H"
#include "transformer.H"
#include "weightsField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace highOrderFit
{

/*---------------------------------------------------------------------------*\
                         Class weightsFieldPair Declaration
\*---------------------------------------------------------------------------*/

class weightsFieldPair
:
    public MeshObject<fvMesh, MoveableMeshObject, weightsFieldPair>
{
    // Private data

        const weightsField owner_;
        const weightsField neighbour_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        weightsFieldPair(const weightsFieldPair&);

        //- Disallow default bitwise assignment
        void operator=(const weightsFieldPair&);


public:

    // Constructors

        //- Construct a pair of \ref weightsField "weightsFields" for the given mesh, stencil, and configuration parameters.
        weightsFieldPair
        (
            const fvMesh& mesh,
            const extendedUpwindCellToFaceStencil& stencils,
            const transformer& transformer,
            const List<order>& moments,
            const multipliers& multipliers
        );


    //- Destructor
    ~weightsFieldPair();


    // Member Functions
    
        //- Obtain a reference to the owner side weightsField.
        const weightsField& owner() const
        {
            return owner_;
        }
    
        //- Obtain a reference to the neighbour side weightsField.
        const weightsField& neighbour() const
        {
            return neighbour_;
        }

        //- Return the appropriately-sided weightsField depending on \t isOwner.
        const weightsField& side(const bool isOwner) const
        {
            return isOwner ? owner_ : neighbour_;
        }

        //- Does nothing since moving meshes are not currently supported.
        bool movePoints();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
