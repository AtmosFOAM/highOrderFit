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
    Foam::highOrderFit

Description
    A high-order approximation for interpolating surface fluxes from an
    upwind-biased stencil of cell averages.

SourceFiles
    highOrderFit.C

\*---------------------------------------------------------------------------*/

#ifndef highOrderFit_H
#define highOrderFit_H

#include "extendedUpwindCellToFaceStencil.H"
#include "upwind.H"
#include "upwindCPCCellToFaceStencilObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class highOrderFit Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class highOrderFit
:
    public upwind<Type>
{
    // Private data

        //- Controls the minimum angle and relative lengths of opposing faces
        //- for constructing the upwind-biased stencil
        const scalar minOpposedness_ = 0.5;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        highOrderFit(const highOrderFit&);

        //- Disallow default bitwise assignment
        void operator=(const highOrderFit&);

public:

    //- Runtime type information
    TypeName("highOrderFit");

    // Constructors

        //- Construct using the configuration provided by the \c Istream \c is.
        // @param[in] is an fvSchemes configuration string.
        //               The first word identifies the \c faceFlux
        //               \c surfaceScalarField.
        highOrderFit(const fvMesh& mesh, Istream& is)
        :
            upwind<Type>
            (
                mesh,
                mesh.lookupObject<surfaceScalarField>(word(is))
            )
        {}

        //- Construct a highOrderFit using the specified \c faceFlux
        // @param[in] faceFlux determines the upwind direction for each face
        highOrderFit
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            upwind<Type>(mesh, faceFlux)
        {}


    // Member Functions
        //- Returns \c true since a high-order correction on upwind is applied.
        virtual bool corrected() const
        {
            return true;
        }

        //- Explicit correction to the surface flux interpolate
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            const fvMesh& mesh = this->mesh();

			const extendedUpwindCellToFaceStencil& stencils = upwindCPCCellToFaceStencilObject::New
            (
                mesh,
                false, // not pure upwind
                minOpposedness_,
                false // not a correction on linear
            );

            scalarListList ownerWeights(mesh.nFaces());
            for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
            {
                const labelList& stencil = stencils.ownStencil()[faceI];
                scalarList& weights = ownerWeights[faceI];
                weights.setSize(stencil.size());
                weights[0] = -0.5;
                weights[1] = 0.5;
                for (label i = 2; i < weights.size(); i++)
                {
                    weights[i] = 0.0;
                }
            }

            scalarListList neighbourWeights(mesh.nFaces());
            for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
            {
                const labelList& stencil = stencils.neiStencil()[faceI];
                scalarList& weights = neighbourWeights[faceI];
                weights.setSize(stencil.size());
                weights[0] = -0.5;
                weights[1] = 0.5;
                for (label i = 2; i < weights.size(); i++)
                {
                    weights[i] = 0.0;
                }
            }

            return stencils.weightedSum
            (
                upwind<Type>::faceFlux_,
                vf,
                ownerWeights,
                neighbourWeights
            );
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
