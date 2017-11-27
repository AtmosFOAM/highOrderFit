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
#include "highOrderFitStencilField.H"
#include "highOrderFitWeightsField.H"
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

        const extendedUpwindCellToFaceStencil& stencils_;
        const highOrderFitStencilField ownerStencils_;
        const highOrderFitWeightsField ownerWeights_;
        const highOrderFitStencilField neighbourStencils_;
        const highOrderFitWeightsField neighbourWeights_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        highOrderFit(const highOrderFit&);

        //- Disallow default bitwise assignment
        void operator=(const highOrderFit&);

        const extendedUpwindCellToFaceStencil& initialiseStencils
        (
            const fvMesh& mesh
        )
        {
            return upwindCPCCellToFaceStencilObject::New
            (
                mesh,
                false, // not pure upwind
                minOpposedness_,
                false // not a correction on linear
            );
        }

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
            ),
            stencils_(initialiseStencils(mesh)),
            ownerStencils_
            (
                stencils_.ownStencil(),
                mesh
            ),
            ownerWeights_(ownerStencils_),
            neighbourStencils_
            (
                stencils_.neiStencil(),
                mesh
            ),
            neighbourWeights_(neighbourStencils_)
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
            upwind<Type>(mesh, faceFlux),
            stencils_(initialiseStencils(mesh)),
            ownerStencils_
            (
                stencils_.ownStencil(),
                mesh
            ),
            ownerWeights_(ownerStencils_),
            neighbourStencils_
            (
                stencils_.neiStencil(),
                mesh
            ),
            neighbourWeights_(neighbourStencils_)
        {}


    // Member Functions
    
        const highOrderFitWeightsField& ownerWeights() const
        {
            return ownerWeights_;
        }

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
            return stencils_.weightedSum
            (
                upwind<Type>::faceFlux_,
                vf,
                ownerWeights_,
                neighbourWeights_
            );
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
