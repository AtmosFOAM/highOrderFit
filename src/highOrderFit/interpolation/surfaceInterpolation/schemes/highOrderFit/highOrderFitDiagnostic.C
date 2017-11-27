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

#include "highOrderFitDiagnostic.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::highOrderFitDiagnostic<Type>::highOrderFitDiagnostic
(
    const Foam::label facei,
    const Foam::GeometricField<Type, fvPatchField, volMesh>& field,
    const Foam::surfaceScalarField& faceFlux,
    const Foam::fvMesh& mesh,
    const Foam::extendedUpwindCellToFaceStencil& stencils,
    const Foam::highOrderFitWeightsField& ownerWeights,
    const Foam::highOrderFitWeightsField& neighbourWeights
)
:
    facei_(facei),
    owner_(faceFlux[facei] >= 0)
{
    stencils.collectData
    (
        owner_ ? stencils.ownMap() : stencils.neiMap(),
        owner_ ? stencils.ownStencil() : stencils.neiStencil(),
        field,
        values_
    );

    stencils.collectData
    (
        owner_ ? stencils.ownMap() : stencils.neiMap(),
        owner_ ? stencils.ownStencil() : stencils.neiStencil(),
        mesh.C(),
        cellCentres_
    );

    weights_ = owner_ ? ownerWeights[facei_] : neighbourWeights[facei_];
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::highOrderFitDiagnostic<Type>::~highOrderFitDiagnostic<Type>()
{}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const highOrderFitDiagnostic<Type>& diagnostic
)
{
    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<(Foam::Ostream&, "
        "const Foam::highOrderFitDiagnostic&)"
    );

    const label facei = diagnostic.facei_;

    os << "highOrderFit[facei=" << facei << ", "
       << (diagnostic.owner_ ? "owner" : "neighbour")
       << ", ";

    for (label i = 0; i < diagnostic.weights_.size(); i++)
    {
        os << (i == 0 ? diagnostic.weights_[i] + 1 : diagnostic.weights_[i])
           << "*" << diagnostic.values_[facei][i]
           << "@" << diagnostic.cellCentres_[facei][i];
        if (i < diagnostic.weights_.size() - 1)
        {
            os << " + ";
        }
    }
       
    os << "]";

    return os;
}


// ************************************************************************* //
