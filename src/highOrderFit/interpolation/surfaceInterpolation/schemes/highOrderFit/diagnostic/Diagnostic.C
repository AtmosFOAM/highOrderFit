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

#include "Diagnostic.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::highOrderFit::Diagnostic<Type>::Diagnostic
(
    const Foam::label facei,
    const Foam::GeometricField<Type, fvPatchField, volMesh>& field,
    const Foam::surfaceScalarField& faceFlux,
    const Foam::fvMesh& mesh,
    const Foam::extendedUpwindCellToFaceStencil& stencils,
    const Foam::highOrderFit::weightsFieldPair& weights
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

    weights_ = weights.side(owner_)[facei_];
    weightsDiagnostic_ = weights.side(owner_).diagnose(facei_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::highOrderFit::Diagnostic<Type>::~Diagnostic<Type>()
{}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::highOrderFit::operator<<
(
    Foam::Ostream& os,
    const Foam::highOrderFit::Diagnostic<Type>& d
)
{
    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::highOrderFit::operator<<(Foam::Ostream&, "
        "const Foam::highOrderFit::Diagnostic&)"
    );

    const label facei = d.facei_;

    os << "highOrderFit[facei=" << facei << ", "
       << (d.owner_ ? "owner" : "neighbour")
       << ", ";

    for (label i = 0; i < d.weights_.size(); i++)
    {
        os << (i == 0 ? d.weights_[i] + 1 : d.weights_[i])
           << "*" << d.values_[facei][i]
           << "@" << d.cellCentres_[facei][i]
           << "{m=" << d.weightsDiagnostic_->multipliers()[i] << "}";
        if (i < d.weights_.size() - 1)
        {
            os << " + ";
        }
    }

    os << ", B=" << d.weightsDiagnostic_->B();
       
    os << "]";

    return os;
}


// ************************************************************************* //
