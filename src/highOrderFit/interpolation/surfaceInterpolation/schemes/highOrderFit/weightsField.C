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

#include "weightsField.H"

#include "stencil.H"
#include "uniformMultipliers.H"
#include "weights.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highOrderFit::weightsField::weightsField
(
    const Foam::highOrderFit::stencilField& stencils
)
:
scalarListList(stencils.size()),
stencils_(stencils)
{
    for (label facei = 0; facei < stencils.mesh().nInternalFaces(); facei++)
    {
        scalarList& w = (*this)[facei];
        const stencil& stencil = stencils_[facei];
        const label size = stencil.size();

        if (size < 2)
        {
            FatalErrorInFunction
                << "stencil for facei " << facei << " has fewer than two cells"
                << abort(FatalError);
        }

        w.setSize(size);

        const uniformMultipliers multipliers(size);

        const weights weights;
        weights.calculate(w, stencil, multipliers);
        w[0] -= 1;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highOrderFit::weightsField::~weightsField()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::autoPtr<Foam::highOrderFit::weightsDiagnostic>
Foam::highOrderFit::weightsField::diagnose(const Foam::label facei) const
{
    const uniformMultipliers multipliers(stencils_[facei].size());

    const weights weights;
    const autoPtr<scalarRectangularMatrix> B = weights.initialiseMatrix
    (
        stencils_[facei],
        multipliers
    );

    return autoPtr<weightsDiagnostic>
    (
        new weightsDiagnostic
        (
            multipliers,
            B
        )
    );
}

// ************************************************************************* //
