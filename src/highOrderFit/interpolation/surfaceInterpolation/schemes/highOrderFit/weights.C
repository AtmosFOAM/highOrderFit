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

#include "weights.H"
#include "scalarMatrices.H"
#include "SVD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highOrderFit::weights::weights()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highOrderFit::weights::~weights()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::highOrderFit::weights::calculate
(
    Foam::scalarList& weights,
    const Foam::highOrderFit::stencil& stencil,
    const Foam::scalarList& multipliers
) const
{
    scalarRectangularMatrix B(weights.size(), 1);
    for (label row = 0; row < B.m(); row++)
    {
        B(row, 0) = 1.0 * multipliers[row];
    }
    const SVD svd(B);
    const scalarRectangularMatrix& Binv = svd.VSinvUt();

    for (label col = 0; col < Binv.n(); col++)
    {
        weights[col] = Binv(0, col) * multipliers[col];
    }
}


// ************************************************************************* //
