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
#include "SVD.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::List<Foam::highOrderFit::order>&
Foam::highOrderFit::weights::momentsFor
(
    const Foam::highOrderFit::stencil& stencil
) const
{
    if (stencil.size() >= 12)
    {
        return moments_;
    }
    else
    {
        return linearMoments_;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highOrderFit::weights::weights
(
    const Foam::List<Foam::highOrderFit::order>& moments
)
:
moments_(moments),
linearMoments_({order(0, 0, 0), order(1, 0, 0), order(0, 0, 1)})
{}


Foam::highOrderFit::weights::weights
(
    const Foam::List<Foam::highOrderFit::order>& moments,
    const Foam::List<Foam::highOrderFit::order>& linearMoments
)
:
moments_(moments),
linearMoments_(linearMoments)
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
    const autoPtr<scalarRectangularMatrix> B =
        createMatrix(stencil, multipliers);
    const SVD svd(B());
    const scalarRectangularMatrix& Binv = svd.VSinvUt();
    const List<order>& moments = momentsFor(stencil);

    for (label col = 0; col < Binv.n(); col++)
    {
        weights[col] = 0.0;
        forAll(moments, row)
        {
            weights[col] += Binv(row, col) *
                stencil.targetFaceMomentAverage(moments[row]);
        }
        weights[col] *= multipliers[col];
    }
}


const Foam::autoPtr<Foam::scalarRectangularMatrix>
Foam::highOrderFit::weights::createMatrix
(
    const Foam::highOrderFit::stencil& stencil,
    const Foam::scalarList& multipliers
) const
{
    const List<order>& moments = momentsFor(stencil);

    scalarRectangularMatrix* B = new scalarRectangularMatrix
    (
        stencil.size(),
        moments.size()
    );

    for (label row = 0; row < B->m(); row++)
    {
        const scalar volume = stencil[row].moment(order(0, 0, 0));

        forAll(moments, col)
        {
            if (moments[col] == order(0, 0, 0))
            {
                (*B)(row, col) = multipliers[row];
            }
            else
            {
                (*B)(row, col) =
                    stencil[row].moment(moments[col]) * multipliers[row];
                (*B)(row, col) /= volume;
            }
        }
    }

    return autoPtr<scalarRectangularMatrix>(B);
}

// ************************************************************************* //
