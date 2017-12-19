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

#include "order.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::highOrderFit::order::calculateRowCandidates
(
    const Foam::label target, 
    Foam::List<Foam::labelVector>& candidates
) const
{
    for (label a = 0; a <= target; a++)
    {
        for (label b = 0; b <= target; b++)
        {
            for (label c = 0; c <= target; c++)
            {
                if (a + b + c == target)
                {
                    candidates.append(labelVector(a, b, c));
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highOrderFit::order::order
(
    const Foam::label x,
    const Foam::label y,
    const Foam::label z
)
:
Foam::labelVector(x, y, z)
{}


Foam::highOrderFit::order::order(const labelVector& v, faceCache& cache)
:
Foam::labelVector(v),
cache_(&cache)        
{}

Foam::highOrderFit::order::order()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highOrderFit::order::~order()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::highOrderFit::order::factorialRatio
(
    const label dimensions
) const
{
    return
        scalar(factorial(x()) * factorial(y()) * factorial(z()))
        /
        scalar(factorial(x() + y() + z() + dimensions));
}


void Foam::highOrderFit::order::calculateExponentTensors
(
    Foam::List<Foam::highOrderFit::exponentTensor>& K
) const
{
    List<labelVector> row0Candidates, row1Candidates, row2Candidates;
    calculateRowCandidates(x(), row0Candidates);
    calculateRowCandidates(y(), row1Candidates);
    calculateRowCandidates(z(), row2Candidates);

    forAll(row0Candidates, i)
    {
        forAll(row1Candidates, j)
        {
            forAll(row2Candidates, k)
            {
                K.append
                (
                    exponentTensor
                    (
                        row0Candidates[i],
                        row1Candidates[j],
                        row2Candidates[k]
                    )
                );
            }
        }
    }
}


bool Foam::highOrderFit::order::cacheEnabled() const
{
    return cache_ != nullptr;
}


/*faceCache& Foam::highOrderFit::order::cache() const
{
    return *cache_;
}*/

// ************************************************************************* //
