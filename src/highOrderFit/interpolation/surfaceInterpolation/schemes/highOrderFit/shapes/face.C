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

#include "face.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highOrderFit::face::face()
{}


Foam::highOrderFit::face::face(std::initializer_list<point> lst)
:
    List<point>(lst)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highOrderFit::face::~face()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::highOrderFit::face::decompose
(
    Foam::List<Foam::highOrderFit::tet>& tets
) const
{
    point centre(0, 0, 0);

    forAll((*this), i)
    {
        centre += (*this)[i];
    }
    centre /= size();

    tets.setSize(size());
    forAll((*this), i)
    {
        tets[i] = tet(centre, (*this)[i], (*this)[(i+1)%size()]);
    }
}


Foam::scalar Foam::highOrderFit::face::moment
(
    const Foam::highOrderFit::order& o
) const
{
    if (o == order(0, 0, 0))
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}

// ************************************************************************* //
