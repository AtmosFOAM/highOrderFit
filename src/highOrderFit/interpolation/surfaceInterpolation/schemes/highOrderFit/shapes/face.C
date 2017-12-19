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
#include "faceCache.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highOrderFit::face::face()
{}


Foam::highOrderFit::face::face
(
    const Foam::fvMesh& mesh,
    const Foam::label facei
)
:
List<point>(mesh.faces()[facei].size())
{
    const Foam::face& f = mesh.faces()[facei];
    forAll(f, pointi)
    {
        (*this)[pointi] = mesh.points()[f[pointi]];
    }
}


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


void Foam::highOrderFit::face::translate(const Foam::vector x)
{
    forAll((*this), pointi)
    {
        (*this)[pointi] += x;
    }
}


void Foam::highOrderFit::face::rotate
(
    const Foam::vector from,
    const Foam::vector to
)
{
    const quaternion q(rotationTensor(from, to));
    rotate(q);
}


void Foam::highOrderFit::face::rotate(const Foam::quaternion& q)
{
    forAll((*this), pointi)
    {
        (*this)[pointi] = q.transform((*this)[pointi]);
    }
}


Foam::scalar Foam::highOrderFit::face::moment
(
    const Foam::highOrderFit::order& o
) const
{
    scalar moment = 0.0;

    /*if (o.cacheEnabled())
    {

        const faceCache::decomposer func = std::bind(&face::decompose, this, std::placeholders::_1);
        o.cache().calculateIfNecessary(*this, func);
    }*/
    List<tet> faceTets;
    decompose(faceTets);

    forAll(faceTets, teti)
    {
        moment += faceTets[teti].surfaceMoment(o);
    }

    return moment;
}


void Foam::highOrderFit::face::flip()
{
    const label n = size();

    for (label i=1; i < (n+1)/2; ++i)
    {
        Swap(operator[](i), operator[](n-i));
    }
}

// ************************************************************************* //
