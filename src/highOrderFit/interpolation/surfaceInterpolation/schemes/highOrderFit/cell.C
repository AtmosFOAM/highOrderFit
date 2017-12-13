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

#include "cell.H"
#include "face.H"
#include "labelList.H"
#include "quaternion.H"
#include "tet.H"
#include "transform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highOrderFit::cell::cell
(
    const Foam::fvMesh& mesh,
    const Foam::label celli
)
:
List<face>(mesh.cells()[celli].size()),
C_(mesh.C()[celli])
{
    const labelList& faces = mesh.cells()[celli];
    forAll(faces, i)
    {
        highOrderFit::face& stencilCellFace = (*this)[i];
        const Foam::face& face = mesh.faces()[faces[i]];
        
        bool owner = (mesh.owner()[faces[i]] == celli);

        stencilCellFace.setSize(face.size());
        forAll(face, pointi)
        {
            if (owner)
            {
                stencilCellFace[pointi] = mesh.points()[face[pointi]];
            }
            else
            {
                stencilCellFace[face.size() - pointi - 1] =
                    mesh.points()[face[pointi]];
            }
        }
    }
}


Foam::highOrderFit::cell::cell
(
    const Foam::List<Foam::highOrderFit::face>& faces,
    const Foam::point centre
)
:
Foam::List<Foam::highOrderFit::face>(faces),
C_(centre)
{}


Foam::highOrderFit::cell::cell(Istream& is)
:
C_(0, 0, 0)
{
    // Check state of Istream
    is.check("Foam::highOrderFit::cell::cell(Foam::Istream&)");

    operator>>(is, *this);
}


Foam::highOrderFit::cell::cell()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highOrderFit::cell::~cell()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::highOrderFit::cell::translate(const Foam::vector x)
{
    C_ += x;

    forAll((*this), facei)
    {
        forAll((*this)[facei], pointi)
        {
            (*this)[facei][pointi] += x;
        }
    }
}


void Foam::highOrderFit::cell::rotate
(
    const Foam::vector from,
    const Foam::vector to
)
{
    const quaternion q(rotationTensor(from, to));

    C_ = q.transform(C_);

    forAll((*this), facei)
    {
        forAll((*this)[facei], pointi)
        {
            (*this)[facei][pointi] = q.transform((*this)[facei][pointi]);
        }
    }
}


Foam::scalar Foam::highOrderFit::cell::moment
(
    const Foam::highOrderFit::order& o
) const
{
    scalar moment = 0.0;

    forAll((*this), facei)
    {
        List<tet> faceTets;
        (*this)[facei].decompose(faceTets);

        forAll(faceTets, teti)
        {
            moment += faceTets[teti].volumeMoment(o);
        }
    }

    const scalar epsilon = 1e-6;

    if (o == order(0, 0, 0))
    {
        return moment;
    }
    else if (o == order(1, 0, 0))
    {
        return C_.x();
    }
    else if (o == order(2, 0, 0))
    {
        const scalar x = C_.x();
        if (fabs(x) >= 2.5 - epsilon)
        {
            return 19.0/3.0;
        }
        else if (fabs(x) >= 1.5 - epsilon)
        {
            return 7.0/3.0;
        }
        else
        {
            return 1.0/3.0;
        }
    }
    else
    {
        const scalar x = C_.x();
        if (x <= -2.5 + epsilon)
        {
            return -65.0/4.0;
        }
        else if (x <= -1.5 + epsilon)
        {
            return -15.0/4.0;
        }
        else if (x <= -0.5 + epsilon)
        {
            return -1.0/4.0;
        }
        else if (x <= 0.5 + epsilon)
        {
            return 1.0/4.0;
        }
        else if (x <= 1.5 + epsilon)
        {
            return 15.0/4.0;
        }
        else
        {
            return 65.0/4.0;
        }
    }
}


Foam::point Foam::highOrderFit::cell::centre() const
{
    return C_;
}


// ************************************************************************* //
