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

#include "mesh.H"

#include "surfaceFields.H"
#include "OStringStream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Test::mesh::mesh(const Foam::fvMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Test::mesh::indexOfCellWithCentreAt
(
    const Foam::point& C,
    const Foam::scalar epsilon
) const
{
    forAll(mesh_.C(), celli)
    {
        if (Foam::magSqr(mesh_.C()[celli] - C) < epsilon) return celli;
    }

	Foam::OStringStream os;
	os << "no cell with centre at " << C;
	throw std::domain_error(os.str());
}


Foam::label Test::mesh::indexOfFaceWithCentreAt
(
    const Foam::point& Cf,
    const Foam::scalar epsilon
) const
{
    forAll(mesh_.Cf(), facei)
    {
        if (Foam::magSqr(mesh_.Cf()[facei] - Cf) < epsilon) return facei;
    }

	Foam::OStringStream os;
	os << "no face with centre at " << Cf;
	throw std::domain_error(os.str());
}

// ************************************************************************* //
