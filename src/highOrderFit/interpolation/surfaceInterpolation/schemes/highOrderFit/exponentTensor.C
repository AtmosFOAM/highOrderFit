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

#include "exponentTensor.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::highOrderFit::exponentTensor::exponentTensor()
{}


Foam::highOrderFit::exponentTensor::exponentTensor
(
            const Foam::label xx,
            const Foam::label xy,
            const Foam::label xz,
            const Foam::label yx,
            const Foam::label yy,
            const Foam::label yz,
            const Foam::label zx,
            const Foam::label zy,
            const Foam::label zz
)
:
    labelTensor(xx, xy, xz, yx, yy, yz, zx, zy, zz)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::highOrderFit::exponentTensor::~exponentTensor()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::highOrderFit::exponentTensor::productOfExponentials
(
    const Foam::tensor A
) const
{
    return pow(A.xx(), xx())
         * pow(A.xy(), xy())
         * pow(A.xz(), xz())
         * pow(A.yx(), yx())
         * pow(A.yy(), yy())
         * pow(A.yz(), yz())
         * pow(A.zx(), zx())
         * pow(A.zy(), zy())
         * pow(A.zz(), zz());
}

Foam::scalar Foam::highOrderFit::exponentTensor::factorialRatio() const
{
    return
    (
        scalar(factorial(xx() + yx() + zx()))
      * scalar(factorial(xy() + yy() + zy()))
      * scalar(factorial(xz() + yz() + zz()))
    )
    /
    ( 
        scalar(factorial(xx()))
      * scalar(factorial(xy()))
      * scalar(factorial(xz()))
      * scalar(factorial(yx()))
      * scalar(factorial(yy()))
      * scalar(factorial(yz()))
      * scalar(factorial(zx()))
      * scalar(factorial(zy()))
      * scalar(factorial(zz()))
    );
}

// ************************************************************************* //
