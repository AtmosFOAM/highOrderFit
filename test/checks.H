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

Class
    Foam::checks

Description
    A collection of assertions for use in the Catch test suite

SourceFiles
    checks.C
    checksTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef checks_H
#define checks_H

#include "catch.hpp"
#include "GeometricField.H"
#include "scalar.H"
#include "scalarList.H"
#include "Vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Test
{

    static Foam::scalar margin = 1e-9;
    static Approx approx = Approx::custom().margin(margin);

    void checkEqual
    (
        const Foam::scalarList& actual,
        const Foam::scalarList& expected,
        Approx approx = Test::approx
    );

    template<class Type>
    void checkEqual
    (
        const Foam::List<Type>& actual,
        const Foam::List<Type>& expected
    );

    template<class Type>
    void checkEqual
    (
        const Foam::Tensor<Type>& actual,
        const Foam::Tensor<Type>& expected
    );

    template<class Type>
    void checkEqual
    (
        const Foam::List<Type>& actual,
        const Type expected
    );

    void checkEqual
    (
        const Foam::scalarList& actual,
        const Foam::scalar expected,
        Approx approx = Test::approx
    );

    void checkEqual
    (
        const Foam::vector actual,
        const Foam::vector expected,
        Approx approx = Test::approx
    );

    void checkEqual
    (
        const Foam::scalarRectangularMatrix& actual,
        const Foam::scalarRectangularMatrix& expected,
        Approx approx = Test::approx
    );

    template<class Type>
    Foam::label countMatches
    (
        const Foam::List<Type>& list,
        const Type item
    );

    template<class Type>
    Foam::label countMatches
    (
        const Foam::Vector<Type>& vector,
        const Type item
    );

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "checksTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
