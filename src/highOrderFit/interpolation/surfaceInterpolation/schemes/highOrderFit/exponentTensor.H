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
    Foam::highOrderFit::exponentTensor

Description
    Represents a 3&times;3 tensor of non-negative integer exponents.

    exponentTensors belong to the set K as defined by 
    <a href="https://doi.org/10.1016/S0031-3203(03)00127-4">
    Tuzikov et al. 2003</a>.

SourceFiles
    exponentTensor.C

\*---------------------------------------------------------------------------*/

#ifndef exponentTensor_H
#define exponentTensor_H

#include "labelTensor.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace highOrderFit
{

/*---------------------------------------------------------------------------*\
                         Class exponentTensor Declaration
\*---------------------------------------------------------------------------*/

class exponentTensor
:
    public labelTensor
{

    // Private Member Functions

        // Disallow default bitwise copy construct
        //exponentTensor(const exponentTensor&);

        // Disallow default bitwise assignment
        //void operator=(const exponentTensor&);


public:

    // Constructors

        //- Construct an null instance.
        //  Necessary only for initialising collections of exponentTensors.
        exponentTensor();

        exponentTensor(const vector& x, const vector& y, const vector& z);

        exponentTensor
        (
            const label xx,
            const label xy,
            const label xz,
            const label yx,
            const label yy,
            const label yz,
            const label zx,
            const label zy,
            const label zz
        );

        //- Construct by reading from an input stream.
        exponentTensor(Istream&);

    //- Destructor
    ~exponentTensor();


    // Member Functions

        //- Calculate the product of A<sub>ij</sub><sup><tt>this</tt><sub>ij</sub></sup> over all elements.
        scalar productOfExponentials(const tensor A) const;

        //- Calculate the ratio of factorials used in calculating a tet moment.
        //  The ratio is the product of <tt>(xj + yj + jh)!</tt> over all
        //  columns <t>j</tt>, divided by the product of <tt>(ij)!</tt>
        //  over all elements.
        //  \sa Foam::highOrderFit::tet::volumeMoment
        //  \sa Foam::highOrderFit::tet::surfaceMoment
        scalar factorialRatio() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
