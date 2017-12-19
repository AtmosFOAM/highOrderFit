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
    Foam::highOrderFit::order

Description
    The <i>(x, y, z)</i> exponents of a polynomial term.

    For example, <tt>order(2, 0, 1)</tt> are the exponents for the
    polynomial term <i>x&sup2;z</i>.

SourceFiles
    order.C

\*---------------------------------------------------------------------------*/

#ifndef order_H
#define order_H

#include "exponentTensor.H"
#include "Hash.H"
#include "labelVector.H"
#include "labelTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace highOrderFit
{

class faceCache;

/*---------------------------------------------------------------------------*\
                         Class order Declaration
\*---------------------------------------------------------------------------*/

class order
    : public labelVector
{
    // Private data

        faceCache* cache_ = nullptr;
 
    // Private Member Functions

        void calculateRowCandidates
        (
            const label target, 
            List<labelVector>& candidates
        ) const;

public:

    //- Hashing function class.
    template<class HashT=Hash<Foam::highOrderFit::order>>
    class Hash
    {
    public:
        Hash()
        {}

        inline unsigned operator()
        (
            const order&,
            unsigned seed = 0
        ) const;
    };

    // Constructors
    
        //- Construct with specified exponents.
        order(const label x, const label y, const label z);

        //- Construct with exponents specified by the labelVector.
        order(const labelVector&, faceCache& cache);

        //- Construct an instance with zero exponents.
        //  Necessary only for initialising collections of orders.
        order();


    //- Destructor
    ~order();


    // Member Functions

        //- Calculate x!y!z!/(x+y+z+dimensions)!.
        scalar factorialRatio(const label dimensions) const;

        //- Populate the list of \ref exponentTensor "exponentTensors" K.
        //  Each exponentTensor has elements k<sub>ij</sub> such that
        //  0 &le; k<sub>ij</sub> &le; k<sub>i</sub> and 
        //  k<sub>i,1</sub> + k<sub>i,2</sub> + k<sub>i,3</sub> = k<sub>i</sub>
        //  where k<sub>i</sub> is the order element x, y or z.
        //  @param [out] K on entry, an empty list; on exit, a list of 
        //               exponentTensors.
        void calculateExponentTensors(List<exponentTensor>& K) const;

        bool cacheEnabled() const;

        //faceCache& cache() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "orderI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
