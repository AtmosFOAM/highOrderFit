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
    Foam::weights

Description

SourceFiles
    weights.C

\*---------------------------------------------------------------------------*/

#ifndef weights_H
#define weights_H

#include "scalarList.H"
#include "stencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace highOrderFit
{

/*---------------------------------------------------------------------------*\
                         Class weights Declaration
\*---------------------------------------------------------------------------*/

class weights
{
    // Private data


    // Private Member Functions

        //- Disallow default bitwise copy construct
        weights(const weights&);

        //- Disallow default bitwise assignment
        void operator=(const weights&);


public:

    // Constructors

        weights();


    //- Destructor
    ~weights();


    // Member Functions
    
        //- Calculate weights for a given stencil
        //  @param [out] weights on entry, a list of weights of the correct
        //                       size with uninitialised values;
        //                       on exit, the weights will be initialised
        void calculate
        (
            scalarList& weights,
            const stencil& stencil,
            const scalarList& multipliers
        ) const;
 

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
