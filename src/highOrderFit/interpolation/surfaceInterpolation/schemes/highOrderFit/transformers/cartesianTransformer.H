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
    Foam::highOrderFit::cartesianTransformer

Description
    Translate and rotate stencils in Cartesian space.
    \sa Foam::cartesianCS

SourceFiles
    cartesianTransformer.C

\*---------------------------------------------------------------------------*/

#ifndef cartesianTransformer_H
#define cartesianTransformer_H

#include "stencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace highOrderFit
{

/*---------------------------------------------------------------------------*\
                         Class cartesianTransformer Declaration
\*---------------------------------------------------------------------------*/

class cartesianTransformer
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        cartesianTransformer(const cartesianTransformer&);

        //- Disallow default bitwise assignment
        void operator=(const cartesianTransformer&);


public:

    // Constructors

        //- Construct a default instance.
        cartesianTransformer();


    //- Destructor
    ~cartesianTransformer();


    // Member Functions

        //- Transform a stencil such that the target face centre is the origin and the primary direction is parallel to the target face normal.
        void transform(stencil&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
