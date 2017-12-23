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
    Foam::highOrderFit::sphericalTransformer

Description
    Translate and rotate stencils on spherical meshes.
    \sa Foam::cartesianCS

SourceFiles
    sphericalTransformer.C

\*---------------------------------------------------------------------------*/

#ifndef sphericalTransformer_H
#define sphericalTransformer_H

#include "stencil.H"
#include "transformer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace highOrderFit
{

/*---------------------------------------------------------------------------*\
                         Class sphericalTransformer Declaration
\*---------------------------------------------------------------------------*/

class sphericalTransformer
:
    public transformer
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        sphericalTransformer(const sphericalTransformer&);

        //- Disallow default bitwise assignment
        void operator=(const sphericalTransformer&);


public:

    //- Runtime type information
    TypeName("spherical");


    // Constructors

        //- Construct a default instance.
        sphericalTransformer();


    //- Destructor
    ~sphericalTransformer();


    // Member Functions

        //- Transform a stencil such that the target face centre is the origin and the primary direction is parallel to the target face normal.
        //  The stencil z-direction is aligned with the radial direction given
        //  by the target face centre vector.
        void transform(stencil&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
