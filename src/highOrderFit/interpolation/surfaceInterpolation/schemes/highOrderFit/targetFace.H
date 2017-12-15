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
    Foam::highOrderFit::targetFace

Description

SourceFiles
    targetFace.C

\*---------------------------------------------------------------------------*/

#ifndef targetFace_H
#define targetFace_H

#include "face.H"
#include "point.H"
#include "stencil.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace highOrderFit
{

// Forward declaration of classes
class stencil;

/*---------------------------------------------------------------------------*\
                         Class targetFace Declaration
\*---------------------------------------------------------------------------*/

class targetFace
{
    // Private data

        const point Cf_;
        const vector Sf_;

    // Private Member Functions

        // Disallow default bitwise copy construct
        //targetFace(const targetFace&);

        // Disallow default bitwise assignment
        //void operator=(const targetFace&);


public:


    // Constructors
    
        targetFace(const fvMesh& mesh, const label facei);

        targetFace(const point Cf, const vector Sf, const face f);


    //- Destructor
    ~targetFace();


    // Member Functions

        void translate(stencil& s) const;

        void rotate(stencil& s) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
