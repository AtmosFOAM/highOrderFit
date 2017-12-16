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
    Represents the target face of a stencil through which a flux is calculated.

    \sa Foam::highOrderFit::stencil

SourceFiles
    targetFace.C

\*---------------------------------------------------------------------------*/

#ifndef targetFace_H
#define targetFace_H

#include "cell.H"
#include "face.H"
#include "fvMesh.H"
#include "point.H"
#include "surfaceFields.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace highOrderFit
{

/*---------------------------------------------------------------------------*\
                         Class targetFace Declaration
\*---------------------------------------------------------------------------*/

class targetFace
    : public face
{
    // Private data

        point Cf_;
        vector unitNormal_;
        scalar area_;

    // Private Member Functions

        // Disallow default bitwise copy construct
        //targetFace(const targetFace&);

        // Disallow default bitwise assignment
        //void operator=(const targetFace&);

        void transform();

public:

    // Constructors

        targetFace();
    
        targetFace(const fvMesh& mesh, const label facei);

        targetFace
        (
            std::initializer_list<point> lst,
            const point Cf,
            const vector Sf
        );


    //- Destructor
    ~targetFace();


    // Member Functions

        //- Translate the stencil such that target face centre is at the origin.
        void translate(List<cell>&) const;

        //- Rotate the stencil such that the primary direction is parallel to the target face area vector.
        void rotate(List<cell>&) const;

        scalar momentAverage(const order&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
