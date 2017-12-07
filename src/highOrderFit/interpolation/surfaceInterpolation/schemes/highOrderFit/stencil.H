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
    Foam::highOrderFit::stencil

Description

SourceFiles
    stencil.C

\*---------------------------------------------------------------------------*/

#ifndef stencil_H
#define stencil_H

#include "cellVertices.H"
#include "label.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace highOrderFit
{

/*---------------------------------------------------------------------------*\
                         Class stencil Declaration
\*---------------------------------------------------------------------------*/

class stencil
:
    public List<cellVertices>
{
    // Private data
    
        label size_;
        scalarList zeroMoment_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        stencil(const stencil&);

        //Disallow default bitwise assignment
        //void operator=(const stencil&);


        void translateVerticesWithOrigin(const point o);

        void rotateVerticesWithPrimaryDirection(const vector x);


public:

    // Constructors

        stencil();

        stencil
        (
            const point targetCf,
            const vector Sf,
            const List<cellVertices>& vertices
        );


    //- Destructor
    ~stencil();


    // Member Functions

        label size() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
