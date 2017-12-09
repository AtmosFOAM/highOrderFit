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
    Foam::highOrderFit::weightsDiagnostic

Description

SourceFiles
    weightsDiagnostic.C

\*---------------------------------------------------------------------------*/

#ifndef weightsDiagnostic_H
#define weightsDiagnostic_H

#include "scalarList.H"
#include "scalarMatrices.H"
#include "stencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace highOrderFit
{

/*---------------------------------------------------------------------------*\
                         Class weightsDiagnostic Declaration
\*---------------------------------------------------------------------------*/

class weightsDiagnostic
{
    // Private data

        const scalarList multipliers_;
        const scalarRectangularMatrix B_;
        const stencil& stencil_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        weightsDiagnostic(const weightsDiagnostic&);

        //- Disallow default bitwise assignment
        void operator=(const weightsDiagnostic&);


public:

    // Constructors

        weightsDiagnostic
        (
            const scalarList& multipliers,
            const scalarRectangularMatrix& B,
            const stencil& stencil
        );


    //- Destructor
    ~weightsDiagnostic();


    // Member Functions

        const scalarList& multipliers() const
        {
            return multipliers_;
        }


        const scalarRectangularMatrix& B() const
        {
            return B_;
        }

        const stencil& getStencil() const
        {
            return stencil_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
