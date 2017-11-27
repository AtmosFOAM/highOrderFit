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
    Foam::highOrderFitDiagnostic

Description
    A diagnostic object for a particular face and face flux which can be
    inserted into an Ostream.

SourceFiles
    highOrderFitDiagnostic.C
    highOrderFitDiagnosticIO.C

\*---------------------------------------------------------------------------*/

#ifndef highOrderFitDiagnostic_H
#define highOrderFitDiagnostic_H

#include "label.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes                                               
 class Ostream;

// Forward declaration of friend functions and operators
class highOrderFitDiagnostic;
Ostream& operator<<(Ostream&, const highOrderFitDiagnostic&);


/*---------------------------------------------------------------------------*\
                         Class highOrderFitDiagnostic Declaration
\*---------------------------------------------------------------------------*/

class highOrderFitDiagnostic
{
    // Private data

        const label facei_;
        const bool owner_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        highOrderFitDiagnostic(const highOrderFitDiagnostic&);

        //- Disallow default bitwise assignment
        void operator=(const highOrderFitDiagnostic&);


public:

    // Constructors

        highOrderFitDiagnostic
        (
            const label facei,
            const bool owner
        );


    //- Destructor
    ~highOrderFitDiagnostic();


    // Member Functions

        // Access
        // Check
        // Edit
        // Write


    // IOstream Operators

        friend Ostream& operator<<
        (
            Ostream& os,
            const highOrderFitDiagnostic& diagnostic
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
