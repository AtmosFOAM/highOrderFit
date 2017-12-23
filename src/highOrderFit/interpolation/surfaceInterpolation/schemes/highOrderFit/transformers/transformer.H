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
    Foam::highOrderFit::transformer

Description

SourceFiles
    transformer.C

\*---------------------------------------------------------------------------*/

#ifndef transformer_H
#define transformer_H

#include "runTimeSelectionTables.H"
#include "stencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace highOrderFit
{


/*---------------------------------------------------------------------------*\
                         Class transformer Declaration
\*---------------------------------------------------------------------------*/

class transformer
{

    // Private Member Functions

        // Disallow default bitwise copy construct
        //transformer(const transformer&);

        // Disallow default bitwise assignment
        //void operator=(const transformer&);


public:

	//- Declare run-time constructor selection table                            

        declareRunTimeSelectionTable
        (
            autoPtr,
            transformer,
            word,
            (),
            ()
        );


    // Selectors

        //- Select transformer instance by type name.
        static autoPtr<transformer> New(const word& type);


    //- Destructor
    virtual ~transformer() {};


    // Member Functions

		//- Translate and rotate the stencil from global to local coordinates.
		virtual void transform(stencil&) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
