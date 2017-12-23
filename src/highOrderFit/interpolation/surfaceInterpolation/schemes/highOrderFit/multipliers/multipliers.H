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
    Foam::highOrderFit::multipliers

Description
    A list of scalars used to weight each cell in the least squares fit.

    \sa Foam::highOrderFit::weights::calculate

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef multipliers_H
#define multipliers_H

#include "runTimeSelectionTables.H"
#include "scalarList.H"
#include "stencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace highOrderFit
{

/*---------------------------------------------------------------------------*\
                         Class multipliers Declaration
\*---------------------------------------------------------------------------*/

class multipliers
{

public:

    //- Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            multipliers,
            word,
            (),
            ()
        );


    // Selectors

        //- Select multipliers instance by type name.
        static autoPtr<multipliers> New(const word& type);


    //- Destructor
    virtual ~multipliers() {};


    // Member Functions

        //- Calculate the list of multipliers using the stencil geometry.
        //  @param [in] stencil
        //  @param [out] m
        virtual void calculate
        (
            const stencil& stencil, scalarList& m
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
