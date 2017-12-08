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
    Foam::highOrderFit::weightsField

Description
    A list of surface interpolation weights indexed by face.

SourceFiles
    weightsField.C

\*---------------------------------------------------------------------------*/

#ifndef weightsField_H
#define weightsField_H

#include "autoPtr.H"
#include "MeshObject.H"
#include "order.H"
#include "stencilField.H"
#include "scalarList.H"
#include "weights.H"
#include "weightsDiagnostic.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace highOrderFit
{

/*---------------------------------------------------------------------------*\
                         Class weightsField Declaration
\*---------------------------------------------------------------------------*/

class weightsField
:
    public scalarListList
{
    // Private data
        const stencilField stencils_;
        const weights weights_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        weightsField(const weightsField&);

        //- Disallow default bitwise assignment
        void operator=(const weightsField&);


public:

    // Constructors

        //- Construct from the specified stencil field.
        weightsField
        (
            const fvMesh& mesh,
            const stencilField& stencils,
            const List<order>& moments
        );


    //- Destructor
    ~weightsField();


    // Member Functions

        const autoPtr<weightsDiagnostic> diagnose(const label facei) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
