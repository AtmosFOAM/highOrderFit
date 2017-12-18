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
    Foam::highOrderFit::faceCache

Description

SourceFiles
    faceCache.C

\*---------------------------------------------------------------------------*/

#ifndef faceCache_H
#define faceCache_H

#include <functional>
#include "face.H"
#include "tet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace highOrderFit
{

/*---------------------------------------------------------------------------*\
                         Class faceCache Declaration
\*---------------------------------------------------------------------------*/

class faceCache
{
    // Private data

        List<face> faces_;
        List<List<tet>> tetsList_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        faceCache(const faceCache&);

        //- Disallow default bitwise assignment
        void operator=(const faceCache&);

        bool approximatelyEqual(const face& a, const face& b) const;

public:

    typedef std::function<void(List<tet>& tets)> decomposer;

    // Constructors

        faceCache();


    //- Destructor
    ~faceCache();


    // Member Functions

        List<tet>& calculateIfNecessary(const face&, const decomposer);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace highOrderFit

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
