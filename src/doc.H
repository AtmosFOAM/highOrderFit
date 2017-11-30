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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

\mainpage A high-order surfaceInterpolationScheme

The \ref Foam::highOrderFitScheme "highOrderFitScheme" is an upwind-biased surfaceInterpolationScheme for solving the transport equation.
The scheme is high-order accurate on arbitrary polyhedral meshes and uses the high-order moments method of <a href="https://doi.org/10.2140/camcos.2017.12.51">Devendran et al. 2017</a>.
The scheme uses an upwind-biased stencil and a weighted least squares fit as documented by <a href="https://doi.org/10.1016/j.jcp.2017.04.061">Shaw et al. 2017</a>.

\namespace Foam::highOrderFit
\brief Contains classes belonging to the highOrderFitScheme

\*---------------------------------------------------------------------------*/
