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

Application
    diagnoseHighOrderFitInstabilities

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "highOrderFitScheme.H"
#include "weightsFieldPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote("Obtain diagnostics on potentially unstable flux approxations using the highOrderFit surface interpolation scheme");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    ITstream& specification
        = mesh.interpolationScheme("diagnoseHighOrderFit");
    const word schemeName(specification);
    if (schemeName != "highOrderFit")
    {
        FatalErrorInFunction
            << "diagnoseHighOrderFit interpolationScheme must be highOrderFit"
            << abort(FatalError);

    }

    const tmp<surfaceInterpolationScheme<scalar>> scheme = 
        surfaceInterpolationScheme<scalar>::New
    (
        mesh,
        phi,
        mesh.interpolationScheme("diagnoseHighOrderFit")
    );

	const highOrderFitScheme<scalar>& highOrderFit =
        dynamic_cast<const highOrderFitScheme<scalar>&>(scheme());

    const highOrderFit::weightsFieldPair& weights = highOrderFit.weights();

    label instabilities = 0;
    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        const scalarList& w = weights.owner()[facei];
        if (w[0] > 0.0)
        {
            Info << "facei=" << facei << ",owner,"
                 << "upwind_weight=" << w[0] + 1.0
                 << ",Cf=" << mesh.Cf()[facei] << endl;
            instabilities++;
        }
        if (w[1] > 0.5 + 1e5)
        {
            Info << "facei=" << facei << ",owner,"
                 << "downwind_weight=" << w[1]
                 << ",Cf=" << mesh.Cf()[facei] << endl;
            instabilities++;
        }
    }

    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        const scalarList& w = weights.neighbour()[facei];
        if (w[0] > 0.0)
        {
            Info << "facei=" << facei << ",neighbour,"
                 << "upwind_weight=" << w[0] + 1.0
                 << ",Cf=" << mesh.Cf()[facei] << endl;
            instabilities++;
        }
        if (w[1] > 0.5 + 1e5)
        {
            Info << "facei=" << facei << ",neighbour,"
                 << "downwind_weight=" << w[1]
                 << ",Cf=" << mesh.Cf()[facei] << endl;
            instabilities++;
        }
    }

    Info << nl << instabilities << " potential instabilities found." << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
