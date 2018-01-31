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
    diagnoseHighOrderFitFace

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "highOrderFitScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote("Obtain diagnostics on the highOrderFit surface interpolation for a given face and upwind direction");
    argList::validArgs.append("facei");
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

	const label facei = std::stoi(args.arg(1));

	const autoPtr<highOrderFit::Diagnostic<scalar>> diagnostic =
        highOrderFit.diagnose(facei, T);

	Info << diagnostic() << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
