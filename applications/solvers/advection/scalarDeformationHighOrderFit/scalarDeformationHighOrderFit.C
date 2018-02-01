/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    scalarDeformationWithGhosts

Description
    Solves a transport equation for a passive scalar using explicit
    time-stepping for the evolving velocity field on a plane
    for deformational flow using ghost cells to handle cyclic boundary conditions

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvGhostMesh.H"
#include "deformationalFlow.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()

    // Create the ghost mesh
    fvGhostMesh ghostMesh
    (
        IOobject
        (
            "ghostMesh", runTime.timeName(), runTime, IOobject::MUST_READ
        ),
        mesh
    );

    // Create the class for the deformational flow
    deformationalFlow defFlow
    (
        IOdictionary
        (
            IOobject
            (
                "deformationalAdvectionDict", "system", runTime,
                IOobject::MUST_READ
            )
        )
    );

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;

        defFlow.update(phi, U, Uf);
        defFlow.update(phiGhost);
        #include "CourantNo.H"

		k1 = -fvc::div(phiGhost, TGhost.oldTime(), "div(phiGhost,TGhost)");
        k1 = ghostMesh.mapToGhost(ghostMesh.mapFromGhost(k1));
		k2 = -fvc::div(phiGhost, TGhost.oldTime() + dt/2 * k1, "div(phiGhost,TGhost)");
        k2 = ghostMesh.mapToGhost(ghostMesh.mapFromGhost(k2));
		k3 = -fvc::div(phiGhost, TGhost.oldTime() + dt/2 * k2, "div(phiGhost,TGhost)");
        k3 = ghostMesh.mapToGhost(ghostMesh.mapFromGhost(k3));
		k4 = -fvc::div(phiGhost, TGhost.oldTime() + dt * k3, "div(phiGhost,TGhost)");

        T = T.oldTime() + dt/6 * 
        (
              ghostMesh.mapFromGhost(k1)
            + 2*ghostMesh.mapFromGhost(k2)
            + 2*ghostMesh.mapFromGhost(k3)
            + ghostMesh.mapFromGhost(k4)
        );
        TGhost = ghostMesh.mapToGhost(T);

        Info << " T goes from " << min(T.internalField()) << " to "
             << max(T.internalField()) << nl << endl;
        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
