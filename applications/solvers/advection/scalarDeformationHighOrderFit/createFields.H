Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading fields phi (compulsory), U and Uf (optional)\n" << endl;

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Read U if present, otherwise create and write (not used)
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::reconstruct(phi)
);

// Read Uf if present, otherwise create and write (not used)
surfaceVectorField Uf
(
    IOobject
    (
        "Uf",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(fvc::reconstruct(phi))
);

// Create fields on the ghost mesh
volScalarField TGhost
(
    IOobject("TGhost", runTime.timeName(), ghostMesh, IOobject::MUST_READ,
    IOobject::AUTO_WRITE),
    ghostMesh
);
TGhost = ghostMesh.mapToGhost(T);
TGhost.write();
TGhost.oldTime();

surfaceScalarField phiGhost
(
    IOobject("phiGhost", runTime.timeName(), ghostMesh, IOobject::MUST_READ),
    ghostMesh
);

defFlow.update(phi, U, Uf);
defFlow.update(phiGhost);
phi.write();
Uf.write();
U.write();

volScalarField k1
(
    IOobject
    (
        "k1",
        runTime.timeName(),
        ghostMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -fvc::div(phiGhost,TGhost)
);

volScalarField k2
(
    IOobject
    (
        "k2",
        runTime.timeName(),
        ghostMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -fvc::div(phiGhost,TGhost)
);

volScalarField k3
(
    IOobject
    (
        "k3",
        runTime.timeName(),
        ghostMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -fvc::div(phiGhost,TGhost)
);

volScalarField k4
(
    IOobject
    (
        "k4",
        runTime.timeName(),
        ghostMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -fvc::div(phiGhost,TGhost)
);
