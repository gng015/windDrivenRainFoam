{
    volScalarField& rDeltaT = trDeltaT.ref();

    const dictionary& pimpleDict = pimple.dict();

    scalar maxCo
    (
        pimpleDict.getOrDefault<scalar>("maxCo", 0.9)
    );

    scalar maxDeltaT
    (
        pimpleDict.getOrDefault<scalar>("maxDeltaT", GREAT)
    );

    scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.getOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.1)
    );

    surfaceScalarField maxPhi("maxPhi", phi);

    forAll (phases, phase_no)
    {
        maxPhi = max(maxPhi, phirain[phase_no]);
    }

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.ref() = max
    (
        1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum(maxPhi)()()
       /((2*maxCo)*mesh.V())
    );

    // Update the boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
}
