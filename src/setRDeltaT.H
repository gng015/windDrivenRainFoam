{

    volScalarField& rDeltaT = trDeltaT.ref();

    const dictionary& pimpleDict = pimple.dict();

    scalar maxCo
    (
        pimpleDict.getOrDefault<scalar>("maxCo", 0.9)
    );

    scalar maxDeltaT
    (
        pimpleDict.getOrDefault<scalar>("maxDeltaT", 10)
    );

    scalar rDeltaTDampingCoeff
    (
        pimpleDict.getOrDefault<scalar>("rDeltaTDampingCoeff", 1.0)
    );

    scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.getOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.1)
    );

    volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.ref() = max
    (
        1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum(phirain[current_phase_no])()()
    /((2*maxCo)*mesh.V())
    );

    // Update the boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    Info<< "Smoothed flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;


    if
    (
        rDeltaTDampingCoeff < 1.0
    && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        rDeltaT =
            rDeltaT0
        *max(rDeltaT/rDeltaT0, scalar(1) - rDeltaTDampingCoeff);

        Info<< "Damped flow time scale min/max = "
            << gMin(1/rDeltaT.primitiveField())
            << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
    }

}


