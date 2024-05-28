fvScalarMatrix alphaEqn
(
    fvm::ddt(alpharain[phase_no])
    + fvm::div(phirain[phase_no],alpharain[phase_no],"div(phirain,alpharain)")
   // + fvm::Sp(fvc::div(phirain[phase_no]),alpharain[phase_no]) //added for bounded div scheme // why add?
    ==
    fvOptions(alpharain[phase_no]) // add fvOptions

    // consider scalartransportfoam equations without the laplacian since no turbulence
    // fvm::ddt(T) + fvm::div(phi, T) - fvm::laplacian(DT, T) == fvOptions(T)

);

alphaEqn.relax();

fvOptions.constrain(alphaEqn); // add fvOptions

alphaEqn.solve();

fvOptions.correct(alpharain[phase_no]); // add fvOptions
