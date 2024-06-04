fvScalarMatrix alphaEqn
(
    fvm::ddt(alpharain[phase_no])
    + fvm::div(phirain[phase_no],alpharain[phase_no],"div(phirain,alpharain)")
   // + fvm::Sp(fvc::div(phirain[phase_no]),alpharain[phase_no]) //added for bounded div scheme
    // why add?
    // bounded schemes will automatically include a source term to correct for div(U,alpha)=alpha*div(U) + U*grad(alpha)
    // material derivative is D(alpha)/dt = ddt(alpha) + U*grad(alpha) -> sub U*grad(alpha) with div(U,alpha)-alpha*div(U) but div(U)=0 at steady state, so U*grad(alpha)=div(U,alpha) when converged, but before convergence, bounded scheme will introduce a source of -alpha*div(U) so material derivative is still correct
    ==
    fvOptions(alpharain[phase_no]) // add fvOptions

    // consider scalartransportfoam equations without the laplacian since no turbulence
    // fvm::ddt(T) + fvm::div(phi, T) - fvm::laplacian(DT, T) == fvOptions(T)

);

alphaEqn.relax();

fvOptions.constrain(alphaEqn); // add fvOptions

alphaEqn.solve();

alpharain[phase_no].correctBoundaryConditions();

fvOptions.correct(alpharain[phase_no]); // add fvOptions
