List<scalar> inletAlpha
(
    transportProperties.lookup("inletAlpha")
);


Info<< nl << "Calculating specific catch ratio fields\n" << endl;
forAll (phases, phase_no)
{
    dimensionedScalar Vt_i("Vt_i", dimLength/dimTime, Vt[phase_no]);

    surfaceScalarField normalvel
    (
         (mesh.Sf()/mesh.magSf()) & fvc::interpolate(Urain[phase_no]) // remove mag since normals will point out of surface
    );
    surfaceScalarField surfaceScr
    (
        (normalvel * fvc::interpolate(alpharain[phase_no])) / (inletAlpha[phase_no] * Vt_i) //Here 1.0 is the alpha values at inlets // set inletAlpha in transportProperties and read from there
    );
    const surfaceScalarField::Boundary& patchSurfaceScr = surfaceScr.boundaryField();

    word namescr_i ("scr" + name(phase_no+1));
    scr.set
    (
    phase_no,
    new volScalarField
    (
        IOobject
        (
            namescr_i,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
         (Urain[phase_no] & g.value()/mag(g.value()) * alpharain[phase_no]) / (inletAlpha[phase_no] * Vt_i) // write volScalarField as (Uz*alpha)/(inletAlpha*inletVt)
        //mesh, // write blank mesh
        // dimensionedScalar(namescr_i, dimless, 0.0)
    )
    );        
    volScalarField::Boundary& scrtempBf = scr[phase_no].boundaryFieldRef();
    
    forAll(scrtempBf,patchi)
    {
        scrtempBf[patchi] = patchSurfaceScr[patchi];
    }    

    scr[phase_no].write();
}

Info << "Calculating catch ratio field\n" << endl;
forAll (Rh, i)
{
    word gcrname ("gcr");
    if (Rh.size()!=1) {
        gcrname = "gcr" + name(Rh[i]);
    }
    volScalarField gcrtemp
    (
        IOobject
        (
            gcrname,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("gcr", dimless, 0.0)
    );

    forAll (phases, phase_no)
    {
        gcrtemp += scr[phase_no]*DFh[i][phase_no];
    }
    gcrtemp.write();
}
