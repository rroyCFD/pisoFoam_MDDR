// PRIME loop Step 6: Using the latest available velocity, flux and pressure
// fields, re-formulate  momentum equation and solve it explicitly.

tmp<fvVectorMatrix> tUPrimeEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
);
fvVectorMatrix& UPrimeEqn = tUPrimeEqn.ref();

volScalarField rAU(1.0/UPrimeEqn.A());

// explicit velocity calculation (eq. 15.176)
U = rAU*(UPrimeEqn.H() - fvc::grad(p));
tUPrimeEqn.clear();


// Rhie-Chow and ddt flux correction
{
    tmp<surfaceVectorField> gradpDiff_f
    (
        (  (fvc::interpolate(fvc::grad(p)) & ed)
          - faceGradient.snGrad(p)
        ) * ed
     );

    phi =
    (   fvc::flux(U)
      + (fvc::interpolate(rAU)* (gradpDiff_f & mesh.Sf()))
      + MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi))
    );

    gradpDiff_f.clear();
}

MRF.makeRelative(phi);

adjustPhi(phi, U, p);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phi, rAU, MRF);


// Step 8: Using the new mass flow rates, assemble the pressure correction
// equation (Eq. 15.183) and solve it to obtain a pressure correction field.

// Non-orthogonal pressure corrector loop
pp.primitiveFieldRef() = 0.0;
pp.correctBoundaryConditions();

while (piso.correctNonOrthogonal())
{
    // Pressure corrector

    fvScalarMatrix ppPrimeEqn
    (
        fvm::laplacian(rAU, pp, "laplacian(p)") == fvc::div(phi)
    );

    // set pp at the refCell such that after solving pp and subsequent update of
    // pressure, p at pRefCell is within solver tolerance of given pRefValue.
    scalar ppRefValue = (pRefValue - getRefCellValue(p, pRefCell));

    ppPrimeEqn.setReference(pRefCell, ppRefValue);

    ppPrimeEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

    // Step 9-1: Update mass flux at the cell-faces to obtain continuity
    // satisfying fields using Eq.(15.101).
    if (piso.finalNonOrthogonalIter())
    {
        phi -= ppPrimeEqn.flux();
    }
}

// #include "continuityErrs.H"

// Step 9-2: Update the pressure and velocity fields at the cell centroids

U -= rAU*fvc::grad(pp);
U.correctBoundaryConditions();
fvOptions.correct(U);

p += pp;
p.correctBoundaryConditions();