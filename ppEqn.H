// Step 3: Update the mass flow rate at the cell faces using the Rhie-Chow
// interpolation technique (Eq. 15.100 & 15.201) to obtain a momentum satisfying
// mass flow rate field phi^∗f.

volScalarField rAU(1.0/UEqn.A());

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
      + (fvc::interpolate(rAU) * (gradpDiff_f & mesh.Sf()))
      + MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi))
    );

    gradpDiff_f.clear();
}

MRF.makeRelative(phi);

adjustPhi(phi, U, p);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phi, rAU, MRF);

// Step 4: Using the new mass flow rates, assemble the pressure correction
// equation (Eq. 15.98) and solve it to obtain a pressure correction field pp.

// Non-orthogonal pressure corrector loop
pp.primitiveFieldRef() = 0.0;
pp.correctBoundaryConditions();

while (piso.correctNonOrthogonal())
{
    // Pressure corrector

    fvScalarMatrix ppEqn
    (
        fvm::laplacian(rAU, pp, "laplacian(p)") == fvc::div(phi)
    );

    // set pp at the refCell such that after solving pp and subsequent update of
    // pressure, p at pRefCell is within solver tolerance of given pRefValue.
    scalar ppRefValue = (pRefValue - getRefCellValue(p, pRefCell));

    ppEqn.setReference(pRefCell, ppRefValue);

    ppEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

    // Step 5-1: Update mass flux at the cell-faces to obtain continuity
    // satisfying fields using Eq.(15.101).
    if (piso.finalNonOrthogonalIter())
    {
        phi -= ppEqn.flux();
    }
}

// #include "continuityErrs.H"

// Step 5-2: Update the pressure and velocity fields at the cell centroids

U -= rAU*fvc::grad(pp);
U.correctBoundaryConditions();
fvOptions.correct(U);

p += pp;
p.correctBoundaryConditions();