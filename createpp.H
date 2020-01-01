//----------------------------------------------------------------------------ss//
Info << "\nCreating pressure-coorection field pp" << endl;
wordList pbt = p.boundaryField().types();

// prepare the boundary condition types for pressure-correction field
wordList ppbt = pbt;
forAll(ppbt, b)
{
    if(ppbt[b] == "timeVaryingMappedFixedValue")
    {
        ppbt[b] = "fixedValue";
    }
}

volScalarField pp
(
    IOobject
    (

        "pp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT, // NO_READ,
        IOobject::AUTO_WRITE // NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", p.dimensions(), 0.0),
    ppbt
);
Info << "Pressure correction field boundary types:\n"
     << pp.boundaryField().types() << endl;


// Set pp boundary values
forAll(pp.boundaryField(), patchi)
{
    if (isType<fixedValueFvPatchScalarField>(pp.boundaryField()[patchi]))
    {
        fixedValueFvPatchScalarField& ppbound =
           refCast<fixedValueFvPatchScalarField>(pp.boundaryFieldRef()[patchi]);

        ppbound == scalarField(ppbound.size(),0.0);
    }

}