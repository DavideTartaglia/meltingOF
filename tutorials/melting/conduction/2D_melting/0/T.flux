/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  7
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      T;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 1 0 0 0];

internalField   uniform 260;

boundaryField
{
    north
    {
        type            fixedGradient;
        gradient        3e3;
    }

    ovest
    {
        type            zeroGradient;
    }


    east
    {
        type            zeroGradient;
    }


    south
    {
        type            fixedValue;
        value           uniform 260;
    }

    frontAndBack
    {
        type            empty;
    }
}


// ************************************************************************* //
