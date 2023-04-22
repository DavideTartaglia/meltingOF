Info<< "\nwriting coefficients of the source, Su + Sp*T" << endl;

const dimensionedScalar dt_inv = 1.0/runTime.deltaT();
const volScalarField Su = - dt_inv*d_rhoh_dT*T_starm;
//Info<< "dimensions Su: " << Su.dimensions() << endl;
//Info<< "Su = " << Su << endl;
const volScalarField Sp = dt_inv*d_rhoh_dT;
//Info<< "dimensions Sp: " << Sp.dimensions() << endl;
//Info<< "Sp = " << Sp << endl;
