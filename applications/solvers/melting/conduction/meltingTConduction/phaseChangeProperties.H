Info<< "Reading phaseChangeProperties\n" << endl;

IOdictionary phaseChangeProperties
(
 IOobject
 (
  "phaseChangeProperties",
  runTime.constant(),
  mesh,
  IOobject::MUST_READ_IF_MODIFIED,
  IOobject::NO_WRITE
  )
 );

Info<< "Reading melting temperature Tm\n" << endl;

dimensionedScalar Tm
(
 phaseChangeProperties.lookup("Tm")
 );


Info<< "Reading delta temperature dT\n" << endl;

dimensionedScalar dT
(
 phaseChangeProperties.lookup("dT")
 );


Info<< "Reading desity rho_l\n" << endl;

dimensionedScalar rho_l
(
 phaseChangeProperties.lookup("rho_l")
 );


Info<< "Reading desity rho_s\n" << endl;

dimensionedScalar rho_s
(
 phaseChangeProperties.lookup("rho_s")
 );


// Info<< "Reading thermal conductivity k_s\n" << endl;

// dimensionedScalar k_s
// (
//  phaseChangeProperties.lookup("k_s")
//  );


// Info<< "Reading thermal conductivity k_l\n" << endl;

// dimensionedScalar k_l
// (
//  phaseChangeProperties.lookup("k_l")
//  );


Info<< "Reading thermal conductivity k\n" << endl;

dimensionedScalar k
(
 phaseChangeProperties.lookup("k")
 );



Info<< "Reading liquid thermal capacity cp_l\n" << endl;

dimensionedScalar cp_l
(
 phaseChangeProperties.lookup("cp_l")
 );


Info<< "Reading solid thermal capacity cp_s\n" << endl;

dimensionedScalar cp_s
(
 phaseChangeProperties.lookup("cp_s")
 );


Info<< "Reading specific enthalpy of melting L\n" << endl;

dimensionedScalar L
(
 phaseChangeProperties.lookup("L")
 );


Info<< "Reading Darcy parameter C\n" << endl;

dimensionedScalar C
(
 phaseChangeProperties.lookup("C")
 );


Info<< "Reading relaxation factor for g, relax_g\n" << endl;

scalar relax_g(readScalar(phaseChangeProperties.lookup("relax_g")));


Info<< "Reading relaxation factor for T, relax_T\n" << endl;

scalar relax_T(readScalar(phaseChangeProperties.lookup("relax_T")));


Info<< "Reading target residual for g, residualTarget\n" << endl;

scalar residualTarget_g(readScalar(phaseChangeProperties.lookup("residualTarget_g")));


Info<< "Reading target residual for T, residualTarget\n" << endl;

scalar residualTarget_T(readScalar(phaseChangeProperties.lookup("residualTarget_T")));


Info<< "Reading maximum number of iterations, nOuterCorrectors\n" << endl;

scalar nOuterCorrectors(readScalar(phaseChangeProperties.lookup("nOuterCorrectors")));


Info<< "Reading type of temperature-enthalpy-liquidfraction correlation, correlation\n" << endl;

label correlation(readInt(phaseChangeProperties.lookup("correlation")));


Info<< "Reading type of temperature-enthalpy-liquidfraction correlation, correlation\n" << endl;

word algorithm(word(phaseChangeProperties.lookup("algorithm")));


dimensionedScalar epsilon = dT/2.;
dimensionedScalar Ts = Tm - epsilon;
dimensionedScalar Tl = Tm + epsilon;
dimensionedScalar deltaCp = cp_l - cp_s;
dimensionedScalar deltaRho = rho_l - rho_s;
dimensionedScalar deltaRhoCp = rho_l*cp_l - rho_s*cp_s;
dimensionedScalar lambda = L-deltaCp*Tm;
dimensionedScalar dHdT_pc = (rho_l*L + deltaRho*cp_s*Tm + epsilon*(rho_s*cp_s+rho_l*cp_l))/(2*epsilon);
dimensionedScalar H_q = (0.5*(rho_s*cp_s+rho_l*cp_l)-dHdT_pc)*Tm + 0.5*rho_l*lambda;

//******************** FIELDS ******************//

//auto isPhaseTransitioning = PackedBoolList(mesh.cells().size(),false); // a list of bool to take into account if a cell has gone from T_star_m < 0 to T_star_m+1 > dT or viceversa, i.e., if a cell has skipped the phase transition.

Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field T_m\n" << endl;

volScalarField T_m
(
    IOobject
    (
        "T_m",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    T
);


Info<< "Initializing field g\n" << endl;

volScalarField g
(
 IOobject
 (
  "g",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 max(scalar(0) , min((T-Tm+epsilon)/dT,scalar(1)) )
 );


Info<< "Initializing field gm\n" << endl;

volScalarField gm
(
 IOobject
 (
  "gm",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 g
 );


Info<< "Initializing field dgdT\n" << endl;

volScalarField dgdT
(
 IOobject
 (
  "dgdT",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 g/dT
 );


Info<< "Reading field H\n" << endl;

volScalarField H
(
    IOobject
    (
        "H",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho_s*cp_s*T + deltaRhoCp*g*T + rho_l*lambda*g
);


Info<< "Initializing field dHdT\n";

volScalarField dHdT
(
 IOobject
 (
  "dHdT",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 (T<(Tm-epsilon)*T_m/T_m)? rho_s*cp_s*T/T : ((T>(Tm+epsilon)*T_m/T_m)? rho_l*cp_l*T/T : dHdT_pc*T/T)
 );


Info<< "Initializing field rho\n" << endl;

volScalarField rho
(
 IOobject
 (
  "rho",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 (1.0-g)*rho_s+g*rho_l
 );


Info<< "Initializing field mean density, f \n" << endl;

volScalarField f
(
 IOobject
 (
  "f",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 g*rho_l/rho
 );


// Info<< "Reading field mean k \n" << endl;

// volScalarField k
// (
//  IOobject
//  (
//   "k",
//   runTime.timeName(),
//   mesh,
//   IOobject::NO_READ,
//   IOobject::NO_WRITE
//   ),
//  k_s*(1.0-g) + k_l*g
//  );
