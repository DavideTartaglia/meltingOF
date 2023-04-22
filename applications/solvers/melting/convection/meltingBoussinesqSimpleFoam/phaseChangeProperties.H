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


Info<< "Reading melting temperature dT\n" << endl;

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


Info<< "Reading thermal conductivity k_s\n" << endl;

dimensionedScalar k_s
(
 phaseChangeProperties.lookup("k_s")
 );


Info<< "Reading thermal conductivity k_l\n" << endl;

dimensionedScalar k_l
(
 phaseChangeProperties.lookup("k_l")
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


Info<< "Reading relaxation factor for f, relax_f\n" << endl;

dimensionedScalar relax_f
(
 phaseChangeProperties.lookup("relax_f")
 );


//******************** FIELDS ******************//


Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field T_star\n" << endl;
volScalarField T_star
(
    IOobject
    (
        "T_star",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    T-Tm
);


Info<< "Reading field f\n" << endl;

volScalarField f
(
 IOobject
 (
  "f",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 max( scalar(0) , min(T_star/dT, scalar(1)) )
 );
    

Info<< "Reading field mean cp, cp \n" << endl;

volScalarField cp
(
 IOobject
 (
  "cp",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
  cp_s*(1.0-f)+cp_l*f
 );


Info<< "Reading field mean enthalpy, h \n" << endl;

volScalarField h
(
 IOobject
 (
  "h",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
  cp*T_star + f*L
 );


Info<< "Reading field mean k, k \n" << endl;

volScalarField k
(
 IOobject
 (
  "k",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
  k_s*(1.0-f)+k_l*f
 );


Info<< "Reading field mean rho, rho \n" << endl;

volScalarField rho
(
 IOobject
 (
  "rho",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
  1.0/((1.0-f)/rho_s+f/rho_l)
 );


Info<< "Initializing field mean rho_cp, rho_cp \n" << endl;

volScalarField rho_cp
(
 IOobject
 (
  "rho_cp",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::NO_WRITE
  ),
 rho*cp
 );
