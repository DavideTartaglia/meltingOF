    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    //volScalarField alphaEff("alphaEff", DT + alphat);

    fvScalarMatrix TEqn
    (
        fvm::div(fvc::interpolate(rho_cp)*phi, T)
      - fvm::laplacian(k,T)
      - fvm::laplacian(rho_l*cp_l*alphat, T) 
     ==
        rho_l*cp_l*radiation->ST(rhoCpRef, T)
      + fvOptions(rho_cp,T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    radiation->correct();

    fvOptions.correct(T);


#include "calculate_f.H"
#include "update_thermalProperties.H"
rhok = 1.0 - beta*(T - TRef);

