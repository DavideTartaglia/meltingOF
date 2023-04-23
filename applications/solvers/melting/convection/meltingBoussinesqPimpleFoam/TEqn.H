phaseChange.updateIteration();
phaseChange.updateFields();

dimensionedScalar delta_t = runTime.deltaT();

// volScalarField T_star = T-Tm;
// Info<<"\nT_star[0]:\t"<<T_star[0];
// Info<<"\nT_star[1]:\t"<<T_star[1];
// Info<<"\ng[0]:\t"<<g[0];
// Info<<"\ng[1]:\t"<<g[1];
// Info<<"\nH[0]:\t" << H[0];
// Info<<"\nH[1]:\t" << H[1];
// Info<<"\ndHdT[0]:\t"<< dHdT[0];
// Info<<"\ndHdT[1]:\t"<< dHdT[1];
// Info<<"\n\n"
// Info<<"\nT_starm:\t"<< T_starm;
// Info<< "\ndgdT:\t"<< dgdT;
// Info<< "\nT:\t"<< T;

{
  volScalarField rho_cp = rho_s*cp_s + liquidFraction*rho_l*deltaCp;

  alphat = turbulence->nut()/Prt;
  alphat.correctBoundaryConditions();

  //volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

  fvScalarMatrix TEqn
    (
     fvm::Sp(phaseChange.Sp()/delta_t,T) + phaseChange.Su()/delta_t // = dH/dt
     + rho_l*cp_l*fvm::div(phi, T)
     - fvm::laplacian(k,T)
     - fvm::laplacian(rho_l*cp_l*alphat, T)
     ==
     /* rho_cp*radiation->ST(rho_l*cp_l,T)
	+ */rho_l*cp_l*fvOptions(T)
     );

  TEqn.relax();
  fvOptions.constrain(TEqn);
  TEqn.solve();
  radiation->correct();
  fvOptions.correct(T);

}

phaseChange.relax_T();

phaseChange.update();

phaseChange.checkResidual();

rhok = 1.0 - beta*(T - TRef);

