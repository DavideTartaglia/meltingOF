/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Fo1undation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    meltingTConduction

Description
    Solves the optimum approach method for melting and solidification problems.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "phaseChangeClassControl.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"
    #include "createControl.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {

      runTime++;
      Info<< "\nTime = " << runTime.timeName() << nl << endl;
 
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
      	  
      while (phaseChange.loop())
	{
	  phaseChange.updateIteration();
	  do
	    {
	      phaseChange.updateFields();
	      dimensionedScalar delta_t = runTime.deltaT();
	      volScalarField T_star = T-Tm;

	      Info<< "\nTstar[0] before solving: " << T_star[0] << endl
		  << "Tstar[1] before solving: " << T_star[1] << endl
		  << "T.oldTime[0]: " << T.oldTime()[0];
	  
	      Info<< "\nCalculating temperature distribution";

	      while (simple.correctNonOrthogonal())
		{
		  fvScalarMatrix TEqn
		    (
		     fvm::SuSp(phaseChange.Sp()/delta_t,T) + phaseChange.Su()/delta_t // dHdt, dependent on the algorithm.
		     // fvm::Sp(dHdT/delta_t,T)*
		     // + (fvc::Su((-dHdT*T_m + H - H.oldTime())/delta_t,T))*(algorithm==(implicit)) //linearization of ddt(H)
		     // + (fvm::Su((-dHdT*T.oldTime())*(algorithm==(apparentCp))
		     - fvm::laplacian(k,T)
		     ==
		     (rho_s*cp_s + g*deltaRhoCp)*fvOptions(T)
		     );

		  //Info<< "\nTeqn:\n" << TEqn;

		  fvOptions.constrain(TEqn);
		  TEqn.solve();
		  fvOptions.correct(T);
		}

	      Info<< "T after solve:\t"<< T[0];
	      T_star = T-Tm;
	      Info<< "\nT_star[0] after solve: "<< T_star[0]
		  << "\nTstar[1] before solving: " << T_star[1] << endl;
	    }
	  while(phaseChange.dynamicChoice());

	  phaseChange.relax_T();
	  
	  Info<< "\ng before update: "<< g[0];
          phaseChange.update();
	  Info<< "\ng_updated: " << g[0];

	  phaseChange.checkResidual();

	  Info<< "\n";
	  
	}
      
      #include "write.H"

      Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	  << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	  << nl << endl;

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
