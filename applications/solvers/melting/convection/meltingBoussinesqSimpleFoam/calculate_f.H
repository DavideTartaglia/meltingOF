Info<< "Updating f" << endl;

  const volScalarField f_correction = relax_f*cp*(T - Tm)/L;
  const volScalarField fnew = f + f_correction;
  f = max( scalar(0) , min(fnew, scalar(1)) ) ;

/*
forAll(mesh.C(),i)
{
  scalar Ti = T[i];
  scalar cp_meani = cp_mean[i];
  scalar fnew = f[i] + cp_meani*(Ti - Tm.value())/L.value();
  f[i] = min(1,max(0, fnew));
}


forAll(mesh.boundary(),patchID)
{

  scalar Ti = T[patchID];
  scalar cp_meani = cp_mean[patchID];
  scalar fnew = f[patchID] + cp_meani*(Ti - Tm.value())/L.value();
  f[patchID] = min(1,max(0, fnew));
*/


/*
  // f[i] = std::min(1,std::max(0,((T[i]-Tm.value())/r.value())));

  if (T[patchID]<Ts.value())
    {
      f[patchID] = 0;
    }
  else if (T[patchID]>Tl.value())
    {
      f[patchID] = 1;
    }
  else if (T[patchID]>=Ts.value() && T[patchID]<=Tl.value())
    {
      f[patchID] = (T[patchID]-Ts.value())/(Tl.value()-Ts.value());
    }
*/

/*
  forAll (mesh.boundary()[patchID],facei)
    {

  	scalar Ti = T.boundaryField()[patchID][facei];
  	scalar cp_meani = cp_mean.boundaryField()[patchID][facei];
  	scalar fnew = f.boundaryFieldRef()[patchID][facei] + cp_meani*(Ti - Tm.value())/L.value();
  	f.boundaryFieldRef()[patchID][facei] = min(1,max(0, fnew));
*/
/*
      if (T.boundaryField()[patchID][facei]<Ts.value())
	{
	  f.boundaryFieldRef()[patchID][facei] = 0;
	}
      else if (T.boundaryField()[patchID][facei]>Tl.value())
	{
	  f.boundaryFieldRef()[patchID][facei] = 1;
	}
      else if (T.boundaryField()[patchID][facei]>=Ts.value() && T.boundaryField()[patchID][facei]<=Tl.value())
	{
	  f.boundaryFieldRef()[patchID][facei] = (T.boundaryField()[patchID][facei]-Ts.value())/(Tl.value()-Ts.value());
	} */
//    }
    
//}

