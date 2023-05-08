#ifndef _PHASECHANGECLASS_H
#define _PHASECHANGECLASS_H
#include "fvCFD.H"

class phaseChangeControl
{  
private:
  dimensionedScalar& my_Tm;
  dimensionedScalar& my_dT;
  dimensionedScalar& my_rho_l;
  dimensionedScalar& my_rho_s;
  dimensionedScalar& my_k;
  dimensionedScalar& my_cp_l;
  dimensionedScalar& my_cp_s;
  dimensionedScalar& my_L;
  dimensionedScalar& my_C;
  scalar& my_relax_g;
  scalar& my_relax_T;
  scalar& my_residualTarget_g;
  scalar& my_residualTarget_T;
  scalar& my_nOuterCorrectors;
  label& my_correlation;
  word& my_algorithm;
  bool dynamicAlgorithm;
  
  dimensionedScalar& my_epsilon;
  dimensionedScalar& my_Ts;
  dimensionedScalar& my_Tl;
  dimensionedScalar& my_deltaCp;
  dimensionedScalar& my_deltaRho;
  dimensionedScalar& my_deltaRhoCp;
  dimensionedScalar& my_lambda;
  dimensionedScalar& my_dHdT_pc;
  dimensionedScalar& my_H_q;

  // auto isPhaseTransitioning = PackedBoolList(mesh.cells().size(),false); a list of bool to take into account if a cell has gone from T_star_m < 0 to T_star_m+1 > dT or viceversa, i.e., if a cell has skipped the phase transition.

  volScalarField& my_T;
  volScalarField& my_T_m;
  volScalarField& my_H;
  volScalarField& my_dHdT;
  volScalarField& my_g;
  volScalarField& my_gm;
  volScalarField& my_dgdT;
  volScalarField& my_rho;
  volScalarField& my_f;
  
  label iterCounter;
  bool isLastIter;
  scalar SMALL;
  dimensionedScalar sumFieldError_g;
  dimensionedScalar measure_g;
  dimensionedScalar residual_g;
  dimensionedScalar sumFieldError_T;
  dimensionedScalar measure_T;
  dimensionedScalar residual_T;

  enum {LINEARg, LINEARH, SIGMOIDg, NUMERICAL, NEWTON};

  void update_derivatives()
  {
    Info<< "\nUpdating derivatives dHdT, dgdT";
    try
      {
	switch(my_correlation)
	  {
	  case LINEARg:
	    {
	      forAll(my_T,i)
		{
		  scalar& Ti = my_T[i];
		  scalar& dgdTi = my_dgdT[i];
		  scalar& Ts = my_Ts.value();
		  scalar& Tl = my_Tl.value();
		  
		  if(Ti<Ts)
		    {
		      dgdTi = 0;
		    }
		  else if (Ti>Tl)
		    {
		      dgdTi = 0;
		    }
		  else if (Ti>=Ts && Ti<=Tl)
		    {
		      dgdTi = 1.0/(Tl-Ts);
		    }
		  else
		    {
		      Info<< "\n\nunkown error in phaseChangeControlClass.H\n\n";
		      std::exit(1);
		    }
		}
	      my_dHdT = my_rho_s*my_cp_s + my_g*my_deltaRhoCp + (my_deltaRhoCp*my_T + my_rho_l*my_lambda)*my_dgdT;
	      break;
	    }
	  case LINEARH:
	    {
	      Info<<"\nLinearH update of derivatives.";
	      forAll(my_T,i)
		{
		  scalar& Ti = my_T[i];
		  scalar& dHdTi = my_dHdT[i];
		  scalar& dgdTi = my_dgdT[i];
		  scalar& cps = my_cp_s.value();
		  scalar& rhos= my_rho_s.value();
		  scalar& cpl = my_cp_l.value();
		  scalar& rhol= my_rho_l.value();
		  scalar& Ts = my_Ts.value();
		  scalar& Tl = my_Tl.value();
		  scalar& dHdT_pc = my_dHdT_pc.value();
		  
		  if(Ti<Ts)
		    {
		      dHdTi = cps*rhos;
		      dgdTi = 0;
		    }
		  else if (Ti>Tl)
		    {
		      dHdTi = cpl*rhol;
		      dgdTi = 0;
		    }
		  else if (Ti>=Ts && Ti<=Tl)
		    {
		      dHdTi = dHdT_pc;
		      dgdTi = 1;
		    }
		  else
		    {
		      Info<< "\n\nunkown error in phaseChangeControlClass.H\n\n";
		      std::exit(1);
		    }
		}
	      break;
	    }
	  case SIGMOIDg:
	    {
	      my_dgdT = my_g*(1.0-my_g)/my_dT;
	      my_dHdT = my_rho_s*my_cp_s + my_g*my_deltaRhoCp + (my_deltaRhoCp*my_T + my_rho_l*my_lambda)*my_dgdT;
	      break;
	    }
	  case NUMERICAL:
	    {
	      if(my_algorithm=="apparentCp")
		{
		  throw my_correlation;
		}
	      my_dHdT = (my_H - my_H.oldTime())/(my_T - my_T.oldTime() + 1e-9*my_dT);
	      my_dgdT = (my_T<my_Ts*(my_T/my_T))? 0/my_dT : (my_T>my_Tl*(my_T/my_T))? 0/my_dT : 1./my_dT;
	      break;
	    }
	  case NEWTON:
	    {
	      forAll(my_H,i)
		{
		  scalar& Hi = my_H[i];
		  scalar& dHdTi = my_dHdT[i];
		  scalar& dgdTi = my_dgdT[i];
		  scalar Hs = my_dHdT_pc.value()*my_Ts.value() + my_H_q.value();
		  scalar Hl = my_dHdT_pc.value()*my_Tl.value() + my_H_q.value();
		  scalar& cps = my_cp_s.value();
		  scalar& rhos= my_rho_s.value();
		  scalar& cpl = my_cp_l.value();
		  scalar& rhol= my_rho_l.value();
		  scalar& dHdT_pc = my_dHdT_pc.value();
		  
		  if(Hi<Hs)
		    {
		      dHdTi = cps*rhos;
		      dgdTi = 0;
		    }
		  else if (Hi>Hl)
		    {
		      dHdTi = cpl*rhol;
		      dgdTi = 0;
		    }
		  else if (Hi>=Hs && Hi<=Hl)
		    {
		      dHdTi = dHdT_pc;
		      dgdTi = 1;
		    }
		  else
		    {
		      Info<< "\n\nunkown error in phaseChangeControlClass.H\n\n";
		      std::exit(1);
		    }
		}
	      break;
	    }
	  default:
	    {
	      throw my_correlation;
	      break;
	    }
	  }
      }
    catch(label correlation)
      {
	this->exceptionHandler(correlation);
      }
  };

  void exceptionHandler(label correlation)
  {
    if(correlation==3 && my_algorithm=="apparentCp")
      {
	Info<< "\n\n\n****************************************\n"
	  << "PHASE CHANGE ERROR: apparentCp algorithm is incompatible with <<3>> correlation\n"
	  << "i.e., apparent Cp can not be used with a numerical estimation of dHdT\n\n\n\n";
	     std::exit(1);
      }
    if((!(correlation==0)) && my_algorithm=="source")
      {
	Info<< "\n\n\n****************************************\n"
	    << "PHASE CHANGE ERROR: source algorithm can only be used with <<0>> correlation\n"
	    << "i.e. source algorithm => linear g-T correlation.\n"
	    << "Please try again by changing either the correlation or the algorithm.\n\n\n\n";
	std::exit(1);
      }
    Info<< "\n\n\n********************************************\n"
	<<"PHASE CHANGE ERROR: INVALID g-H-T correlation. Valid Correlations:\n"
	<<"0 -> linear g-T\n"
	<<"1 -> linear H-T\n"
	<<"2 -> sigmoid g-T\n"
	<<"3 -> numerical calc of dHdT\n"
	<<"You tried: correlation = "<< correlation
	<<"\nThe number of correlation does not exist.\n"
	<<"\nPlease check phaseChangeProperties dict or the source code.\n\n\n" << endl;
    std::exit(1);
  };

  void exceptionHandler(word algorithm)
  {
    if(!(algorithm=="source" || algorithm=="apparentCp" || algorithm=="implicit"))
      {
	Info<< "\n\n\n****************************************\n"
	    << "PHASE CHANGE ERROR: invalid algorithm. Valid algorithms are:\n"
	    << "source\n"
	    << "apparentCp\n"
	    << "implicit\n"
	    << "you tried:\t"<< algorithm
	    << endl;
	std::exit(1);
      }
  };
  
  //********************** PUBLIC MEMBER FUNCTIONS ***********************************//
public:
  phaseChangeControl(dimensionedScalar& Tm,
		     dimensionedScalar& dT,
		     dimensionedScalar& rho_l,
		     dimensionedScalar& rho_s,
		     dimensionedScalar& k,
		     dimensionedScalar& cp_l,
		     dimensionedScalar& cp_s,
		     dimensionedScalar& L,
		     dimensionedScalar& C,
		     scalar& relax_g,
		     scalar& relax_T,
		     scalar& residualTarget_g,
		     scalar& residualTarget_T,
		     scalar& nOuterCorrectors,
		     label& correlation,
		     word& algorithm,
		     
		     dimensionedScalar& epsilon,
		     dimensionedScalar& Ts,
		     dimensionedScalar& Tl,
		     dimensionedScalar& deltaCp,
		     dimensionedScalar& deltaRho,
		     dimensionedScalar& deltaRhoCp,
		     dimensionedScalar& lambda,
		     dimensionedScalar& dHdT_pc,
		     dimensionedScalar& H_q,
		     
		     volScalarField& T,
		     volScalarField& T_m,
		     volScalarField& H,
		     volScalarField& dHdT,
		     volScalarField& g,
		     volScalarField& gm,
		     volScalarField& dgdT,
		     volScalarField& rho,
		     volScalarField& f):
    my_Tm(Tm),
    my_dT(dT),
    my_rho_l(rho_l),
    my_rho_s(rho_s),
    my_k(k),
    my_cp_l(cp_l),
    my_cp_s(cp_s),
    my_L(L),
    my_C(C),
    my_relax_g(relax_g),
    my_relax_T(relax_T),
    my_residualTarget_g(residualTarget_g),
    my_residualTarget_T(residualTarget_T),
    my_nOuterCorrectors(nOuterCorrectors),
    my_correlation(correlation),
    my_algorithm(algorithm),
    dynamicAlgorithm((my_algorithm=="dynamic")? true : false),

    my_epsilon(epsilon),
    my_Ts(Ts),
    my_Tl(Tl),
    my_deltaCp(deltaCp),
    my_deltaRho(deltaRho),
    my_deltaRhoCp(deltaRhoCp),
    my_lambda(lambda),
    my_dHdT_pc(dHdT_pc),
    my_H_q(H_q),

    my_T(T),
    my_T_m(T_m),
    my_H(H),
    my_dHdT(dHdT),
    my_g(g),
    my_gm(gm),
    my_dgdT(dgdT),
    my_rho(rho),
    my_f(f),

    iterCounter(0),
    isLastIter(false),
    SMALL(1e-10),
    sumFieldError_g(sumMag(g.internalField() - gm.internalField())\
		  + sumMag(g.boundaryField() - gm.boundaryField())),
    measure_g(sumMag(g.internalField()) + sumMag(g.boundaryField()) + SMALL),
    residual_g(sumFieldError_g/measure_g),
    sumFieldError_T(sumMag(T.internalField() - T_m.internalField())),
    measure_T(sumMag(T.internalField())),
    residual_T(sumFieldError_T/measure_T)
  {
    if(dynamicAlgorithm)
      {
	this->setAlgorithm("apparentCp");
      }
    
    if(my_correlation==SIGMOIDg)
      {
	my_dT = my_dT/10.;
	Info<< "\ndT: " << my_dT << endl;
	my_epsilon = my_dT/2.;
	Info<< "\nepsilon: "<< my_epsilon << endl;
	scalar my_offset = 0.0;
	my_g = 1.0/(1.0 + exp(-(my_T-my_Tm)/my_dT + my_offset));
	my_H = my_rho_s*my_cp_s*my_T + my_deltaRhoCp*my_g*my_T + my_rho_l*my_lambda*my_g;
      }
    
    this->exceptionHandler(my_algorithm);
  };


  bool loop()
  {
    if(isLastIter==true && iterCounter==0)
      {
	Info<< "\nPhase change loop ended.\n\n";
	isLastIter = false;
	return false;
      }
    else
      {
	return true;
      }
  };
  

  void checkResidual()
  {
    sumFieldError_g = sumMag(my_g.internalField() - my_gm.internalField())\
      + sumMag(my_g.boundaryField() - my_gm.boundaryField());
    measure_g = sumMag(my_g.internalField()) + sumMag(my_g.boundaryField()) + SMALL;
    residual_g = sumFieldError_g/measure_g;
    sumFieldError_T = sumMag(my_T.internalField() - my_T_m.internalField());
    measure_T = sumMag(my_T.internalField());
    residual_T = sumFieldError_T/measure_T;
    Info<< "\ng residual: " << residual_g;
    Info<< "\nT residual: " << residual_T;
    
    if(isLastIter==true)
      {
	iterCounter = 0;
      }
    if(isLastIter==false && (residual_g.value()<my_residualTarget_g && residual_T.value()<my_residualTarget_T))
      {
	Info<< "\n\nPhase change loop converged. Doing final iteration.";
	isLastIter = true;
      }
    if(isLastIter==false && iterCounter>=my_nOuterCorrectors)
      {
	Info<<"\n\nPhase change loop did not converge in "<< my_nOuterCorrectors << " iterations. Doing final iteration.";
	isLastIter = true;
      }
    if(my_algorithm=="apparentCp")
      {
	Info<<"\napparentCp model does not need any iteration.";
	isLastIter = true;
	iterCounter = 0;
      }
  };

  
  void updateIteration()
  {
    iterCounter++;
    Info<< "\nPhase change iter: "<< iterCounter;

    if(dynamicAlgorithm && iterCounter==1)
      {
	this->setAlgorithm("apparentCp");
      }
  };


  void relax_T()
  {
    my_T = my_relax_T*my_T + (1.0 - my_relax_T)*my_T_m;
  };


  volScalarField Sp()
  {
	// if T.mesh().steady() is false (ie if the ddt fvScheme is not transient), then
	if(!my_T.mesh().steady())
	{
		if(my_algorithm==("implicit") || my_algorithm==("apparentCp"))
		{
		Info<< "\nSp: implicit OR apparentCp algorithm\n";
		return my_dHdT;
		}
		else if(my_algorithm==("source"))
		{
		Info<< "\nSp: source method\n";
		return (my_rho_s*my_cp_s + my_deltaRhoCp*my_g);
		}
	}
	else
	{ return my_dHdT*0.0; }
  }

  volScalarField Su()
  {
	if(!my_T.mesh().steady())
	{
		if(my_algorithm==("implicit"))
		{
		Info<< "\nSu: implicit algorithm";
		return (-my_dHdT*my_T_m + my_H - my_H.oldTime());
		}
		else if (my_algorithm==("apparentCp"))
		{
		Info<< "\nSu: apparentCp algorithm";
		return (-my_dHdT*my_T.oldTime());
		}
		else if (my_algorithm==("source"))
		{
		Info<< "\nSu: source algorithm";
		return (- my_rho_s*my_cp_s*my_T.oldTime()\
			- my_deltaRhoCp*my_g.oldTime()*my_T.oldTime()\
			+ my_rho_l*my_lambda*(my_g - my_g.oldTime()));
		}
	}
	else
	{ return my_H*0.0; }
  }
      
  

  void update()
  {
    Info<< "\nupdating g";
    if(!(my_algorithm=="source"))
      {
	try
	  {
	    switch(my_correlation)
	      {
	      case LINEARg:
		{
		  Info<< "\nLinear g-T correlation updating";
		  if(isLastIter) //if we are at the last iteration, the relaxation factor is set to 0.
		    {
		      my_g = max( scalar(0) , min((my_T - my_Tm + my_epsilon)/my_dT, scalar(1)) ) ; // limiting 0 < g < 1
		      my_H = my_rho_s*my_cp_s*my_T + my_deltaRhoCp*my_g*my_T + my_rho_l*my_lambda*my_g;
		    }
		  else
		    {
		      my_g = (1.0-my_relax_g)*my_g + my_relax_g*max( scalar(0) , min((my_T - my_Tm + my_epsilon)/my_dT, scalar(1)) ) ;
		      my_H = my_rho_s*my_cp_s*my_T + my_deltaRhoCp*my_g*my_T + my_rho_l*my_lambda*my_g;
		    }
		  break;
		}
	      case LINEARH:
		{
		  Info<< "\nLinear H-T correlation updating";
		  volScalarField H_tmp = max(my_rho_s*my_cp_s*my_T , min(my_rho_l*my_cp_l*my_T + my_rho_l*my_lambda , my_dHdT_pc*my_T + my_H_q));
		  my_H = H_tmp;
		  my_g = max( scalar(0) , min((my_T - my_Tm + my_epsilon)/my_dT, scalar(1)) ) ; 

		  if(!isLastIter) //if we are at the last iteration, the relaxation factor is set to 0.
		    {
		      my_H = (1.0-my_relax_g)*my_H + my_relax_g*H_tmp;
		    }
		  break;
		}
	      case SIGMOIDg: 
		{
		  Info<< "\nSigmoid correlation updating";
		  scalar my_offset = 0.0;
		  volScalarField g_tmp = my_g;
		  my_g = 1.0/(1.0 + exp(-(my_T-my_Tm)/my_dT + my_offset));
		  my_g = g_tmp*(1.0-my_relax_g) + my_relax_g*my_g;
		  my_H = my_rho_s*my_cp_s*my_T + my_deltaRhoCp*my_g*my_T + my_rho_l*my_lambda*my_g;
		  break;
		}
	      case NUMERICAL:
		{
		  Info<< "\nNumerical correlation updating";
		  volScalarField H_tmp = max(my_rho_s*my_cp_s*my_T , min(my_rho_l*my_cp_l*my_T + my_rho_l*my_lambda , my_dHdT_pc*my_T + my_H_q));
		  my_H = (1.0-my_relax_g)*my_H + my_relax_g*H_tmp;
		  my_g = max( scalar(0) , min((my_T - my_Tm + my_epsilon)/my_dT, scalar(1)) ) ;
		  break;
		}
	      case NEWTON:
		{
		  Info<< "\nNewton updating";
		  dimensionedScalar Hs = my_dHdT_pc*my_Ts + my_H_q;
		  dimensionedScalar Hl = my_dHdT_pc*my_Tl + my_H_q;
		  
		  forAll(my_H,i)
		    {
		      scalar& Hi = my_H[i];
		      scalar& T_mi = my_T_m[i];
		      scalar Hs = my_dHdT_pc.value()*my_Ts.value() + my_H_q.value();
		      scalar Hl = my_dHdT_pc.value()*my_Tl.value() + my_H_q.value();
		      scalar& cps = my_cp_s.value();
		      scalar& rhos= my_rho_s.value();
		      scalar& cpl = my_cp_l.value();
		      scalar& rhol= my_rho_l.value();
		      scalar& dHdT_pc = my_dHdT_pc.value();
		      scalar& Hq = my_H_q.value();
		      scalar& lambda = my_lambda.value();
		  
		      if(Hi<Hs)
			{
			  T_mi = Hi/cps/rhos;
			}
		      else if (Hi>Hl)
			{
			  T_mi = (Hi-rhol*lambda)/rhol/cpl;
			}
		      else if (Hi>Hs && Hi<Hl)
			{
			  T_mi = (Hi-Hq)/dHdT_pc;
			}
		      else
			{
			  Info<< "\n\nunkown error in phaseChangeControlClass.H\n\n";
			  std::exit(1);
			}
		    }
		  Info<< "\ndeltaH[0] correction:\t" << my_relax_g*my_dHdT[0]*(my_T[0]-my_T_m[0]);
		  Info<< "\ndeltaH[1] correction:\t" << my_relax_g*my_dHdT[1]*(my_T[1]-my_T_m[1]);
		  my_H = my_H + my_relax_g*my_dHdT*(my_T - my_T_m);
		  my_g = max( scalar(0) , min((my_H-Hs)/(Hl-Hs), scalar(1)) ) ;
		  break;
		}
	      default:
		{
		  throw my_correlation;
		  break;
		}
	      }
	  }
	catch(label correlation)
	  {
	    this->exceptionHandler(correlation);
	  }
      }
    else if(my_algorithm=="source")
      {
	try
	  {
	    if(my_correlation==0)
	      {
		volScalarField my_cp = my_cp_s*(1.-my_f)+my_cp_l*my_f;
		my_g = my_g + my_relax_g*my_cp*(my_T - my_Tm + my_epsilon - 2.0*my_epsilon*my_g)/(my_L);
		my_g = max( scalar(0) , min(my_g, scalar(1)) ) ;
	      }
	    else
	      {
		throw my_correlation;
	      }
	  }
	catch(label correlation)
	  {
	    this->exceptionHandler(correlation);
	  }
      }
  };


  void updateFields()
  {
    Info<< "\nupdating variable fields";
    dimensionedScalar Hs = my_dHdT_pc*my_Ts + my_H_q;
    dimensionedScalar Hl = my_dHdT_pc*my_Tl + my_H_q;
    Info<< "\nHs = " << Hs.value();
    Info<< "\nHl = " << Hl.value();
    Info<< "\nHq = " << my_H_q.value();
    //f = rho_l*rhof/(rhof*(rho_l - rho_s)+rho_s*rho_l);
    //T = T_star + Tm;
    // dfdT = f*(1.0-f)/dT;

    if(!(my_algorithm=="source"))
      {
	if(my_algorithm=="implicit")
	  {
	    this->update_derivatives();
	  }
	if(my_algorithm=="apparentCp")
	  {
	    if(iterCounter==1)
	      {
		this->update_derivatives();
	      };
	  }
      }
    my_rho = (1.0-my_g)*my_rho_s+my_g*my_rho_l;
    //h = cp*T_star + f*L;
    my_f = my_rho_l*my_g/my_rho;
    //k = k_s*rho*(1.0-f)/rho_s + k_l*f*rho/rho_l;
    my_gm = my_g;
    my_T_m = my_T;
  };


  /***************************************************************************\
   *                      DYNAMIC ALGORITHM METHODS                          */

  bool dynamicChoice()
  {
    if(dynamicAlgorithm)
      {
	if(this->skip() && my_algorithm=="apparentCp")
	  {
	    my_T = my_T_m;
	    this->setAlgorithm("source");
	    return true;
	  }
	else
	  {
	    return false;
	  }
      }
    else
      {
	return false;
      }
  }

  bool skip()
  {
    // does the temperature at the new time skip the latent heat peak?
    forAll(my_T,i)
      {
	scalar& Ti = my_T[i];
	scalar& Toi = my_T.oldTime()[i];
	scalar& Ts = my_Ts.value();
	scalar& Tl = my_Tl.value();
		  
	if( ((Toi<Ts) && (Ti>Tl))
	    ||
	    ((Toi>Tl) && (Ti<Ts))
	    )
	  {
	    if(my_algorithm=="apparentCp")
	      {
		Info<< "\nSKIPPING phase change interval!!!\n"
		    << "adopting dynamic algorithm adjustment";
	      }
	    return true;
	  }
      }
	      
    return false;
  }

  void setAlgorithm(const word & algorithm)
  {
    Info<< "\nDynamic algorithm adjustment:\n"
      << "Changing from "
      << my_algorithm
      << " to ";
    my_algorithm = algorithm;
    Info<< my_algorithm;
  }
  
  
};
  
#endif

// sumFieldError =
// sumMag
// (
// effe_n.internalField()
// - effe_m.internalField()
// );

// measure = sumMag(effe.internalField()) + SMALL;

// residual_f = sumFieldError/measure;

