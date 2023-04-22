Info<< "Updating f" << endl;

if(pimple.finalIter()) //if we are at the last iteration, the relaxation factor is set to 0.
  {
    f = 1.0/(1.0 + exp(-T_star/dT + offset));
    f = max( scalar(0) , min(f, scalar(1)) ) ; // limiting 0 < f < 1
  }
 else
   {
     f = f + relax_f*dfdT*(T_star-T_starm);
     f = max( scalar(0) , min(f, scalar(1)) ) ;
   }
