Info<< "updating thermal properties" << endl;

T = Tm + T_star;

//f = (rhoh/rho_s - cp_s*T_star)/(rhoh*(1.0/rho_s - 1.0/rho_l) + (cp_l - cp_s)*T_star + L);

//#include "calculate_dfdT.H"

dfdT = f*(1.0-f)/dT;

rho = 1.0/((1.0-f)/rho_s + f/rho_l);
cp = cp_s*(1.0-f)+cp_l*f;
h = cp*T_star + f*L;
d_rhoh_dT = rho*cp + rho*rho*(1.0/rho_s - 1.0/rho_l)*dfdT*h + rho*((cp_l - cp_s)*T_star + L)*dfdT;

rhoh = rho*h;

k = k_s*rho*(1.0-f)/rho_s + k_l*f*rho/rho_l;
