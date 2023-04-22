Info<< "updating thermal properties" << endl;

T_star = T - Tm;

cp = cp_s*(1.0-f)+cp_l*f;

h = cp*T_star +f*L;

rho = 1.0/((1.0-f)/rho_s + f/rho_l);

rho_cp = rho*cp;
k = k_s*rho*(1.0-f)/rho_s + k_l*f*rho/rho_l;
