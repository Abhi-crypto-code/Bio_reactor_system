function [ xdot] = CSTR_Model_JacobFn( Z_vec )

global CSTR_mod

CSTR_mod.Ca = Z_vec(1);
CSTR_mod.T = Z_vec(2);
CSTR_mod.Fc = Z_vec(3);
CSTR_mod.F = Z_vec(4);
CSTR_mod.Cao = Z_vec(5);
 

K = CSTR_mod.Ko*exp(-CSTR_mod.E/CSTR_mod.T);
Gc = CSTR_mod.a/(CSTR_mod.V*CSTR_mod.rho*CSTR_mod.Cp);
B1 = CSTR_mod.a*CSTR_mod.Fc^CSTR_mod.b/(2*CSTR_mod.rhoc*CSTR_mod.Cpc);
B = CSTR_mod.Fc^(CSTR_mod.b+1)/(CSTR_mod.Fc+B1);

Q_ratio = CSTR_mod.delH/(CSTR_mod.rho*CSTR_mod.Cp);

% derivatives

dCa_by_dt = CSTR_mod.F*(CSTR_mod.Cao-CSTR_mod.Ca)/CSTR_mod.V-K*CSTR_mod.Ca;

dT_by_dt = CSTR_mod.F*(CSTR_mod.To-CSTR_mod.T)/CSTR_mod.V;
dT_by_dt = dT_by_dt - Gc*B*(CSTR_mod.T-CSTR_mod.Tcin);
dT_by_dt = dT_by_dt + Q_ratio*K*CSTR_mod.Ca;
xdot = [dCa_by_dt  dT_by_dt]';

end

