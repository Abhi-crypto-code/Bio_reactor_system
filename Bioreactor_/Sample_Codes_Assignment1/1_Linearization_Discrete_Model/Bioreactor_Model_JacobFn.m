function [xdot] = Bioreactor_Model_JacobFn(Z_vec)
    global BioReactor_mod
    
    BioReactor_mod.X = Z_vec(1);
    BioReactor_mod.S = Z_vec(2); 
    BioReactor_mod.P = Z_vec(3);  
    BioReactor_mod.D = Z_vec(4);
    BioReactor_mod.S_f = Z_vec(5);
    BioReactor_mod.Disturbance = Z_vec(6);
    
    u_SP = (BioReactor_mod.u_m)*(1- (BioReactor_mod.P/BioReactor_mod.P_m))*(BioReactor_mod.S)/(BioReactor_mod.K_m + BioReactor_mod.S + ((BioReactor_mod.S)^2)/BioReactor_mod.K_i);
    
    dX_dt = -BioReactor_mod.D * BioReactor_mod.X + u_SP*BioReactor_mod.X;
    dS_dt = BioReactor_mod.D*(BioReactor_mod.S_f - BioReactor_mod.S) * u_SP * BioReactor_mod.X / BioReactor_mod.Yxs + + BioReactor_mod.Disturbance;
    dP_dt = -BioReactor_mod.D*BioReactor_mod.P + (BioReactor_mod.alpha*u_SP + BioReactor_mod.Beta)*BioReactor_mod.X;
    
    xdot = [dX_dt dS_dt dP_dt]';
end