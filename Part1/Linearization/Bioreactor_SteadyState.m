function [xdot] = Bioreactor_SteadyState(Xs0)
    global BioReactor_mod
    
    BioReactor_mod.X = Xs0(1);
    BioReactor_mod.S = Xs0(2);
    BioReactor_mod.P = Xs0(3);
    
    u_SP = (BioReactor_mod.u_m)*(1- (BioReactor_mod.P/BioReactor_mod.P_m))*(BioReactor_mod.S)/(BioReactor_mod.K_m + BioReactor_mod.S + ((BioReactor_mod.S)^2)/BioReactor_mod.K_i);
    
    dX_dt = -BioReactor_mod.D * BioReactor_mod.X + u_SP*BioReactor_mod.X;
    dS_dt = BioReactor_mod.D*(BioReactor_mod.S_f - BioReactor_mod.S) * u_SP * BioReactor_mod.X / BioReactor_mod.Yxs;
    dP_dt = -BioReactor_mod.D*BioReactor_mod.P + (BioReactor_mod.alpha*u_SP + BioReactor_mod.Beta)*BioReactor_mod.X;
    
    xdot = [dX_dt dS_dt dP_dt]';
end

