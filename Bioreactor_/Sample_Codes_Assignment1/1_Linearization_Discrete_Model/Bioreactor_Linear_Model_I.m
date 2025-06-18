



clear all ; clc ;  close all  



%--------------------------------------------------------------Part 1------------------------------------------------------------------------------------------



global BioReactor_mod; 
Init_Graphics_Style ; 
  BioReactor_mod.Yxs = 0.4; %0.4 g/L
  BioReactor_mod.Beta = 0.2; %0.2 h-1
  BioReactor_mod.P_m = 50; %50 g/L
  BioReactor_mod.K_i = 22; % 22 g/L
  BioReactor_mod.alpha = 2.2; % 2.2 g/L
  BioReactor_mod.u_m = 0.48; %0.48 h-1
  BioReactor_mod.K_m = 1.2; % 1.2 g/L
  
  BioReactor_mod.S_f = 20; % 20 g/L
  BioReactor_mod.D = 0.202; % 0.202
 %  BioReactor_mod.Disturbance = 2.5397;
   
   fprintf( '\n\n Development of Linear Perturbation Models for Bioreactor system' )
   
   n_st = 3;                    %  No. of states                   
   n_ip =  2;                     %  No. of inputs
   n_op = 3;                    %  No. of outputs
   n_ud = 1;                    %  No. of unmeasured disturbances
 
   
   Us = [ BioReactor_mod.D BioReactor_mod.S_f]' ; % manipulated variables.
   
    Ws = 2.5397;                   %  disturbance
   
   fprintf( '\n\n Bioreactor: Steady State Operating Point' )
   fprintf( '\n\n Steady State Inputs \n\n\t Dilution Rate : %6.4f m^3/s and Feed Substrate Concentration: %6.4f m^3/s', BioReactor_mod.D, BioReactor_mod.S_f )
   fprintf( '\n\n Steady State Disturbance \n\n\t Disturbance.: %6.4f% mol/m3 \n\n', 0 )
   fprintf( '\n\n\t\t\t HIT ANY KEY TO CONTINUE \n'), pause
   fprintf( '\n\n Finding steadty state ......')
   
  
   
   BioReactor_mod.X =  0.5;%0.25 
   BioReactor_mod.S = 20;       
   BioReactor_mod.P = 100;
   Xs0 = [ BioReactor_mod.X  BioReactor_mod.S BioReactor_mod.P]'; % output 
 
options = optimoptions('fsolve', 'Display', 'iter', 'TolFun', 1e-10, 'TolX', 1e-10, 'MaxIter', 1000);

Xs = fsolve(@Bioreactor_SteadyState, Xs0, options); % to find steady state
   
   
   BioReactor_mod.X = Xs(1) ;
   BioReactor_mod.S  = Xs(2) ;
   BioReactor_mod.P = Xs(3);



         
   fprintf( '\n\n Steady State Operating Point')
   fprintf( '\n\n\t X: %6.4f  S: %6.4f  P: %6.4f', BioReactor_mod.X, BioReactor_mod.S,BioReactor_mod.P) % at steady state operating points
   fprintf( '\n\n\t\t\t HIT ANY KEY TO CONTINUE'), pause
   
   
      Xs = [BioReactor_mod.X  BioReactor_mod.S BioReactor_mod.P]'  ;      
      Xp = Xs ;
     C_mat = [1 0 0; 0 1 0 ; 0 0 1];                                 
     Ys = C_mat .* Xs ;                           
   
   Z_vec = [ Xs' Us' Ws' ]' ;  
   
   Jacob_mat = Num_Jacobian(@Bioreactor_Model_JacobFn, Z_vec);
   
   
    
    A_mat = Jacob_mat(:,1:n_st) ;
    B_mat  = Jacob_mat(:,n_st+1:n_st+n_ip);
    H_mat  = Jacob_mat(:,n_st+n_ip+1:n_st+n_ip+n_ud);      
    
    
    fprintf( '\n\n Continuous time state space model: ')
    fprintf( '\n A matrix =\n'), disp(A_mat)
    fprintf( '\n B matrix =\n'), disp(B_mat)
    fprintf( '\n H matrix =\n'), disp(H_mat) 
    fprintf( '\n C matrix =\n'), disp(C_mat)
    fprintf( '\n\t\t\t HIT ANY KEY TO CONTINUE \n'), pause
    
   fprintf( '\n\n Eigenvalues of A matrix (continuous time state space model): \n\n')
   eig_A = eig(A_mat);
   disp(eig_A) ;
   fprintf( '\n\n\t\t\t HIT ANY KEY TO CONTINUE \n'), pause
   
    
   samp_T = input('\n\n     Select Sampling Interval (T) for model discretization : '); 
   
   

    phy = expm( samp_T .* A_mat ) ;
    gama_u =  (phy - eye(n_st)) * inv(A_mat) * B_mat ; 
    gama_d = (phy - eye(n_st)) * inv(A_mat) * H_mat ;
 
     [  phy, gam ] = c2d( A_mat,  B_mat, samp_T ) ;   
   
    fprintf( '\n\n Discrete time state space model: ')
    fprintf( '\n phy matrix =\n'), disp(phy)
    fprintf( '\n gama_u matrix =\n'), disp(gama_u)  
    fprintf( '\n gama_d matrix =\n'), disp(gama_d)  
    fprintf( '\n\n\t\t\t HIT ANY KEY TO CONTINUE \n'), pause
    
    
    dmod_lin.phy = phy ;   dmod_lin.gama_u = gama_u ;  dmod_lin.gama_d = gama_d ;
    dmod_lin.C = C_mat ; dmod_lin.T = samp_T ; 
    dmod_lin.n_st = n_st ; dmod_lin.n_op = n_op ; 
    dmod_lin.n_ip = n_ip ; dmod_lin.n_ud = n_ud ; 
    dmod_lin.Xs = Xs ; dmod_lin.Ys = Ys ; dmod_lin.Us = Us ; dmod_lin.Ws = Ws ; 
  
   fprintf( '\n Eigenvalues of phy matrix (discrete state space model): \n\n') ;
   
   eig_phy = eig( phy ) ;
   
   disp( eig_phy ) ;
  


   save Bioreactor_Linear_Model_I.mat dmod_lin    
   save Bioreactor_para.mat BioReactor_mod