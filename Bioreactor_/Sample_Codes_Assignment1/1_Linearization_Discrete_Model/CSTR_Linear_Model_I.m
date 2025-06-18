 % <---- Main Program for Development of Linear perturbation models for Continuously Stirred Tank Reactor (CSTR) System ----> 
     
  clear all ; clc ;  close all  
   
   global CSTR_mod ;                                                                                                                          %  Global Data structure containing System related parameters 
 
   % <---- Program Initialization  ---->

    % Initialize CSTR_mod data structure   
     
   % <---- Initialization of Process parameters ---->

   CSTR_mod.Ko = 1.0e10 ;                  % ( min-1 )
   CSTR_mod.E  = 8330.1 ;                  % ( oK )
   CSTR_mod.V  = 1.0    ;                     % ( m3 )
   CSTR_mod.Cp = 1.0    ;                  % ( cal / g K )
   CSTR_mod.rho = 1e6   ;                  % ( g / m3 )
   CSTR_mod.delH = 1.3e8 ;                 % ( cal / kmol )
   
   CSTR_mod.Cpc  = 1.0   ;                 % ( cal / g K )
   CSTR_mod.rhoc = 1e6   ;                 % ( g / m3 )
   CSTR_mod.a  = 1.291e6 ;                 % ( cal / min K )
   CSTR_mod.b  = 0.5     ;                 % constant
   CSTR_mod.To  = 323.0 ;
 
   % -- Initialize steady state values of the man_ip and dist variables
   
   CSTR_mod.F   = 1.0  ;
   CSTR_mod.Cao = 2.0  ;
   CSTR_mod.Fc  = 15.0 ;
   CSTR_mod.a    = 1.678e6         ; % (cal/min)/K
   CSTR_mod.To   = 323             ; % k
   CSTR_mod.Tcin = 365             ; % K
    
   
   fprintf( '\n\n Development of Linear Perturbation Models for Continuously Stirred Tank Reactor (CSTR)' )
   
   % <---- Initialization of Inputs and States ---->
   
   n_st = 2 ;                    %  No. of states                   
   n_ip = 2 ;                     %  No. of inputs
   n_op = 1 ;                    %  No. of outputs
   n_ud = 1 ;                    %  No. of unmeasured disturbances
 
   
   Us = [ CSTR_mod.Fc CSTR_mod.F]' ;  % Steady state manipulated input vector 
   
   Ws = CSTR_mod.Cao ;                          % Steady state disturbance vector 
   
   fprintf( '\n\n CSTR: Steady State Operating Point' )
   fprintf( '\n\n Steady State Inputs \n\n\t Inlet Reactant Flow: %6.4f m^3/s and Coolant Flow: %6.4f m^3/s', CSTR_mod.F, CSTR_mod.Fc )
   fprintf( '\n\n Steady State Disturbance \n\n\t Inlet Reactant Conc.: %6.4f% mol/m3 \n\n', CSTR_mod.Cao )
   fprintf( '\n\n\t\t\t HIT ANY KEY TO CONTINUE \n'), pause
   fprintf( '\n\n Finding steadty state ......')
   
   % Solve for steady state using Newton's method  / optimization 
   
   CSTR_mod.Ca =  0.25 ;     %  Initial Guess for Steady state concenration 
   CSTR_mod.T = 390  ;       %  Initial Guess for Steady state temparature 
   
   Xs0 = [ CSTR_mod.Ca  CSTR_mod.T ]'; 
 
   % Solving nonlinear algebraic equations simultaneously using Newton's method / optimization  
   
   Xs = fsolve ('CSTR_SteadyState', Xs0);    
   
   CSTR_mod.Ca = Xs(1) ;    %  Steady state concenration 
   CSTR_mod.T  = Xs(2)  ;     %  Steady state temparature 
         
   fprintf( '\n\n Steady State Operating Point')
   fprintf( '\n\n\t Concentration: %6.4f mol/m3 and Temperature: %6.4f K', CSTR_mod.Ca, CSTR_mod.T)
   fprintf( '\n\n\t\t\t HIT ANY KEY TO CONTINUE'), pause
   
   
   Xs = [ CSTR_mod.Ca  CSTR_mod.T ]'  ;   % Steady state corresponding to specified inputs     
   Xp = Xs ;
   C_mat = [ 0 1 ]    ;                                   %  output / observation matrix when only Temperature is measured 
   Ys = C_mat * Xs ;                                   % steady state output 
   

   % Linearlized state space model using numerical diffrentiation 
   
   Z_vec = [ Xs' Us' Ws' ]' ;  
   
   Jacob_mat = Num_Jacobian('CSTR_Model_JacobFn', Z_vec );
   
    % Extract Matrices Linear state space model (continuus time)
    
    A_mat = Jacob_mat(:,1:n_st) ;              % Continuous time local linear perturbation model 
    B_mat  = Jacob_mat(:,n_st+1:n_st+n_ip);
    H_mat  = Jacob_mat(:,n_st+n_ip+1:n_st+n_ip+n_ud);      
    
    % create data structure for linear continuous time model and display
  
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
   
    
   % Discretize the continuous time model 
    
   %  Set sampling time (suggested value for the system: 0.1 min) 
    
   samp_T = input('\n\n     Select Sampling Interval (T) for model discretization : ')    ; 
   
    % generate discrete time linear state space model 
    phy = expm( samp_T * A_mat ) ;
    gama_u =  (phy - eye(n_st)) * inv(A_mat) * B_mat ; 
    gama_d = (phy - eye(n_st)) * inv(A_mat) * H_mat ;

%     % Alternative approach: Use 'c2d' function from control system tool-box  
%      [  phy, gam ] = c2d( A_mat,  Bu_mat, samp_T ) ;   
   
    fprintf( '\n\n Discrete time state space model: ')
    fprintf( '\n phy matrix =\n'), disp(phy)
    fprintf( '\n gama_u matrix =\n'), disp(gama_u)  
    fprintf( '\n gama_d matrix =\n'), disp(gama_d)  
    fprintf( '\n\n\t\t\t HIT ANY KEY TO CONTINUE \n'), pause
   

   % Create a data structure to save discrete time model 
    
    dmod_lin.phy = phy ;   dmod_lin.gama_u = gama_u ;  dmod_lin.gama_d = gama_d ;
    dmod_lin.C = C_mat ; dmod_lin.T = samp_T ; 
    dmod_lin.n_st = n_st ; dmod_lin.n_op = n_op ; 
    dmod_lin.n_ip = n_ip ; dmod_lin.n_ud = n_ud ; 
    dmod_lin.Xs = Xs ; dmod_lin.Ys = Ys ; dmod_lin.Us = Us ; dmod_lin.Ws = Ws ; 
  
   fprintf( '\n Eigenvalues of phy matrix (discrete state space model): \n\n') ;
   
   eig_phy = eig( phy ) ;
   
   disp( eig_phy ) ;
  

   % Save data of interest in a MATLAB data file 

   save CSTR_LinMod_I.mat dmod_lin    % Linear perturbation model 
   save CSTR_para.mat CSTR_mod         % CSTR parameters and operating conditions 