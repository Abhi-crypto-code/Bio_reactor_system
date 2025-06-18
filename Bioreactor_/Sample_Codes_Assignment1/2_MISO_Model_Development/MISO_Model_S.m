% <---- Main Program for Development 2,nd order ARMAX model 
%      model for Continuously Stirred Tank Reactor (CSTR) System  ---> 
     
clear all ; clc ;  close all  

global BioReactor_mod;        %  Global Data structure containing System related parameters 

load Bioreactor_para    % Initialize CSTR_mod data structure and operating conditions
load Bioreactor_Linear_Model_I  % Load discrete linear model obtained through linearization 

% Following local variables are created only for improving readability of the program 
n_st = dmod_lin.n_st ; n_op = dmod_lin.n_op  ;
n_ip = dmod_lin.n_ip ; n_ud = dmod_lin.n_ud ;
Xs = dmod_lin.Xs ; Ys =  dmod_lin.Ys  ;    % Steady state operating conditions 
Us = dmod_lin.Us  ; Ws  = dmod_lin.Ws  ; 
phy = dmod_lin.phy ; gama_u = dmod_lin.gama_u ;
gama_d = dmod_lin.gama_d ; C_mat = dmod_lin.C ;

% Note: It is possible to work directly with elements of dmod_lin object
% without requiring creation of above local variables 

samp_T = dmod_lin.T ;       % Sampling interval  
N_samples = 501 ;   % Number of samples in open loop simulation run     

% <---- Program Initialization  ---->

fprintf( '\n\n Development of ARMAX for Bioreactor' )
fprintf( '\n\n\t'), 
mod_type = input('Plant Simulation using (0) : Linear (1) : Nonlinear model? ')  ;

%----Initialization for absolute and dev state variables -------

% Create dummy arrays for simulation 
% k'th column of these arrays corresponds to vector at k'th sampling instant 

xk = zeros(n_st, N_samples) ;     % Matrices for saving deviation variables             
uk = zeros(n_ip, N_samples) ;                 
yk = zeros(1, N_samples) ;        % CHANGE: Using only one output

state_sigma = [0.001]' ;                   % Generate state noise sequence for simulation 
wk =  state_sigma * randn(n_ud, N_samples) ;
meas_sigma  = [0.001]';                      % Generate state noise sequence for simulation 
vk = meas_sigma * randn(1, N_samples) ;        

% Use only the first row of C_mat if you want only the first output
C_mat_single = C_mat(2,:);        % CHANGE: Using only first row of output matrix
yk(1) = C_mat_single * xk(:,1) + vk(1);  % Initial dev. Measurement (at k = 0)  

%    Generation of Random Binary Input Sequences for System  Identification
      
ip1 = idinput( N_samples, 'rbs', [0 0.005] ) ;
ip2 = 0.1 *idinput( N_samples,'rbs', [0 0.005] ) ;
uk = [ ip1' ; ip2' ]  ; 

%  Matrices for saving Absolute variables 
Xk_abs = zeros(n_st, N_samples) ;           
Uk_abs = zeros(n_ip, N_samples) ;             
Yk_abs = zeros(1, N_samples) ;      % CHANGE: Using only one output
Wk_abs = zeros(n_ud, N_samples) ;

Xk_abs(:,1) = Xs + xk(:,1) ;                 % Plant: Initial abs. state (at k = 0)
Yk_abs(1) = C_mat_single * Xs + yk(1) ;     % CHANGE: Using only first output
Uk_abs(:,1) = Us + uk(:,1) ; 
Wk_abs(1) = Ws + wk(1) ;

% <------------- Generate the plant data: Open Loop Dynamic Simulation ------------------>

kT = zeros(N_samples,1) ;
kT(1) = 0 * samp_T ;   %  k = 1 corresponds to time = 0 

for k = 2 : N_samples 
                                            % Print sampling time on screen 
    kT(k) = (k-1) * samp_T ; 


       
    %   <---- Plant simulation form instnat (k-1) to (k) ----> 
    %   Integrate the model equations for U(k-1) and W(k-1) to compute X(k) and Y(k) ---->
   
    if ( mod_type == 0 )     % Process simulation using discrete linear perturbation model
        xk(:,k) = phy * xk(:,k-1) + gama_u * uk(:,k-1) + gama_d * wk(:,k-1) ;
        yk(k) = C_mat_single * xk(:,k) + vk(k) ;   % CHANGE: Using only first output

        Xk_abs(:,k) = Xs + xk(:,k) ;              % Save simulation data in absolute variables 
        Yk_abs(k) = Ys(2) + yk(k) ;               % CHANGE: Using only first element of Ys
      
    else                             % Process simulation using nonlinear ODE model 
      
       
        
        %  Runge Kutta integration over interval [(k-1)T, kT] using MATLAB ODE solver 
        [t,Xt] = ode15s( 'BioReactor_Dynamics', [0 samp_T] , Xk_abs(:,k-1) ) ;

        %odefun = @(t, x) BioReactor_Dynamics(t, x, BioReactor_mod);
        %[t, Xt] = ode15s(odefun, [0 samp_T], Xk_abs(:,k-1));

        Xk_abs(:,k) = Xt( length(t),:)' ;         % State at instant (k+1)  
            
        yk(k) = C_mat_single * Xk_abs(:,k) + vk(k);   % CHANGE: Using only first output
        Yk_abs(k) = yk(k) + Ys(2) ;                    % CHANGE: Using only first element of Ys
        xk(:,k) = Xk_abs(:,k) - Xs ;                  % Generate Perturbation variables 
    end

    % <-----Specify Inputs at k'th sampling instant ------>    
    Uk_abs(:,k) = Us + uk(:,k) ; 
    Wk_abs(k) = Ws + wk(k) ;
end


%  <---- Display simulation results graphically ----> 

Init_Graphics_Style ;     % Set parameters for graphics (Optional) 

%figure(1),subplot(211), plot( kT, xk(1,:) ) ;   
%xlabel('Sampling Instant'), ylabel('Conc.(mod/m3)'), title( 'State Perturbations') ;
%figure(1),subplot(212), plot( kT , xk(2,:) ) ;
%xlabel('Sampling Instant'), ylabel('Temp.(K)') ;
   
%figure(2),subplot(211), stairs( kT , uk(1,:) ) ;   
%xlabel('Sampling Instant'), ylabel('Coolent Flow'), title( 'Man. Input Perturbations') ;
%figure(2),subplot(212), stairs(kT , uk(2,:) ) ;
%xlabel('Sampling Instant'), ylabel('Inflow')
   
%figure(3),stairs( kT, wk ) ;   
%xlabel('Sampling Instant'), ylabel('Inlet Conc. (mol/m3)'), title( 'Unmeas. Dist. Perturbations') ;
   
%figure(4),stairs( kT, yk ) ;   
%xlabel('Sampling Instant'), ylabel('Output'), title( 'Output Perturbations') ;
   
% <---- MISO ARMAX model identification -----> 
   
% split data set into identification and validation sets  
%Bioreactor_id_data = iddata(yk(1:400)', uk(:,1:400)', samp_T);      % Training data
%Bioreactor_val_data = iddata(yk(400:501)', uk(:,400:501)', samp_T); % Testing data
 

% Add 1% disturbance
disturbance_level = 0.05;  % 1% disturbance

% Apply disturbance to outputs (yk) and/or inputs (uk) as needed
yk_disturbed = yk + disturbance_level * randn(size(yk)) .* yk; 
uk_disturbed = uk + disturbance_level * randn(size(uk)) .* uk;

% Create IDDATA objects with disturbed data
Bioreactor_id_data = iddata(yk_disturbed(1:400)', uk_disturbed(:,1:400)', samp_T); % Training data
Bioreactor_val_data = iddata(yk_disturbed(1:501)', uk_disturbed(:,1:501)', samp_T); % Testing data


















% Third order model
na = 3;            % order of A polynomial 
nb = [3 3];        % order of B polynomials w.r.t. man. inputs
nc = 3;            % order of C polynomial 
nd = [1 1];        % time delays w.r.t. man. inputs  
Bioreactor_armax_2 = armax(Bioreactor_id_data, [na nb nc nd]);

% Fourth order model
na = 4;            % order of A polynomial 
nb = [4 4];        % order of B polynomials w.r.t. man. inputs
nc = 4;            % order of C polynomial 
nd = [1 1];        % time delays w.r.t. man. inputs  
Bioreactor_armax_3 = armax(Bioreactor_id_data, [na nb nc nd]);

% Analysis of model residuals / innovations {e(k)} (3rd order model)
ek_obj = pe(Bioreactor_armax_2, Bioreactor_id_data);  % FIXED: Variable name corrected
ek = get(ek_obj, 'OutputData');
%figure(5), plot(ek) 
%title('Model Residuals / Innovations'), xlabel('Sampling Instant'), ylabel('e(k)')
   
% Auto-correlation and cross correlation analysis 
%figure(6), cra([ek ek], 20, 0, 1) 
%figure(7), subplot(211), cra([ek uk(1,1:400)'], 20, 0, 1)
%figure(7), subplot(212), cra([ek uk(2,1:400)'], 20, 0, 1)

% Model validation
figure(8), compare(Bioreactor_val_data, Bioreactor_armax_2, Bioreactor_armax_3)  
   
% Generate state realization of the identified models
dmod_id_2 = ss(Bioreactor_armax_2) ;   % FIXED: Variable name corrected
dmod_id_3 = ss(Bioreactor_armax_3)  ;  % FIXED: Variable name corrected
   
save Bioreactor_IdMod.mat Bioreactor_armax_2 Bioreactor_armax_3 dmod_id_2 dmod_id_3  % Save identification results

% <---- Additional System Identification Models for Bioreactor ---->
% This code implements OE, ARX, and BJ models similar to the ARMAX model

% The following code assumes that the data preparation is already done
% by the main script (with single output as requested)

% <---- ARX Model Identification ---->
fprintf('\n\n ARX Model Identification for Bioreactor \n');

% Third order ARX model
na_arx = 3;            % order of A polynomial 
nb_arx = [3 3];        % order of B polynomials w.r.t. man. inputs
nk_arx = [1 1];        % time delays w.r.t. man. inputs  
Bioreactor_arx_3 = arx(Bioreactor_id_data, [na_arx nb_arx nk_arx]);

% Fourth order ARX model
na_arx = 4;            % order of A polynomial 
nb_arx = [4 4];        % order of B polynomials w.r.t. man. inputs
nk_arx = [1 1];        % time delays w.r.t. man. inputs  
Bioreactor_arx_4 = arx(Bioreactor_id_data, [na_arx nb_arx nk_arx]);

% Model validation for ARX models
figure(9), compare(Bioreactor_val_data, Bioreactor_arx_3, Bioreactor_arx_4);
title('ARX Model Validation');

% Analysis of ARX model residuals
ek_arx = pe(Bioreactor_arx_3, Bioreactor_id_data);
ek_arx_data = get(ek_arx, 'OutputData');
figure(10), plot(ek_arx_data);
title('ARX Model Residuals'), xlabel('Sampling Instant'), ylabel('e(k)');

% <---- Output-Error (OE) Model Identification ---->
fprintf('\n\n OE Model Identification for Bioreactor \n');

% Third order OE model
nb_oe = [3 3];         % order of B polynomials
nf_oe = [3 3];         % order of F polynomials
nk_oe = [1 1];         % time delays
Bioreactor_oe_3 = oe(Bioreactor_id_data, [nb_oe nf_oe nk_oe]);

% Fourth order OE model
nb_oe = [4 4];         % order of B polynomials
nf_oe = [4 4];         % order of F polynomials
nk_oe = [1 1];         % time delays
Bioreactor_oe_4 = oe(Bioreactor_id_data, [nb_oe nf_oe nk_oe]);

% Model validation for OE models
figure(11), compare(Bioreactor_val_data, Bioreactor_oe_3, Bioreactor_oe_4);
title('OE Model Validation');

% Analysis of OE model residuals
ek_oe = pe(Bioreactor_oe_3, Bioreactor_id_data);
ek_oe_data = get(ek_oe, 'OutputData');
figure(12), plot(ek_oe_data);
title('OE Model Residuals'), xlabel('Sampling Instant'), ylabel('e(k)');

% <---- Box-Jenkins (BJ) Model Identification ---->
fprintf('\n\n BJ Model Identification for Bioreactor \n');

% Third order BJ model
nb_bj = [3 3];         % order of B polynomials 
nc_bj = 3;             % order of C polynomial
nd_bj = 3;             % order of D polynomial
nf_bj = [3 3];         % order of F polynomials
nk_bj = [1 1];         % time delays
Bioreactor_bj_3 = bj(Bioreactor_id_data, [nb_bj nc_bj nd_bj nf_bj nk_bj]);

% Fourth order BJ model
nb_bj = [4 4];         % order of B polynomials
nc_bj = 4;             % order of C polynomial
nd_bj = 4;             % order of D polynomial
nf_bj = [4 4];         % order of F polynomials
nk_bj = [1 1];         % time delays
Bioreactor_bj_4 = bj(Bioreactor_id_data, [nb_bj nc_bj nd_bj nf_bj nk_bj]);

% Model validation for BJ models
figure(13), compare(Bioreactor_val_data, Bioreactor_bj_3, Bioreactor_bj_4);
title('BJ Model Validation');

% Analysis of BJ model residuals
ek_bj = pe(Bioreactor_bj_3, Bioreactor_id_data);
ek_bj_data = get(ek_bj, 'OutputData');
figure(14), plot(ek_bj_data);
title('BJ Model Residuals'), xlabel('Sampling Instant'), ylabel('e(k)');

% <---- Model Comparison Across Different Types ---->
fprintf('\n\n Comparing Different Model Types \n');

% Compare the best models of each type
figure(15), compare(Bioreactor_val_data, Bioreactor_arx_3, Bioreactor_oe_3, Bioreactor_bj_3, Bioreactor_armax_2);
title('Comparison of Different Model Types (3rd Order)');

% Generate state space realizations
dmod_arx = ss(Bioreactor_arx_3);
dmod_oe = ss(Bioreactor_oe_3);
dmod_bj = ss(Bioreactor_bj_3);

% Save all models
save Bioreactor_AllModels.mat Bioreactor_arx_3 Bioreactor_arx_4 Bioreactor_oe_3 Bioreactor_oe_4 Bioreactor_bj_3 Bioreactor_bj_4 dmod_arx dmod_oe dmod_bj

% <---- Calculate Fit Metrics For Model Comparison ---->
fprintf('\n\n Model Fit Percentages: \n');

% Calculate and display fit percentages
fit_arx_3 = compare(Bioreactor_val_data, Bioreactor_arx_3);
fit_arx_4 = compare(Bioreactor_val_data, Bioreactor_arx_4);
fit_oe_3 = compare(Bioreactor_val_data, Bioreactor_oe_3);
fit_oe_4 = compare(Bioreactor_val_data, Bioreactor_oe_4);
fit_bj_3 = compare(Bioreactor_val_data, Bioreactor_bj_3);
fit_bj_4 = compare(Bioreactor_val_data, Bioreactor_bj_4);
fit_armax_2 = compare(Bioreactor_val_data, Bioreactor_armax_2);
fit_armax_3 = compare(Bioreactor_val_data, Bioreactor_armax_3);










