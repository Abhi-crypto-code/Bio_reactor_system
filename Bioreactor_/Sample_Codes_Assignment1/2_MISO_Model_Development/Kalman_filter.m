clear all ; clc ;  close all  

global BioReactor_mod;

load Bioreactor_para
load Bioreactor_Linear_Model_I


n_st = dmod_lin.n_st ; n_op = dmod_lin.n_op  ;
n_ip = dmod_lin.n_ip ; n_ud = dmod_lin.n_ud ;
Xs = dmod_lin.Xs ; Ys =  dmod_lin.Ys  ;
Us = dmod_lin.Us  ; Ws  = dmod_lin.Ws  ; 
phy = dmod_lin.phy ; gama_u = dmod_lin.gama_u ;
gama_d = dmod_lin.gama_d ; C_mat = dmod_lin.C ;


samp_T = dmod_lin.T ;    
N_samples = 501 ;    

mod_type = input('Plant Simulation using (0) : Linear (1) : Nonlinear model? ')  ;


xk = zeros(n_st, N_samples) ;                  
uk = zeros(n_ip, N_samples) ;                 
yk = zeros(1, N_samples) ;        

state_sigma = (0.1)' ;                   
wk =  state_sigma * randn(n_ud, N_samples) ;
meas_sigma  = (0.1)';                       
vk = meas_sigma * randn(1, N_samples) ;        


C_mat_single = C_mat(1,:); 
yk(1) = C_mat_single * xk(:,1) + vk(1); 
      
ip1 = idinput( N_samples, 'rbs', [0 0.5] ) ;
ip2 = 0.1 *idinput( N_samples,'rbs', [0 0.5] ) ;
uk = [ ip1' ; ip2' ]  ; 

 
Xk_abs = zeros(n_st, N_samples) ;           
Uk_abs = zeros(n_ip, N_samples) ;             
Yk_abs = zeros(1, N_samples) ;    
Wk_abs = zeros(n_ud, N_samples) ;

Xk_abs(:,1) = Xs + xk(:,1) ;                 
Yk_abs(1) = C_mat_single * Xs + yk(1) ;     
Uk_abs(:,1) = Us + uk(:,1) ; 
Wk_abs(1) = Ws + wk(1) ;

kT = zeros(N_samples,1) ;
kT(1) = 0 * samp_T ; 


%kalman filter

kT = zeros(N_samples,1) ;  

kT(1) = 0 * samp_T ;  


R  = diag(meas_sigma(1)^2);                                     
Q = diag((state_sigma(1)^2));


xk_pred = zeros(n_st, N_samples);
xkhat(:,1) = zeros(n_st, 1);

Xkhat = zeros(n_st, N_samples);


xkhat(:,1) = [0.5;0;0.9];  

if size(Xs, 1) ~= size(xkhat, 1)
    error('Dimensions of Xs and xkhat must match. Check the loaded data.');
end

Xkhat(:,1) = xkhat(:,1) + Xs;      

Pk = 5*Q;

res_matrix = zeros(N_samples, 10);
res_matrix(1,:) = [1, Xk_abs(1,1), Xk_abs(2,1),Xk_abs(3,1), Xkhat(1,1) Xkhat(2,1),Xkhat(3,1), Yk_abs(1,1), Uk_abs(1,1), Uk_abs(2,1)];


% luenberger observer

%poles
desired_poles = eig(phy);

%gain
Lp = place(phy',C_mat_single',desired_poles)';

%inital val
xhat_obs = zeros(n_st,N_samples);
xhat_obs(:,1) = xkhat(:,1);





for k = 2 : N_samples 
    kT(k) = (k-1) * samp_T ; 
       
    % Plant simulation
    if (mod_type == 0)     % Linear model
        xk(:,k) = phy * xk(:,k-1) + gama_u * uk(:,k-1) + gama_d * wk(:,k-1);
        yk(k) = C_mat_single * xk(:,k) + vk(k);
        Xk_abs(:,k) = Xs + xk(:,k);
        Yk_abs(k) = Ys(1) + yk(k);
    else                   % Nonlinear model
        [t, Xt] = ode113('BioReactor_Dynamics', [0 samp_T], Xk_abs(:,k-1));
        Xk_abs(:,k) = Xt(length(t),:)';
        xk(:,k) = Xk_abs(:,k) - Xs;  % Compute deviation state first
        yk(k) = C_mat_single * xk(:,k) + vk(k);  % Use deviation state
        Yk_abs(k) = Ys(1) + yk(k);
    end

    % Inputs
    Uk_abs(:,k) = Us + uk(:,k);
    Wk_abs(:,k) = Ws + wk(:,k);

    %---------------------------- Kalman Filter----------------------------
    % Prediction step
    xk_pred(:,k) = phy * xkhat(:,k-1) + gama_u * uk(:,k-1);
    yk_pred = C_mat_single * xk_pred(:,k);
    Pk = phy * Pk * phy' + Q;

    % Kalman Gain
    Vk = R + C_mat_single * Pk * C_mat_single';
    Lk = Pk * C_mat_single' / Vk;

    % Update step
    xkhat(:,k) = xk_pred(:,k) + Lk * (yk(k) - yk_pred);
    Pk = (eye(n_st) - Lk * C_mat_single) * Pk;

    % Absolute state estimate
    Xkhat(:,k) = xkhat(:,k) + Xs;

    % Store results
    res_matrix(k,:) = [k, Xk_abs(1,k), Xk_abs(2,k), Xk_abs(3,k), Xkhat(1,k), Xkhat(2,k),Xkhat(3,k), Yk_abs(k), Uk_abs(1,k), Uk_abs(2,k)];

    %---------------- Luenberger ------------------------------------------
    xhat_obs(:,k) = phy*xhat_obs(:,k-1) + gama_u*uk(:,k-1) +  Lp*(yk(k) - C_mat_single*xhat_obs(:,k-1));

end

figure(1)
subplot(211)
plot(res_matrix(:,1), res_matrix(:,2), 'b', res_matrix(:,1), res_matrix(:,5), 'r--')
ylabel('X')
legend('True State', 'Estimated State')

subplot(212)
plot(res_matrix(:,1), res_matrix(:,8), 'k', res_matrix(:,1), res_matrix(:,2), 'b--', res_matrix(:,1), res_matrix(:,5), 'r.')
xlabel('Sampling Instant')
ylabel('S')
legend('Measured State', 'True State', 'Estimated State')


figure(1)
subplot(211)
plot(res_matrix(:,1), res_matrix(:,2), 'b', res_matrix(:,1), res_matrix(:,5), 'r--')
ylabel('X')
legend('True State', 'Estimated State')
figure(2)

plot(res_matrix(:,1), res_matrix(:,2) - res_matrix(:,5), 'b')
ylabel('X Error')
title('State Estimation Error')

figure(3)
plot(res_matrix(:,1), res_matrix(:,3) - res_matrix(:,6), 'b')
xlabel('Sampling Instant')
ylabel('S Error')


figure(4)
plot(res_matrix(:,1), res_matrix(:,4) - res_matrix(:,7), 'b')
xlabel('Sampling Instant')
ylabel('P Error')
   
% Calculate estimation errors
esterr1 = res_matrix(2:end,2) - res_matrix(2:end,5);
esterr2 = res_matrix(2:end,3) - res_matrix(2:end,6);
esterr3 = res_matrix(2:end,4) - res_matrix(2:end,7);


SSE1 = esterr1' * esterr1;
SSE2 = esterr2' * esterr2;
SSE3 = esterr3' * esterr3;
SSEa = [SSE1 SSE2 SSE3]; 
   
fprintf('SSE for X: %f\n', SSEa(1));
fprintf('SSE for S: %f\n', SSEa(2));
fprintf('SSE for P : %f\n',SSEa(3));


%---luenberger
figure;
plot(kT, xk(1,:), 'b',kT, xhat_obs(1,:), 'r--');
legend('True value','Luenberger Observer');


save result_Kalman_Filter

