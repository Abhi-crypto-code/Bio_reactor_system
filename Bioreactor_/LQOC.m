%% 
clear all; clc; close all  

global BioReactor_mod;

load Bioreactor_para
load Bioreactor_Linear_Model_I

n_st = dmod_lin.n_st; n_op = dmod_lin.n_op;
n_ip = dmod_lin.n_ip; n_ud = dmod_lin.n_ud;
Xs = dmod_lin.Xs; Ys = dmod_lin.Ys;
Us = dmod_lin.Us; Ws = dmod_lin.Ws; 
phy = dmod_lin.phy; gama_u = dmod_lin.gama_u;
gama_d = dmod_lin.gama_d; C_mat = dmod_lin.C;

samp_T = dmod_lin.T;    
N_samples = 501;    

disp('Assignment 3: Linear Quadratic Optimal Controller Implementation');
disp('1. Linear plant simulation with setpoint tracking');
disp('2. Linear plant simulation with disturbance rejection');
disp('3. Nonlinear plant simulation with setpoint tracking');
disp('4. Nonlinear plant simulation with disturbance rejection');
mod_type = input('Select simulation type (1-4): ');

% Converting input to appropriate simulation parameters
is_nonlinear = mod_type > 2;
is_setpoint_tracking = mod_type == 1 || mod_type == 3;

xk = zeros(n_st, N_samples);                  
uk = zeros(n_ip, N_samples);                 
yk = zeros(n_op, N_samples);        

state_sigma = 0.0001 * ones(n_ud, 1);                   
wk = state_sigma .* randn(n_ud, N_samples);
meas_sigma = 0.0001 * ones(n_op, 1);                       
vk = meas_sigma .* randn(n_op, N_samples);        

% Output matrix for each measured output
for i = 1:n_op
    yk(i, 1) = C_mat(i,:) * xk(:,1) + vk(i,1); 
end

% Initialize input (not used for closed-loop control but needed for simulation setup)
ip1 = 0.01 * idinput(N_samples, 'rbs', [0 0.05]);
ip2 = 0.1 * idinput(N_samples, 'rbs', [0 0.05]);
uk = [ip1'; ip2'];

% Absolute state, input, output, and disturbance
Xk_abs = zeros(n_st, N_samples);           
Uk_abs = zeros(n_ip, N_samples);             
Yk_abs = zeros(n_op, N_samples);    
Wk_abs = zeros(n_ud, N_samples);

% Initialize at steady state (k=1)
Xk_abs(:,1) = Xs;                 
for i = 1:n_op
    Yk_abs(i,1) = Ys(i) + yk(i,1);     
end
Uk_abs(:,1) = Us;
Wk_abs(:,1) = Ws + wk(:,1);

% Time array
kT = zeros(N_samples,1);
kT(1) = 0;

% Kalman filter setup
R = diag(meas_sigma.^2);                                     
Q = diag(state_sigma.^2);

xk_pred = zeros(n_st, N_samples);
xkhat = zeros(n_st, N_samples);
% Initial state estimate (can be different from true initial state to test convergence)
xkhat(:,1) = [0.5; 0; 0.9];  

Xkhat = zeros(n_st, N_samples);
Xkhat(:,1) = xkhat(:,1) + Xs;      

Pk = 5*Q;

% Result matrix for storing simulation data
res_matrix = zeros(N_samples, 13);

% LQOC design with DARE
% Choose weights for states and inputs based on their magnitudes
% Trial and error procedure as mentioned in assignment
state_magnitudes = abs(Xs);
input_magnitudes = abs(Us);

% Make weights inversely proportional to state magnitudes to normalize effects
Wx_diag = 1./(state_magnitudes + 0.1).^2;  % Adding 0.1 to avoid division by zero
Wx = diag(Wx_diag);

% Input weights - can be adjusted based on control objectives
Wu_diag = 1./(input_magnitudes + 0.1).^2;
Wu = diag(Wu_diag);

% Solve Discrete Algebraic Riccati Equation
[X_mat, L_vec, G_inf] = dare(phy, gama_u, Wx, Wu);

% For setpoint tracking, define setpoint changes
if is_setpoint_tracking
    % Create step changes in setpoint at different times
    Xsetpoint = repmat(Xs, 1, N_samples);
    
    % Step change in X at sample 100
    if N_samples > 100
        Xsetpoint(1, 100:end) = Xs(1) * 1.1;  % 10% increase in X
    end
    
    % Step change in S at sample 300
    if N_samples > 300
        Xsetpoint(2, 300:end) = Xs(2) * 0.9;  % 10% decrease in S
    end
else
    % For disturbance rejection, keep setpoint at steady state
    Xsetpoint = repmat(Xs, 1, N_samples);
end

% Initialize first row of results matrix
res_matrix(1,:) = [1, Xk_abs(1,1), Xk_abs(2,1), Xk_abs(3,1), Xkhat(1,1), Xkhat(2,1), Xkhat(3,1), Yk_abs(1,1), Uk_abs(1,1), Uk_abs(2,1), Xsetpoint(1,1), Xsetpoint(2,1), Xsetpoint(3,1)];

% Main simulation loop
for k = 2:N_samples 
    kT(k) = (k-1) * samp_T; 
       
    % Plant simulation
    if ~is_nonlinear  % Linear model
        xk(:,k) = phy * xk(:,k-1) + gama_u * uk(:,k-1) + gama_d * wk(:,k-1);
        for i = 1:n_op
            yk(i,k) = C_mat(i,:) * xk(:,k) + vk(i,k);
        end
        Xk_abs(:,k) = Xs + xk(:,k);
        for i = 1:n_op
            Yk_abs(i,k) = Ys(i) + yk(i,k);
        end
    else  % Nonlinear model
        BioReactor_mod.D = Uk_abs(1,k-1);
        BioReactor_mod.S_f = Uk_abs(2,k-1);
        BioReactor_mod.Disturbance = Wk_abs(1,k-1);
     
    [t, Xt] = ode15s(@BioReactor_Dynamics, [0 samp_T], Xk_abs(:,k-1));
        Xk_abs(:,k) = Xt(end,:)';

        xk(:,k) = Xk_abs(:,k) - Xs;
        % Add disturbance for disturbance rejection test
        if ~is_setpoint_tracking && k == 200  % Add disturbance at sample 200
            Xk_abs(:,k) = Xk_abs(:,k) + [0.01; 0.01; 0.01];  % Disturbance in all states
        end
        
        for i = 1:n_op
            Yk_abs(i,k) = C_mat(i,:) * Xk_abs(:,k) + vk(i,k);
        end
        xk(:,k) = Xk_abs(:,k) - Xs;  % Compute deviation state
        for i = 1:n_op
            yk(i,k) = Yk_abs(i,k) - Ys(i);  % Compute deviation output
        end
    end

    % Kalman Filter implementation
    % Prediction step
    xk_pred(:,k) = phy * xkhat(:,k-1) + gama_u * uk(:,k-1);
    yk_pred = C_mat * xk_pred(:,k);
    Pk = phy * Pk * phy' + Q;

    % Kalman Gain calculation
    Vk = R + C_mat * Pk * C_mat';
    Lk = Pk * C_mat' / Vk;

    % Update step using actual measurement
    xkhat(:,k) = xk_pred(:,k) + Lk * (yk(:,k) - yk_pred);
    Pk = (eye(n_st) - Lk * C_mat) * Pk;

    % Absolute state estimate
    Xkhat(:,k) = xkhat(:,k) + Xs;

    % LQOC control law
    % For setpoint tracking, calculate error
    error_k = Xsetpoint(:,k) - Xkhat(:,k);
    
    % Calculate control input considering offset
    if is_setpoint_tracking
        % For setpoint tracking
        uk(:,k) = -G_inf * xkhat(:,k) + G_inf * (Xsetpoint(:,k) - Xs);
    else
        % For disturbance rejection (regulatory control)
        uk(:,k) = -G_inf * xkhat(:,k);
    end

    % Apply input and store absolute values
    Uk_abs(:,k) = Us + uk(:,k);
    Wk_abs(:,k) = Ws + wk(:,k);

    % Store results for this time step
    res_matrix(k,:) = [k, Xk_abs(1,k), Xk_abs(2,k), Xk_abs(3,k), Xkhat(1,k), Xkhat(2,k), Xkhat(3,k), Yk_abs(1,k), Uk_abs(1,k), Uk_abs(2,k), Xsetpoint(1,k), Xsetpoint(2,k), Xsetpoint(3,k)];
end

% Plots for analysis

% Plot 1: Comparison of true state and estimated state
figure(1);
subplot(311);
plot(res_matrix(:,1), res_matrix(:,2), 'b', res_matrix(:,1), res_matrix(:,5), 'r', res_matrix(:,1), res_matrix(:,11), 'g--');
ylabel('X');
legend('True State', 'Estimated State', 'Setpoint');
title('State 1 (X) - Comparison of True, Estimated, and Setpoint');

subplot(312);
plot(res_matrix(:,1), res_matrix(:,3), 'b', res_matrix(:,1), res_matrix(:,6), 'r', res_matrix(:,1), res_matrix(:,12), 'g--');
ylabel('S');
legend('True State', 'Estimated State', 'Setpoint');
title('State 2 (S) - Comparison of True, Estimated, and Setpoint');

subplot(313);
plot(res_matrix(:,1), res_matrix(:,4), 'b', res_matrix(:,1), res_matrix(:,7), 'r', res_matrix(:,1), res_matrix(:,13), 'g--');
xlabel('Sampling Instant');
ylabel('P');
legend('True State', 'Estimated State', 'Setpoint');
title('State 3 (P) - Comparison of True, Estimated, and Setpoint');

% Plot 2: Estimation error
figure(2);
subplot(311);
plot(res_matrix(:,1), res_matrix(:,2) - res_matrix(:,5), 'b');
ylabel('X Error');
title('Estimation Error for State X');

subplot(312);
plot(res_matrix(:,1), res_matrix(:,3) - res_matrix(:,6), 'b');
ylabel('S Error');
title('Estimation Error for State S');

subplot(313);
plot(res_matrix(:,1), res_matrix(:,4) - res_matrix(:,7), 'b');
xlabel('Sampling Instant');
ylabel('P Error');
title('Estimation Error for State P');

% Plot 3: Control inputs
figure(3);
subplot(211);
plot(res_matrix(:,1), res_matrix(:,9), 'b');
ylabel('D');
title('Manipulated Variable: Dilution Rate (D)');

subplot(212);
plot(res_matrix(:,1), res_matrix(:,10), 'b');
xlabel('Sampling Instant');
ylabel('S_f');
title('Manipulated Variable: Feed Substrate Concentration (S_f)');

% Calculate performance metrics
esterr1 = res_matrix(2:end,2) - res_matrix(2:end,5);  % X estimation error
esterr2 = res_matrix(2:end,3) - res_matrix(2:end,6);  % S estimation error
esterr3 = res_matrix(2:end,4) - res_matrix(2:end,7);  % P estimation error

% Sum of squared errors for estimation
SSE1 = esterr1'*esterr1;
SSE2 = esterr2'*esterr2;
SSE3 = esterr3'*esterr3;

% Control performance - tracking error
if is_setpoint_tracking
    track_err1 = res_matrix(2:end,11) - res_matrix(2:end,2);  % X tracking error
    track_err2 = res_matrix(2:end,12) - res_matrix(2:end,3);  % S tracking error
    
    % Sum of squared errors for tracking
    SSE_track1 = track_err1'*track_err1;
    SSE_track2 = track_err2'*track_err2;
    
    disp('Setpoint Tracking Performance:');
    disp(['SSE for X tracking: ', num2str(SSE_track1)]);
    disp(['SSE for S tracking: ', num2str(SSE_track2)]);
end

disp('State Estimation Performance:');
disp(['SSE for X estimation: ', num2str(SSE1)]);
disp(['SSE for S estimation: ', num2str(SSE2)]);
disp(['SSE for P estimation: ', num2str(SSE3)]);

% Save results based on simulation type
if is_nonlinear
    if is_setpoint_tracking
        save LQOC_Nonlinear_Setpoint_Tracking
        disp('Results saved as LQOC_Nonlinear_Setpoint_Tracking');
    else
        save LQOC_Nonlinear_Disturbance_Rejection
        disp('Results saved as LQOC_Nonlinear_Disturbance_Rejection');
    end
else
    if is_setpoint_tracking
        save LQOC_Linear_Setpoint_Tracking
        disp('Results saved as LQOC_Linear_Setpoint_Tracking');
    else
        save LQOC_Linear_Disturbance_Rejection
        disp('Results saved as LQOC_Linear_Disturbance_Rejection');
    end
end