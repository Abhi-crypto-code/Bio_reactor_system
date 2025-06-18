% Bioreactor Simulation with Various Initial Conditions
clear all; clc; close all;

global BioReactor_mod;

% Set system parameters
BioReactor_mod.Yxs = 0.4;  % 0.4 g/L
BioReactor_mod.Beta = 0.2; % 0.2 h-1
BioReactor_mod.P_m = 50;   % 50 g/L
BioReactor_mod.K_i = 22;   % 22 g/L
BioReactor_mod.alpha = 2.2; % 2.2 g/L
BioReactor_mod.u_m = 0.48; % 0.48 h-1
BioReactor_mod.K_m = 1.2;  % 1.2 g/L
BioReactor_mod.S_f = 20;   % 20 g/L
BioReactor_mod.D = 0.2;    % 0.2 h-1

fprintf('\n\n Development of Linear Perturbation Models for Bioreactor system')

n_st = 3;  % No. of states                   
n_ip = 2;  % No. of inputs
n_op = 3;  % No. of outputs
n_ud = 0;  % No. of unmeasured disturbances

Us = [BioReactor_mod.D BioReactor_mod.S_f]'; % manipulated variables
Ws = 0; % no disturbance

fprintf('\n\n Bioreactor: Testing Different Initial Conditions');
fprintf('\n\n Steady State Inputs \n\n\t Dilution Rate: %6.4f h^-1 and Feed Substrate Concentration: %6.4f g/L', BioReactor_mod.D, BioReactor_mod.S_f);

% Define a set of initial guesses to test
initialGuesses = [
    0.4, 0.2, 80;   % Slightly below steady state
    0.6, 0.4, 120;  % Slightly above steady state
    0.4, 0.4, 120;  % X below, S and P above
    0.6, 0.2, 80;   % X above, S and P below
    0.3, 0.1, 60;   % All values more significantly below
    0.8, 0.6, 140;  % All values more significantly above
];

% Options for fsolve
options = optimoptions('fsolve', 'Display', 'iter', 'TolFun', 1e-10, 'TolX', 1e-10, 'MaxIter', 1000);

% Create a figure for plotting
figure('Position', [100, 100, 1200, 600]);

% Process each initial guess
for i = 1:size(initialGuesses, 1)
    fprintf('\n\n\nCASE %d: Testing Initial Guess [X=%g, S=%g, P=%g]', i, initialGuesses(i,1), initialGuesses(i,2), initialGuesses(i,3));
    
    % Set initial guess
    BioReactor_mod.X = initialGuesses(i,1);
    BioReactor_mod.S = initialGuesses(i,2);       
    BioReactor_mod.P = initialGuesses(i,3);
    Xs0 = [BioReactor_mod.X BioReactor_mod.S BioReactor_mod.P]'; 
    
    % Find steady state
    try
        fprintf('\n Finding steady state...');
        Xs = fsolve(@Bioreactor_SteadyState, Xs0, options);
        
        % Update model with steady state values
        BioReactor_mod.X = Xs(1);
        BioReactor_mod.S = Xs(2);
        BioReactor_mod.P = Xs(3);
        
        fprintf('\n Steady State Operating Point');
        fprintf('\n\t X: %6.4f  S: %6.4f  P: %6.4f', BioReactor_mod.X, BioReactor_mod.S, BioReactor_mod.P);
        
        % Compute Jacobian matrix and get eigenvalues
        Z_vec = [Xs' Us' Ws']';  
        Jacob_mat = Num_Jacobian(@Bioreactor_Model_JacobFn, Z_vec);
        A_mat = Jacob_mat(:, 1:n_st);
        B_mat = Jacob_mat(:, n_st+1:n_st+n_ip);
        H_mat = Jacob_mat(:, n_st+n_ip+1:n_st+n_ip+n_ud);
        
        fprintf('\n\n A matrix for Case %d =\n', i);
        disp(A_mat);
        
        eig_A = eig(A_mat);
        fprintf('\n Eigenvalues of A matrix for Case %d:\n', i);
        disp(eig_A);
        
        % Check stability
        if all(real(eig_A) < 0)
            fprintf('\n The system is STABLE for Case %d (all eigenvalues have negative real parts)\n', i);
        else
            fprintf('\n The system is UNSTABLE for Case %d (at least one eigenvalue has non-negative real part)\n', i);
        end
        
        % Simulate the system for this case
        tspan = [0 50]; % Simulate for 50 time units
        [t, y] = ode45(@(t,y) Bioreactor_Dynamics(t,y,Us), tspan, Xs0);
        
        % Plot the results
        subplot(2,3,i);
        plot(t, y(:,1), 'r-', t, y(:,2), 'g-', t, y(:,3)/100, 'b-');
        title(sprintf('Case %d: [X₀=%g, S₀=%g, P₀=%g]', i, initialGuesses(i,1), initialGuesses(i,2), initialGuesses(i,3)));
        xlabel('Time (h)');
        ylabel('Concentration');
        legend('X (Biomass)', 'S (Substrate)', 'P/100 (Product/100)');
        grid on;
        
    catch ME
        fprintf('\n Error in Case %d: %s\n', i, ME.message);
    end
end

% Add a common title
sgtitle('Bioreactor Simulation with Different Initial Conditions');

% Adjust the layout
set(gcf, 'Color', 'w');

% Print summary
fprintf('\n\n==== SUMMARY OF ALL CASES ====\n');
for i = 1:size(initialGuesses, 1)
    fprintf('Case %d: Initial Guess [X=%g, S=%g, P=%g]\n', i, initialGuesses(i,1), initialGuesses(i,2), initialGuesses(i,3));
end

% Helper functions for the simulation

function dXdt = Bioreactor_Dynamics(t, X, Us)
    global BioReactor_mod;
    
    % Extract states
    x = X(1);    % Biomass
    s = X(2);    % Substrate
    p = X(3);    % Product
    
    % Extract inputs
    D = Us(1);   % Dilution rate
    S_f = Us(2); % Feed substrate concentration
    
    % Calculate growth rate (Monod kinetics with inhibition)
    mu = BioReactor_mod.u_m * s / (BioReactor_mod.K_m + s) * (1 - p/BioReactor_mod.P_m);
    
    % State equations
    dx_dt = mu * x - D * x;
    ds_dt = D * (S_f - s) - (1/BioReactor_mod.Yxs) * mu * x;
    dp_dt = BioReactor_mod.alpha * mu * x + BioReactor_mod.Beta * x - D * p;
    
    dXdt = [dx_dt; ds_dt; dp_dt];
end

function F = Bioreactor_SteadyState(X)
    global BioReactor_mod;
    
    % Extract states
    x = X(1);
    s = X(2);
    p = X(3);
    
    % Calculate growth rate
    mu = BioReactor_mod.u_m * s / (BioReactor_mod.K_m + s) * (1 - p/BioReactor_mod.P_m);
    
    % Set the equations equal to zero for steady state
    F = zeros(3, 1);
    F(1) = mu * x - BioReactor_mod.D * x;
    F(2) = BioReactor_mod.D * (BioReactor_mod.S_f - s) - (1/BioReactor_mod.Yxs) * mu * x;
    F(3) = BioReactor_mod.alpha * mu * x + BioReactor_mod.Beta * x - BioReactor_mod.D * p;
end

function Jacob = Num_Jacobian(fun, Z)
    % Calculate numerical Jacobian matrix
    n = length(Z);
    Jacob = zeros(n-2, n);
    delta = 1e-6;
    
    % Calculate baseline function value
    F0 = feval(fun, Z);
    
    for i = 1:n
        % Perturb the i-th variable
        Z_perturb = Z;
        Z_perturb(i) = Z_perturb(i) + delta;
        
        % Calculate perturbed function value
        F_perturb = feval(fun, Z_perturb);
        
        % Compute the i-th column of Jacobian
        Jacob(:, i) = (F_perturb - F0) / delta;
    end
end

function F = Bioreactor_Model_JacobFn(Z)
    global BioReactor_mod;
    
    % Extract states, inputs and disturbances from Z
    x = Z(1);
    s = Z(2);
    p = Z(3);
    D = Z(4);
    S_f = Z(5);
    
    % Calculate growth rate
    mu = BioReactor_mod.u_m * s / (BioReactor_mod.K_m + s) * (1 - p/BioReactor_mod.P_m);
    
    % Set the equations
    F = zeros(3, 1);
    F(1) = mu * x - D * x;
    F(2) = D * (S_f - s) - (1/BioReactor_mod.Yxs) * mu * x;
    F(3) = BioReactor_mod.alpha * mu * x + BioReactor_mod.Beta * x - D * p;
end