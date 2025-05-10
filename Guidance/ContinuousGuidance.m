
clear; close all; clc;
format long g
plotStyle;
K = constant();
%% Ex 3.1: Plot debris density, compute initial and final state

% Define the debris density function q(rho), which depends on the radial distance (rho)
q = @(rho) K.k1./(K.k2 + ((rho - K.rho0) ./ K.DU).^2);

% Define the altitude range for the plot (from hi-100 to hf+100)
h_vec = linspace(K.hi - 100, K.hf + 100, 1e3);  % Altitude range for plotting
rho_vec = K.Re + h_vec;  % Compute radial distances from the center of Earth

% Plot debris density as a function of altitude
plot(h_vec, q(rho_vec), "Color", [0, 0.4470, 0.7410], 'DisplayName', 'q($\rho$)')
xlabel('Altitude [km]')  % Label for the x-axis
ylabel('Debris density [$DU^{-3}$]')  % Label for the y-axis, with units in DU^(-3)
xline(K.hi, '--', 'LineWidth',1.5,'Color',[0.4660, 0.6740, 0.1880], 'DisplayName','departure altitude');
xline(K.hf, '--', 'LineWidth',1.5,'Color',[0.8500, 0.3250, 0.0980], 'DisplayName','arrival altitude');
legend('Location', 'north', 'Orientation','vertical')
% Compute initial state (at altitude K.hi + K.Re)
ri = K.hi + K.Re;  % Initial radial distance from Earth's center
vci = sqrt(K.mu / ri);  % Initial velocity (circular orbit) using mu (standard gravitational parameter)
xxi = [ri; 0; 0; 0; vci; 0; K.m0];  % Initial state vector: [radius; velocity; mass]
disp('xxi')
disp(xxi)

% Compute final state (at altitude K.hf + K.Re)
rf = K.hf + K.Re;  % Final radial distance from Earth's center
vcf = sqrt(K.mu / rf);  % Final velocity (circular orbit) at the final altitude
rvf = [rf; 0; 0; 0; vcf * cos(K.di); vcf * sin(K.di)];  % Final state vector with inclination angle di
disp('rvf')
disp(rvf)
%% Ex 3.2: Adimensionalize the problem

% Convert initial state (xxi) to dimensionless units using characteristic units (DU, VU, MU)
rri_ad = xxi(1:3) / K.DU;  % Convert position (in km) to dimensionless units by dividing by DU
vvi_ad = xxi(4:6) / K.VU;  % Convert velocity (in km/s) to dimensionless units by dividing by VU
mi_ad  = xxi(7) / K.MU;    % Convert mass (in kg) to dimensionless units by dividing by MU

% Convert final state (rvf) to dimensionless units
rrf_ad = rvf(1:3) / K.DU;  % Convert position (in km) to dimensionless units by dividing by DU
vvf_ad = rvf(4:6) / K.VU;  % Convert velocity (in km/s) to dimensionless units by dividing by VU

% Create the dimensionless initial and final state vectors
xxi_ad = [rri_ad; vvi_ad; mi_ad];  % Initial state vector in dimensionless units
rvf_ad = [rrf_ad; vvf_ad];         % Final state vector in dimensionless units

%% Ex 3.4: solve the problem
% This section was commented out due to the high computational cost of 
% computing the initial conditions for the Lagrange multipliers. 
% Instead, the solution used in the report is provided. 
% If the user wants to verify the correctness of the code, he can 
% simply uncomment this section.
% 
% 
% X_opt = nan(8, 1);          % Placeholder for the optimized solution
% exitflag = -2;                  % Initialize exit flag to check solver status
% Nmax = 10;                      % Maximum number of iterations
% iter = 0;                        % Initialize iteration counter
% val = 1e9;                       % Initial large value for the function norm (used for convergence)
% 
% % Define solver options for fsolve
% %options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter', 'FunctionTolerance', 1e-7, 'MaxFunctionEvaluations', 5e4, 'MaxIterations', 2e4, 'FiniteDifferenceType', 'central');
% options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter', 'FunctionTolerance', 1e-7, 'MaxFunctionEvaluations', 2e4, 'MaxIterations', 1e3, 'FiniteDifferenceType', 'central');
% options.StepTolerance = 1e-10;  % Set step tolerance for solver
% 
% % Optimization loop: solve the system iteratively
% while exitflag <= 0 && iter < Nmax
% % Initialize random guess for the Lagrange multipliers and final time
%     ll0 = -250 + rand(7,1) * 500;
%     ll0(end) = rand() * 250;  
%     tf = 20 * pi - pi + 2 * pi * rand();           
% 
%     X0 = [ll0; tf];           % Create the new guess for the optimization variables
% 
%     % Call fsolve to solve the system of equations
%     [X, F, exitflag] = fsolve(@zero_finding_fun, X0, options);
% 
%     % Update the best solution if the current one is better (lower norm)
%     if norm(F) < val
%         X_opt = X;             % Store the new optimized values
%         val = norm(F);         % Update the best function norm
%     end
%     iter = iter + 1;          % Increment the iteration counter
% end
% 
% 
% %%
% 
% % Reset variables for the refinement process
% exitflag = -2;
% Nmax = 10;
% iter = 0;
% val = 1e9;
% perc = 1e-5;  % Perturbation fraction for refining the guess
% X_opt_old = X_opt;
% % Define options for fsolve with tighter tolerances for further refinement
% options = optimoptions('fsolve', 'Display', 'iter', 'FunctionTolerance', 1e-7, 'MaxFunctionEvaluations', 2e4, 'MaxIterations', 1e3, 'FiniteDifferenceType', 'central');
% options.StepTolerance = 1e-10;  % Set step tolerance for solver
% 
% % Optimization loop for refinement
% while exitflag <= 0 && iter < Nmax
%     X0 = nan(size(X_opt_old));  % Initialize new guess
%     for k = 1 : length(X_opt)
%         X0(k) = X_opt_old(k) + perc * X_opt_old(k) * (-1 + 2 * rand());  % Apply random perturbation to each element
%     end
% 
%     % Call fsolve again to refine the solution
%     [X, F, exitflag] = fsolve(@zero_finding_fun, X0, options);
% 
%     % Update the best solution if the norm of F is improved
%     if norm(F) < val
%         X_opt = X;
%         val = norm(F);
%     end
%     iter = iter + 1;  % Increment the iteration counter
% end

%% Simulation and Error Calculation

X_opt=[  -214.9812206843300;
        -10.365872735288638;
          0.885567766633861;
        -10.392920734051977;
         -214.6104511518981;
         -112.9453556547766;
          2.596449162488288;
        64.480106244461300]; 

ll0 = X_opt(1:7);
disp('ll0')
disp(ll0)

% Compute the residuals (errors) from the zero-finding function for the optimized solution
F = zero_finding_fun(X_opt);

% Retrieve constants from the K structure
K = constant();

% Initial state vector (dimensionless units)
Y0 = [K.xxi_ad; X_opt(1:end-1)];
tf = X_opt(end);
tspan = [0, tf];  % Time span for the simulation
disp('tf [min]')
disp(tf*K.TU/60)

% Define options for ODE solver with tight tolerances
options = odeset('RelTol', 3e-14, 'AbsTol', 3e-14);

% Solve the system of equations using ode113 (adaptive step-size solver)
[tt, YY] = ode113(@(t, Y) aug_state_rhs(t, Y), tspan, Y0, options);
YY = YY';  % Transpose the output matrix

m = YY(7, :);     % Extract the mass from the solution
disp('mf [kg]')
disp(m(end)*K.MU)

% Extract the position vector (3D) from the solution
rr = YY(1:3, :);

% Compute errors in position, velocity, Lagrange multiplier, and Hamiltonian
err_pos = norm(F(1:3)) * K.DU;  % Position error in km
err_vel = norm(F(4:6)) * K.VU * 1000;  % Velocity error in m/s

% Display the errors
fprintf('position error: %fe-8 [km] \n', err_pos*1e8)
fprintf('velocity error: %fe-8 [m/s] \n', err_vel*1e8)

%% Plotting Results
max_x = max(abs(rr(1, :)));
max_y = max(abs(rr(2, :)));
max_z = max(abs(rr(3, :)));

% Plot the trajectory in 3D space
figure
hold on
plot3(rr(1, :), rr(2, :), rr(3, :))  % Trajectory in 3D space
plot_circles(K.ri_ad)  % Plot initial orbit circle
plot_circles(K.rf_ad, K.di)  % Plot final orbit circle
xlabel('X [DU]')
ylabel('Y [DU]')
zlabel('Z [DU]')
quiver3(0, 0, 0, max_x, 0, 0, 'k', 'LineWidth', 1.5);  % X-axis
quiver3(0, 0, 0, 0, max_y, 0, 'k', 'LineWidth', 1.5);  % Y-axis
quiver3(0, 0, 0, 0, 0, max_z, 'k', 'LineWidth', 1.5);  % Z-axis
text(max_x * 1.1, 0, 0, 'X', 'FontSize', 12, 'Color', 'k'); % Label for X-axis
text(0, max_y * 1.1, 0, 'Y', 'FontSize', 12, 'Color', 'k'); % Label for Y-axis
text(0, 0, max_z * 1.1, 'Z', 'FontSize', 12, 'Color', 'k'); % Label for Z-axis
view(50, 25)

%Polar plot
r = vecnorm(rr);
right_asc = atan2(rr(2, :), rr(1, :));
decl   = asin(rr(3, :)./r);
r_log = log(r);
figure
polarplot(right_asc, r_log, 'Color', '#0072BD', 'DisplayName', '$log(r)$')
hold on
polarplot(0, r_log(1), 'diamond', 'MarkerSize', 10, 'MarkerEdgeColor', '#0072BD', 'MarkerFaceColor', '#0072BD', 'DisplayName', 'Departure')
polarplot(0, r_log(end), '^', 'MarkerSize', 10, 'MarkerEdgeColor', '#0072BD', 'MarkerFaceColor', '#0072BD', 'DisplayName', 'Arrival')
legend('Location','northeast', 'Orientation','vertical', 'FontSize',15)
rlim = [min(r_log) max(r_log)];
figure
plot(tt, r)
% Compute the Hamiltonian for each time step
H = nan(size(tt));
for k = 1 : length(H)
    H(k) = hamiltonian(YY(:, k));  % Compute the Hamiltonian at each time step
end

% Plot the Hamiltonian over time
figure
hold on
plot(tt, H)
xlabel('time [TU]')
ylabel('H [-]')

indices = find((decl(1:end-1) < 0 & decl(2:end) > 0) | (decl(1:end-1) > 0 & decl(2:end) < 0));

% Compute the primer vector (N, T, W) components at each time step
alpha_NTW = primer_vector(YY);
dark_grey = [0.3 0.3 0.3];

% Plot the individual components of the primer vector over time
figure
subplot(3, 1, 1)
hold on
plot(tt, alpha_NTW(1, :))
for k = 1 : length(indices)
    ind = indices(k);
    xline(tt(ind), '--', 'LineWidth', 1.5, 'Color', dark_grey)
end
xlabel('t [TU]')
ylabel('$\hat{\alpha}^*_N$ [-]', 'Rotation', 0)

subplot(3, 1, 2)
plot(tt, alpha_NTW(2, :))
for k = 1 : length(indices)
    ind = indices(k);
    xline(tt(ind), '--', 'LineWidth', 1.5, 'Color', dark_grey)
end
xlabel('t [TU]')
ylabel('$\hat{\alpha}^*_T$ [-]', 'Rotation', 0)

subplot(3, 1, 3)
plot(tt, alpha_NTW(3, :))
for k = 1 : length(indices)
    ind = indices(k);
    xline(tt(ind), '--', 'LineWidth', 1.5, 'Color', dark_grey)
end
xlabel('t [TU]')
ylabel('$\hat{\alpha}^*_W$ [-]', 'Rotation', 0)
%% Ex 3.5: Numerical Continuation

% Initialize constants
K = constant();  % Load constants such as spacecraft parameters and units
Xopt_num = nan(length(X_opt), length(K.Tvec_ad));  % Matrix to store optimal solutions for each thrust value
Xopt_num(:, length(K.Tvec_ad)) = X_opt;  % Set the last column of Xopt_num to the initial solution X_opt
X0 = Xopt_num(:, 1);  % Set the initial guess for the first iteration as the first column of Xopt_num

% Define solver options for fsolve with specific tolerances and iteration limits
options = optimoptions('fsolve', 'Display', 'iter',  ...
                       'OptimalityTolerance', 1e-10, 'FunctionTolerance', 1e-8, ...
                       'StepTolerance', 1e-10, 'FiniteDifferenceType', 'central');

% Perform the numerical continuation process by solving the system for decreasing thrust values
for k = length(K.Tvec_ad):-1:2
    % Use the current solution as the initial guess for the next thrust level
    X0 = Xopt_num(:, k);  
    
    % Define the function to minimize for the current thrust value (zero_finding_num)
    fun = @(X) zero_finding_num(X, k-1);  
    
    % Solve the system using fsolve
    [X, F, exitflag] = fsolve(fun, X0, options);  
    
    % Store the solution for the current thrust value
    Xopt_num(:, k-1) = X;  % Store the optimized values for the previous thrust level
end

%% Plotting the Results

X_opt =  Xopt_num(:, 1);
ll0 = X_opt(1:7);
disp('ll0')
disp(ll0)
F = zero_finding_num(X_opt, 1);
% Set initial conditions for the simulation based on the optimized solution
Y0 = [K.xxi_ad; X_opt(1:end-1)];  % Initial state vector from optimized parameters
tf = X_opt(end);  % Time of flight from the optimized solution
tspan = [0, tf];  % Define the time span for the simulationclear
disp('tf [min]')
disp(tf*K.TU/60)

% Define solver options for ode113 with high precision tolerances
options = odeset('RelTol', 3e-14, 'AbsTol', 3e-14);

% Solve the low-thrust dynamics using ode113
[tt, YY] = ode113(@(t, Y) low_thrust_dynamics(t, Y, 1), tspan, Y0, options);

YY = YY';  % Transpose the solution for easier handling
rr = YY(1:3, :);  % Extract the position from the solution
m = YY(7, :);     % Extract the mass from the solution
lm = YY(14, :);   % Extract the mass Lagrange multiplier from the solution

err_pos = norm(F(1:3)) * K.DU;  % Position error in km
err_vel = norm(F(4:6)) * K.VU * 1000; 

% Display the errors
fprintf('position error: %fe-8 [km] \n', err_pos*1e8)
fprintf('velocity error: %fe-8 [m/s] \n', err_vel*1e8)

disp('mf [kg]')
disp(m(end)*K.MU)

%Polar plot
r = vecnorm(rr);
right_asc = atan2(rr(2, :), rr(1, :));
decl   = asin(rr(3, :)./r);
r_log = log(r);
figure
polarplot(right_asc, r_log, 'Color', '#0072BD', 'DisplayName', '$log(r)$')
hold on
polarplot(0, r_log(1), 'diamond', 'MarkerSize', 10, 'MarkerEdgeColor', '#0072BD', 'MarkerFaceColor', '#0072BD', 'DisplayName', 'Departure')
polarplot(0, r_log(end), '^', 'MarkerSize', 10, 'MarkerEdgeColor', '#0072BD', 'MarkerFaceColor', '#0072BD', 'DisplayName', 'Arrival')
legend('Location','northeast', 'Orientation','vertical', 'FontSize',15)
rlim = [min(r_log) max(r_log)];
figure
plot(tt, r)
% Compute the Hamiltonian for each time step
H = nan(size(tt));
for k = 1 : length(H)
    H(k) = hamiltonian_num(YY(:, k), 1);  % Compute the Hamiltonian at each time step
end

% Plot the Hamiltonian over time
figure
hold on
plot(tt, H)
xlabel('time [TU]')
ylabel('H [-]')

indices = find((decl(1:end-1) < 0 & decl(2:end) > 0) | (decl(1:end-1) > 0 & decl(2:end) < 0));
%%
% Compute the primer vector (N, T, W) components at each time step
alpha_NTW = primer_vector(YY);
dark_grey = [0.3 0.3 0.3];

% Plot the individual components of the primer vector over time
figure
subplot(3, 1, 1)
hold on
plot(tt, alpha_NTW(1, :))
for k = 1 : length(indices)
    ind = indices(k);
    xline(tt(ind), '--', 'LineWidth', 1.5, 'Color', dark_grey)
end
xlabel('t [TU]')
ylabel('$\hat{\alpha}_N^*$ [-]', 'Rotation', 0)

subplot(3, 1, 2)
plot(tt, alpha_NTW(2, :))
for k = 1 : length(indices)
    ind = indices(k);
    xline(tt(ind), '--', 'LineWidth', 1.5, 'Color', dark_grey)
end
xlabel('t [TU]')
ylabel('$\hat{\alpha}_T^*$ [-]', 'Rotation', 0)

subplot(3, 1, 3)
plot(tt, alpha_NTW(3, :))
for k = 1 : length(indices)
    ind = indices(k);
    xline(tt(ind), '--', 'LineWidth', 1.5, 'Color', dark_grey)
end
xlabel('t [TU]')
ylabel('$\hat{\alpha}_W^*$ [-]', 'Rotation', 0)

%%
function plotStyle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set figure properties for better looking plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interpreter:
set(0, 'defaultTextInterpreter', 'Latex')
set(0, 'defaultAxesTickLabelInterpreter', 'Latex')
set(0, 'defaultLegendInterpreter', 'Latex')
set(0,'defaultAxesXGrid','on')
set(0,'defaultAxesYGrid','on')
% lines:
set(0,'defaultLineLineWidth', 1.5);
set(0,'defaultLineMarkerSize',6) ;
set(0,'defaultLineMarkerEdgeColor','k')
set(0,'defaultLineMarkerFaceColor','auto')
% legend:
set(0, 'defaultLegendLocation','northoutside');
set(0, 'defaultLegendOrientation','horizontal');
set(0, 'defaultLegendFontSize',12);
% axes:
set(0,'defaultAxesFontSize',16);
end

function K = constant()

% Initialize a structure to store constants and parameters
K = struct();

% Physical and Orbital Parameters
K.hi = 800;                         % [km] Initial orbit altitude
K.hf = 1000;                        % [km] Final orbit altitude
K.di = deg2rad(0.75);               % [rad] Inclination change (converted from degrees)
K.Re = 6378.1366;                   % [km] Earth's radius
K.mu = 398600.435;                  % [km^3/s^2] Gravitational parameter of Earth
K.rho0 = 750 + K.Re;                % [km] Reference radius for debris flux

% Reference Units
K.DU = 7178.1366;                   % [km] Distance Unit (reference length)
K.k1 = 1e-5;                        % [DU^-1] Debris density constant 1
K.k2 = 1e-4;                        % [DU^2] Debris density constant 2
K.m0 = 1000;                        % [kg] Initial spacecraft mass
K.Tmax = 3;                         % [N] Maximum thrust
K.Is = 3120;                        % [s] Specific impulse of the thruster

% Derived Units and Adimensionalized Parameters
K.MU = K.m0;                        % [kg] Mass Unit (set to the initial mass)
K.TU = sqrt(K.DU^3 / K.mu);         % [s] Time Unit (derived from scaling gravitational parameter)
K.VU = K.DU / K.TU;                 % [km/s] Velocity Unit
K.g0 = 9.81;                        % [m/s^2] Standard gravitational acceleration
K.Tmax_ad = K.Tmax / (K.MU * K.VU * 1e3 / K.TU); % [adim] Adimensionalized thrust
K.Is_ad = K.Is / K.TU;              % [adim] Adimensionalized specific impulse
K.g0_ad = K.g0 / (K.VU * 1e3 / K.TU); % [adim] Adimensionalized gravitational acceleration
K.rho0_ad = K.rho0 / K.DU;          % [adim] Adimensionalized reference radius for debris flux

% Initial and Final Orbital States
K.ri = K.hi + K.Re;                 % [km] Initial orbit radius
K.vci = sqrt(K.mu / K.ri);          % [km/s] Circular velocity at the initial orbit
K.xxi = [K.ri; 0; 0; 0; K.vci; 0; K.m0]; % [km, km/s, kg] Initial state vector [position; velocity; mass]

K.rf = K.hf + K.Re;                 % [km] Final orbit radius
K.vcf = sqrt(K.mu / K.rf);          % [km/s] Circular velocity at the final orbit
K.rvf = [K.rf; 0; 0; 0; K.vcf*cos(K.di); K.vcf*sin(K.di)]; % Final state vector in rotated frame

% Adimensionalized State Vectors
K.rri_ad = K.xxi(1:3) / K.DU;       % [adim] Initial position vector
K.vvi_ad = K.xxi(4:6) / K.VU;       % [adim] Initial velocity vector
K.mi_ad = K.xxi(7) / K.MU;          % [adim] Initial mass
K.rrf_ad = K.rvf(1:3) / K.DU;       % [adim] Final position vector
K.vvf_ad = K.rvf(4:6) / K.VU;       % [adim] Final velocity vector
K.xxi_ad = [K.rri_ad; K.vvi_ad; K.mi_ad]; % [adim] Adimensionalized initial state vector
K.rvf_ad = [K.rrf_ad; K.vvf_ad];    % [adim] Adimensionalized final state vector

K.ri_ad = K.ri / K.DU;              % [adim] Adimensionalized initial orbit radius
K.rf_ad = K.rf / K.DU;              % [adim] Adimensionalized final orbit radius

% Optimization Bounds
K.lb = [-250 * ones(6,1); 0];       % Lower bounds for optimization variables
K.ub = 250 * ones(7,1);             % Upper bounds for optimization variables

% Thrust Levels for Analysis
K.Tvec = linspace(2.860, 3.000, 5); % [N] Range of thrust levels for analysis
K.Tvec_ad = K.Tvec / (K.MU * K.VU * 1e3 / K.TU); % [adim] Adimensionalized thrust range

end

function Ydot = aug_state_rhs(~, Y)
% Computes the time derivative of the augmented state vector.
% Includes state dynamics and costate dynamics for the low-thrust orbit-raising problem.

% Initialize the derivative vector
Ydot = nan(size(Y));

% Load constants
K = constant;

% Decompose the state vector
rr = Y(1:3);      % Position vector
vv = Y(4:6);      % Velocity vector
m = Y(7);         % Spacecraft mass
llr = Y(8:10);    % Costate associated with position
llv = Y(11:13);   % Costate associated with velocity

% Compute derived quantities
r = norm(rr);     % Magnitude of the position vector
lv = norm(llv);   % Magnitude of the velocity costate

% State dynamics
Ydot(1:3) = vv;                              % Position derivative
Ydot(4:6) = -rr / r^3 - K.Tmax_ad / m * llv / lv; % Velocity derivative (acceleration due to gravity and thrust)
Ydot(7) = -K.Tmax_ad / (K.Is_ad * K.g0_ad);  % Mass derivative (due to propellant usage)

% Costate dynamics
Ydot(8:10) = ...
    K.k1 * 2 * (r - K.rho0_ad) / (K.k2 + (r - K.rho0_ad)^2)^2 * rr / r ... % Debris density gradient term
    - 3 / r^5 * dot(rr, llv) * rr ... % Gravitational gradient term
    + llv / r^3;                     % Coupling with velocity costate
Ydot(11:13) = -llr;                  % Velocity costate derivative
Ydot(14) = -lv * K.Tmax_ad / m^2;    % Mass costate derivative

end

function H = hamiltonian(Y)
    % Constants and Parameters
    K = constant();  % Retrieve constants from the constant function
    rr = Y(1:3);      % Position vector 
    vv = Y(4:6);      % Velocity vector 
    m = Y(7);         % Mass 
    llr = Y(8:10);    % Lagrange multipliers for position 
    llv = Y(11:13);   % Lagrange multipliers for velocity 
    lm = Y(14);       % Lagrange multiplier for mass 

    % Magnitude of position vector
    r = norm(rr);    

    % debris density function 
    l = K.k1 / (K.k2 + (r - K.rho0_ad)^2);

    % Magnitude of the Lagrange multipliers for velocity
    lv = norm(llv);

    % Hamiltonian formula, considering various energy terms and constraints
    H = l + dot(llr, vv) - dot(llv, rr) / r^3 - K.Tmax_ad / m * lv - ...
        lm * K.Tmax_ad / (K.Is_ad * K.g0_ad);
end

function H = hamiltonian_num(Y, k)

% Constants
K = constant();
Tmax_ad = K.Tvec_ad(k);  % Adapted maximum thrust for current iteration

% State variables
rr = Y(1:3);   % Position vector
vv = Y(4:6);   % Velocity vector
m = Y(7);      % Mass
llr = Y(8:10); % Lagrange multipliers for position 
llv = Y(11:13);% Lagrange multipliers for velocity 
lm = Y(14);    % Lagrange multiplier for mass 

% Calculations
r = norm(rr);  % Magnitude of position vector
l   = K.k1/(K.k2+(r-K.rho0_ad)^2);  % Some potential function
lv = norm(llv);  % Magnitude of the velocity Lagrange multipliers

% Hamiltonian
H = l + dot(llr, vv) - dot(llv, rr)/r^3 - Tmax_ad/m * lv - ...
    lm * Tmax_ad / (K.Is_ad * K.g0_ad);  

end

function F = zero_finding_fun(X)
    % Constants
    K = constant();  
    
    % Initial conditions setup
    Y0 = [K.xxi_ad; X(1:end-1)]; % Initial state (admissible values and variables from X)
    tf = X(end);  % Final time (last element of X)
    tspan = [0, tf];  % Time span for the ODE solver

    % ODE solver options (high precision)
    options = odeset('RelTol', 2.5e-14, 'AbsTol', 2.5e-14);

    % Solve the augmented system of ODEs
    [~, YY] = ode113(@(t, Y) aug_state_rhs(t, Y), tspan, Y0, options);

    % Final state vector (state at final time)
    Yf = YY(end, :).';

    % Extract final state components
    rv_tf = Yf(1:6);  % Final position and velocity
    lm    = Yf(end);   % Final Lagrange multiplier for mass
    H     = hamiltonian(Yf);  % Hamiltonian evaluated at final state

    % Function value used for zero-finding
    F = [rv_tf - K.rvf_ad; lm; H];

    % Boundary check: enforce bounds on the variables
    if any(X(1:7) - K.lb < 0) || any(X(1:7) - K.ub > 0)
        F = F + 1e6 * ones(size(F));  % Add large penalty if out of bounds
    end
end

function plot_circles(r, di)
    % r: radius of the orbit
    % di: inclination angle (angle between orbital plane and reference plane)

    if nargin < 2
        % If no inclination angle is provided, default inclination to 0
        di = 0;
    end

    % Create an array of 100 equally spaced angles between 0 and 2*pi (full circle)
    angles = linspace(0, 2*pi, 1e2);

    % Initialize a 3xN array to store the rotated circle coordinates in 3D
    rr = nan(3, length(angles));

    % Loop through each angle to calculate the positions
    for k = 1 : length(angles)
        % Position on the circle in the XY-plane
        theta = angles(k);
        pos = r * [cos(theta); sin(theta); 0];  % Circle in the XY plane

        % Rotation matrix R to rotate the circle by the inclination angle (di) around the X-axis
        R = [1, 0, 0; 
             0, cos(di), -sin(di); 
             0, sin(di), cos(di)];

        % Apply the rotation to the position vector
        rr(:, k) = R * pos;
    end

    % Plot the circle in 3D space
    plot3(rr(1, :), rr(2, :), rr(3, :));

end

function alpha_NTW = primer_vector(YY)
% primer_vector computes the primer vector for a set of state vectors (YY)
% The primer vector is derived from the Lagrange multipliers and represents 
% the direction of the thrust in the NTW frame.

% Initialize the output matrix to store the primer vectors
alpha_NTW = nan(3, size(YY, 2));

% Loop through each state vector in the input matrix YY
for k = 1 : size(alpha_NTW, 2)

    % Extract the position (rr) and velocity (vv) from the state vector rv
    rv  = YY(1:6,k);  % First 6 elements are position and velocity
    rr = rv(1:3);     % Position vector (3 elements)
    vv = rv(4:6);     % Velocity vector (3 elements)
    
    % Compute the angular momentum vector (hh) as the cross product of position and velocity
    hh = cross(rr, vv);

    % Define the unit vectors of the orbital frame:
    w = -rr/norm(rr); % Radial unit vector (opposite direction of rr)
    t = vv/norm(vv);  % Tangential unit vector (direction of velocity)
    n = hh/norm(hh);  % Normal unit vector (perpendicular to the orbital plane)

    % Construct the rotation matrix R from the orbital frame to the inertial frame
    R = [n,t,w]\eye(3); % Inverse of the matrix formed by [n t w]

    % Extract the Lagrange multiplier for velocity (llv) from the state vector
    llv = YY(11:13,k);  % Lagrange multipliers for velocity components (3 elements)
    
    % Normalize the Lagrange multiplier vector and negate it to get the direction
    alpha = -llv/norm(llv);

    % Transform the Lagrange multiplier into the inertial frame using the rotation matrix
    alpha = R * alpha;

    % Normalize the primer vector and store it in the output matrix
    alpha_NTW(:, k) = alpha / norm(alpha);  % Normalize and store the result
end

end

function Ydot = low_thrust_dynamics(~, Y, k)
% low_thrust_dynamics computes the state derivative for a system with low-thrust dynamics.
% This function is used in numerical continuation for the integration of the system's motion.
% The system includes the position, velocity, mass, and Lagrange multipliers.

% Initialize the output derivative vector
Ydot = nan(size(Y));

% Load constants from the constant function
K = constant;

% Get the thrust magnitude for the current step in the numerical continuation process
Tmax_ad = K.Tvec_ad(k);

% Extract state variables from the input vector Y
rr = Y(1:3);   % Position vector (3 elements)
vv = Y(4:6);   % Velocity vector (3 elements)
m  = Y(7);     % Mass of the spacecraft
llr = Y(8:10); % Lagrange multipliers for position (3 elements)
llv = Y(11:13);% Lagrange multipliers for velocity (3 elements)

% Compute the distance from the origin
r = norm(rr);

% Compute the norm (magnitude) of the Lagrange multiplier for velocity
lv = norm(llv);

% Compute the state derivatives:

% Position dynamics (velocity is the derivative of position)
Ydot(1:3) = vv;

% Velocity dynamics, including gravitational force and low-thrust acceleration
Ydot(4:6) = -rr/r^3 - Tmax_ad/m * llv/lv;  % Gravitational acceleration and low-thrust contribution

% Mass dynamics, driven by the rate of fuel usage (dependent on Tmax_ad)
Ydot(7) = -Tmax_ad / (K.Is_ad * K.g0_ad);  % Fuel depletion rate

% Lagrange multiplier dynamics for position, derived from the specific system's equations
Ydot(8:10) = K.k1 * 2 * (r - K.rho0_ad) / (K.k2 + (r - K.rho0_ad)^2)^2 * rr / r ...
    - 3 / r^5 * dot(rr, llv) * rr + llv / r^3;

% Lagrange multiplier dynamics for velocity (negative of llr)
Ydot(11:13) = -llr;

% Lagrange multiplier dynamics for mass, influenced by the rate of thrust application
Ydot(14) = -lv * Tmax_ad / m^2;  % Mass Lagrange multiplier dynamics
end

function F = zero_finding_num(X, k)
% zero_finding_num computes the residual of the system state at the final time step
% for numerical continuation. It integrates the low-thrust dynamics from the initial state
% and compares the final state with a desired reference (rvf_ad), along with the Lagrange multiplier 
% and Hamiltonian. The result is used for root-finding (e.g., in numerical continuation).
%
% Input:
%   - X: A vector containing the current guess for the system state and final time.
%   - k: The index for the current numerical continuation step (used to get specific thrust data).
%
% Output:
%   - F: The residual vector, consisting of the differences between the final state and desired values
%        (desired velocity, Lagrange multiplier, and Hamiltonian).

% Load constants from the constant function
K = constant();

% Construct the initial state vector Y0 by combining the desired initial conditions (xxi_ad) and the 
% guessed control variables (X). The last element of X represents the final time.
Y0 = [K.xxi_ad; X(1:end-1)];

% Extract the final time tf from the input vector X
tf = X(end);

% Define the time span for the ODE solver, from 0 to tf
tspan = [0, tf];

% Set the options for the ODE solver with very high precision (relative and absolute tolerances)
options = odeset('RelTol', 3e-14, 'AbsTol', 3e-14);

% Integrate the system dynamics using ode113 (numerical ODE solver)
[~, YY] = ode113(@(t, Y) low_thrust_dynamics(t, Y, k), tspan, Y0, options);

% Extract the final state (last column of the solution)
Yf = YY(end, :).';

% Extract the final position and velocity vector (rv_tf)
rv_tf = Yf(1:6);

% Extract the final value of the Lagrange multiplier (lm) from the state
lm = Yf(end);

% Compute the Hamiltonian for the final state
H = hamiltonian_num(Yf, k);

% Compute the residual vector F:
% - The first part compares the final state velocity with the desired reference (rvf_ad).
% - The second part is the residual of the Lagrange multiplier (lm).
% - The third part is the residual of the Hamiltonian (H).
F = [rv_tf - K.rvf_ad; lm; H];

end
