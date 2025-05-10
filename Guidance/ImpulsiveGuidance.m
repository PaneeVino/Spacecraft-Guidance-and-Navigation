clear; close all; clc;
format long g

plotStyle
% addpath("kernels/")
% addpath("../mice/src/mice/")
% addpath("../mice/lib/")
% % Clear SPICE kernel and load necessary SPICE kernels
cspice_kclear();  % Clear any existing kernels
cspice_furnsh('ex02.tm');  % Load SPICE kernel file 'ex02.tm'

%In the multiple shooting the commented options are those used to obtain the graphs in point 3. 
% They were commented out due to the computational cost of multiple shooting.
% The uncommented options (default MaxIterations and MaxFunctionEvaluations) 
% are those used to obtain the results presented in APPENDIX F of the report.
% To visualize the results presented in the report, simply uncomment the 
% currently commented options.
 %% Ex2.1: Plot Initial Guess 

% Retrieve problem constants
K = constant();

% Compute the initial state and time span from parameters
% Y0 contains the initial state in the synodic frame and time span [x0; y0; vx0; vy0; ti; tf]
Y0 = par2car(K.alpha, K.beta, K.ti, K.delta, K.ri, K.mu);
%Y0 = par2car(1.5*pi, 1.41, 0, 7, K.ri, K.mu);
% Extract initial state (xx0) and time span (tspan)
xx0 = Y0(1:4);       % Initial state: [x0; y0; vx0; vy0]
tspan = Y0(5:end);   % Time span: [ti; tf]

% Set options for the ODE solver
% High precision is ensured using RelTol and AbsTol
options = odeset('RelTol', 2.5e-14, 'AbsTol', 2.5e-14);

% Integrate the equations of motion in the synodic frame
[tt_0, xx_0] = ode113(@(t,xx) PBRFBP_rhs(t,xx, K.mu, K.rho, K.oms, K.ms), ...
    tspan, xx0, options);
xx_M = zeros(length(tt_0), 4);
xx_M(:,1) = 1-K.mu;

disp('x_0_guess:')
disp(xx0)

% Compute trajectory in the inertial (ECI) frame
% Use the rotation transformation for the synodic-to-inertial conversion
XX_0 = rot2ECI(xx_0, tt_0, K.mu);
XX_M = rot2ECI(xx_M, tt_0, K.mu);
XX_M_0 = XX_M(1,:);

% Plot the trajectory in the synodic frame
figure
hold on
axis equal
grid minor
box on
plot(xx_0(:, 1), xx_0(:, 2),'Color', "#0072BD", 'DisplayName', 'S/C Trajectory')  % Plot x-y trajectory in the synodic frame
plot_circle(K.ri, -K.mu, 0, 'Parking Orbit', "#77AC30")
xlabel('x [-]')  % Label x-axis
ylabel('y [-]')  % Label y-axis
legend('S/C Trajectory', 'Parking Orbit')
axes('position',[.4 .35 .19 .19])
box on % put box around new pair of axes
hold on
plot(xx_0(1:150, 1), xx_0(1:150, 2),'Color', "#0072BD", 'DisplayName', 'none')
plot_circle(K.ri, -K.mu, 0, 'Parking Orbit', "#77AC30")
axis equal

% Plot the trajectory in the ECI frame
figure
hold on
axis equal
grid minor
box on
plot(XX_0(:, 1), XX_0(:, 2),'Color', "#0072BD", 'DisplayName', 'S/C Trajectory')
plot(XX_M(:, 1), XX_M(:, 2),'--', 'Color', "#D95319", 'DisplayName', 'Moon Orbit')
plot(XX_0(1, 1), XX_0(1, 2), '^', 'MarkerSize', 10, 'MarkerEdgeColor', "#0072BD", 'MarkerFaceColor', "#0072BD",'DisplayName', 'Starting Point')
plot(XX_0(end, 1), XX_0(end, 2), 'diamond', 'MarkerSize', 10,  'MarkerEdgeColor', "#0072BD", 'MarkerFaceColor', "#0072BD",'DisplayName', 'Arrirval Point')
plot(XX_M_0(1), XX_M_0(2), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'DisplayName', 'Initial Moon position')
xlabel('X [-]')  % Label X-axis
ylabel('Y [-]')  % Label Y-axis
legend('Location','northeast', 'Orientation','vertical')
%% Ex2.2.a: Simple shooting without providing any derivative

clc
% Define bounds for the decision variables
% [x0, y0, vx0, vy0, t0, tf]
lb = [-(K.mu+K.ri); -K.ri; -sqrt(2)*K.vi; -sqrt(2)*K.vi; 0; 0];  % Lower bounds
ub = [-K.mu+K.ri; K.ri; sqrt(2)*K.vi; sqrt(2)*K.vi; -2*pi/K.oms; K.T];  % Upper bounds

% Define linear inequality constraints: A*Y = b
A = [0, 0, 0, 0, 1, -1];  % Ensures t0 - tf <= 0
b = 0;

% Set optimization options for `fmincon`
options = optimoptions('fmincon', ...
    'Display', 'iter-detailed', ...              
    'Algorithm', 'active-set', ...              
    'SpecifyObjectiveGradient', false, ...     
    'SpecifyConstraintGradient', false, ...     
    'ConstraintTolerance', 1e-10, ...
    'MaxFunctionEvaluations', 1e3);  

% Run optimization using `fmincon`
[Y_simple, dv_simple, exitflag_simple, output_simple, ~, ~, hessian_simple] = ...
    fmincon(@obj_function, Y0, A, b, [], [], lb, ub, @myconstrains, options);

Iter    = output_simple.iterations;
FunEval = output_simple.funcCount;
MaxCon  = output_simple.constrviolation;
dt = (Y_simple(end)-Y_simple(end-1))*K.TU;
dv = dv_simple*K.VU;

components_out = {'dv [km/s]', 'dt [days]', 'iterations', 'funeval', 'MaxCon'}';
vec_out = [dv; dt; Iter; FunEval; MaxCon];
tab_out = table(components_out, vec_out, 'VariableNames', {'Parameter', 'Value'});
disp(tab_out)

% Evaluate constraints for the optimal solution
[c, ceq] = myconstrains(Y_simple);
% Create the table
components = {'x', 'y', 'v_x', 'v_y', 't_i', 't_f'}';
Y_simple_table = table(components, Y_simple, 'VariableNames', {'Component', 'Value'});

% Display the table
disp(Y_simple_table);

% Compute and display additional outputs
disp(['Delta-v: ', num2str(dv_simple)]);
disp(['Exit flag: ', num2str(exitflag_simple)]);
disp(['Minimum eigenvalue of Hessian: ', num2str(min(eig(hessian_simple)))]);

% Display the results
disp('Inequality constraints (c):');
disp(c);

disp('Equality constraints (ceq):');
disp(ceq);


% Plot simple shooting optimal solution
% Extract the optimized initial state and time span
xx0_simple = Y_simple(1:4);  % Initial state [x0, y0, vx0, vy0]
tspan = Y_simple(5:end);     % Optimized time span [t0, tf]

% Set ODE solver options
options = odeset('RelTol', 2.5e-14, 'AbsTol', 2.5e-14);

% Integrate the equations of motion in the synodic frame
[tt_simple, xx_simple] = ode113(@(t,xx) PBRFBP_rhs(t,xx, K.mu, K.rho, K.oms, K.ms), ...
    tspan, xx0_simple, options);

% Compute trajectory in the inertial (ECI) frame
XX_simple = rot2ECI(xx_simple, tt_simple, K.mu);
xx_M = zeros(length(tt_simple), 4);
xx_M(:,1) = 1-K.mu;
XX_M = rot2ECI(xx_M, tt_simple, K.mu);
XX_M_0 = XX_M(1,:);

% Plot the optimized trajectory in the synodic frame
figure
hold on
axis equal
grid minor
box on
plot(xx_simple(:, 1), xx_simple(:, 2),'Color', "#0072BD", 'DisplayName', 'S/C Trajectory')  % Plot x-y trajectory in synodic frame
plot_circle(K.ri, -K.mu, 0, 'Parking Orbit', [0.4660, 0.6740, 0.1880])
plot_circle(K.rf, 1-K.mu, 0, 'Arrival Orbit', [0.8500, 0.3250, 0.0980])
legend('S/C Trajectory', 'Parking Orbit', 'Arrival Orbit')
axes('position',[.3 .3 .19 .19])
box on % put box around new pair of axes
hold on
plot(xx_simple(1:150, 1), xx_simple(1:150, 2),'Color', "#0072BD", 'DisplayName', 'none')
plot_circle(K.ri, -K.mu, 0, 'Parking Orbit', [0.4660, 0.6740, 0.1880])
axis equal
axes('position',[.6 .3 .19 .19])
box on % put box around new pair of axes
hold on
plot(xx_simple(end-80:end, 1), xx_simple(end-80:end, 2),'Color', "#0072BD", 'DisplayName', 'none')
plot_circle(K.rf, 1-K.mu, 0, 'Arrival Orbit', [0.8500, 0.3250, 0.0980])
axis equal

% Plot the optimized trajectory in the ECI frame
figure
hold on
axis equal
grid minor
box on
plot(XX_simple(:, 1), XX_simple(:, 2),'Color', "#0072BD", 'DisplayName', 'S/C Trajectory') 
plot(XX_M(:, 1), XX_M(:, 2),'--', 'Color', "#D95319", 'DisplayName', 'Moon Orbit')
plot(XX_simple(1, 1), XX_simple(1, 2), '^', 'MarkerSize', 10, 'MarkerEdgeColor', "#0072BD", 'MarkerFaceColor', "#0072BD",'DisplayName', 'Starting Point')
plot(XX_simple(end, 1), XX_simple(end, 2), 'diamond', 'MarkerSize', 10,  'MarkerEdgeColor', "#0072BD", 'MarkerFaceColor', "#0072BD",'DisplayName', 'Arrirval Point')
plot(XX_M_0(1), XX_M_0(2), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'DisplayName', 'Initial Moon position')
xlabel('X [-]')  % Label X-axis
ylabel('Y [-]')  % Label Y-axis
legend('Location','northeast', 'Orientation','vertical')
%% Ex2.2.b: Simple shooting providing derivatives
clc

% Define the bounds for the decision variables in the optimization
% Initial state variables (position and velocity) and final time.
lb = [-(K.mu+K.ri); -K.ri; -sqrt(2)*K.vi; -sqrt(2)*K.vi; 0; 0]; % Lower bounds
ub = [-K.mu+K.ri; K.ri; sqrt(2)*K.vi; sqrt(2)*K.vi; -2*pi/K.oms; K.T]; % Upper bounds

% Define the linear constraint: A*Y = b
A  = [0, 0, 0, 0, 1, -1]; % Constraint ensures that t0 - tf = 0
b  = 0;  % b = 0 for the linear constraint

% Set optimization options for `fmincon`
options = optimoptions('fmincon', ...
    'Display', 'iter-detailed', ...              
    'Algorithm', 'active-set', ...              
    'SpecifyObjectiveGradient', true, ...     
    'SpecifyConstraintGradient', true, ...     
    'ConstraintTolerance', 1e-10);  

% Perform optimization using fmincon
[Y_simple, dv_simple, exitflag_simple, output_simple, ~, ~, hessian_simple] = ...
    fmincon(@obj_function, Y0, A, b, [], [], lb, ub, @myconstrains, options);

Iter    = output_simple.iterations;
FunEval = output_simple.funcCount;
MaxCon  = output_simple.constrviolation;
dt = (Y_simple(end)-Y_simple(end-1))*K.TU;
dv = dv_simple*K.VU;

Deltav_simple = dv*1e3;
Deltat_simple = dt;

components_out = {'dv [km/s]', 'dt [days]', 'iterations', 'funeval', 'MaxCon'}';
vec_out = [dv; dt; Iter; FunEval; MaxCon];
tab_out = table(components_out, vec_out, 'VariableNames', {'Parameter', 'Value'});
disp(tab_out)

% Define the component names
components = {'x', 'y', 'v_x', 'v_y', 't_i', 't_f'}';

% Evaluate constraints for the optimal solution
[c, ceq] = myconstrains(Y_simple);
% Create the table
Y_simple_table = table(components, Y_simple, 'VariableNames', {'Component', 'Value'});

% Display the table
disp(Y_simple_table);

% Compute and display additional outputs
disp(['Delta-v: ', num2str(dv_simple)]);
disp(['Exit flag: ', num2str(exitflag_simple)]);
disp(['Minimum eigenvalue of Hessian: ', num2str(min(eig(hessian_simple)))]);

% Display the results
disp('Inequality constraints (c):');
disp(c);

disp('Equality constraints (ceq):');
disp(ceq);

% Plot the simple shooting optimal solution
% Extract the optimized initial state and time span from the solution vector Y_simple
xx0_simple = Y_simple(1:4);  % Initial position and velocity
tspan = Y_simple(5:end);     % Optimized time span [t0, tf]

% Set ODE solver options with high accuracy
options = odeset('RelTol', 2.5e-14, 'AbsTol', 2.5e-14);

% Solve the ODE system using ode113 to get the trajectory of the spacecraft
[tt_simple, xx_simple] = ode113(@(t, xx) PBRFBP_rhs(t, xx, K.mu, K.rho, K.oms, K.ms), ...
    tspan, xx0_simple, options);

% Convert trajectory from Synodic to ECI frame
XX_simple = rot2ECI(xx_simple, tt_simple, K.mu);
xx_M = zeros(length(tt_simple), 4);
xx_M(:,1) = 1-K.mu;
XX_M = rot2ECI(xx_M, tt_simple, K.mu);
XX_M_0 = XX_M(1,:);

% Plot the optimized trajectory in the synodic frame
figure
hold on
axis equal
grid minor
box on
plot(xx_simple(:, 1), xx_simple(:, 2),'Color', "#0072BD", 'DisplayName', 'S/C Trajectory')  % Plot x-y trajectory in synodic frame
plot_circle(K.ri, -K.mu, 0, 'Parking Orbit', [0.4660, 0.6740, 0.1880])
plot_circle(K.rf, 1-K.mu, 0, 'Arrival Orbit', [0.8500, 0.3250, 0.0980])
xlabel('x [-]')  % Label x-axis
ylabel('y [-]')  % Label y-axis
legend('S/C Trajectory', 'Parking Orbit', 'Arrival Orbit')
axes('position',[.3 .3 .19 .19])
box on % put box around new pair of axes
hold on
plot(xx_simple(1:150, 1), xx_simple(1:150, 2),'Color', "#0072BD", 'DisplayName', 'none')
plot_circle(K.ri, -K.mu, 0, 'Parking Orbit', [0.4660, 0.6740, 0.1880])
axis equal
axes('position',[.6 .3 .19 .19])
box on % put box around new pair of axes
hold on
plot(xx_simple(end-80:end, 1), xx_simple(end-80:end, 2),'Color', "#0072BD", 'DisplayName', 'none')
plot_circle(K.rf, 1-K.mu, 0, 'Arrival Orbit', [0.8500, 0.3250, 0.0980])
axis equal

% Plot the optimized trajectory in the ECI frame
figure
hold on
axis equal
grid minor
box on
plot(XX_simple(:, 1), XX_simple(:, 2),'Color', "#0072BD", 'DisplayName', 'S/C Trajectory') 
plot(XX_M(:, 1), XX_M(:, 2),'--', 'Color', "#D95319", 'DisplayName', 'Moon Orbit')
plot(XX_simple(1, 1), XX_simple(1, 2), '^', 'MarkerSize', 10, 'MarkerEdgeColor', "#0072BD", 'MarkerFaceColor', "#0072BD",'DisplayName', 'Starting Point')
plot(XX_simple(end, 1), XX_simple(end, 2), 'diamond', 'MarkerSize', 10,  'MarkerEdgeColor', "#0072BD", 'MarkerFaceColor', "#0072BD",'DisplayName', 'Arrirval Point')
plot(XX_M_0(1), XX_M_0(2), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'DisplayName', 'Initial Moon position')
xlabel('X [-]')  % Label X-axis
ylabel('Y [-]')  % Label Y-axis
legend('Location','northeast', 'Orientation','vertical')
%% Ex2.3: Multiple Shooting
K = constant();  % Define constants
Y0 = par2car(K.alpha, K.beta, K.ti,K.delta, K.ri, K.mu);  % Initial guess in state vector format
% Set ODE options
options = odeset('AbsTol', 2.5e-14, 'RelTol', 2.5e-14);

tt = linspace(Y0(5), Y0(6), 4);
[tt1, xx] = ode113(@(t, xx) PBRFBP_rhs(t, xx, K.mu, K.rho, K.oms, K.ms), ...
    tt, Y0(1:4), options);
xx = reshape(xx', 16, 1);
%Boundary conditions and initial guess for optimization (Y0_multiple)
Y0_multiple = [xx; Y0(5); Y0(6)];

% Define the linear inequality constraint: A*Y <= b
A = zeros(1, 18);  % Constraint ensures that t0 - tf = 0
A(1, 17) = 1;  
A(1, 18) = -1;
b  = 0;  % b = 0 for the linear inequality constraint

% Set optimization options for fmincon (for multiple shooting)

% Here, the commented options are those used to obtain the graphs in point 3. 
% They were commented out due to the computational cost of multiple shooting.
% The uncommented options (default MaxIterations and MaxFunctionEvaluations) 
% are those used to obtain the results presented in APPENDIX F of the report.
% To visualize the results presented in the report, simply uncomment the 
% currently commented options.

options = optimoptions('fmincon', ...
    'Display', 'iter-detailed', ...              
    'Algorithm', 'active-set', ...              
    'SpecifyObjectiveGradient', true, ...     
    'SpecifyConstraintGradient', true, ...     
    'ConstraintTolerance', 1e-10);

% options = optimoptions('fmincon', ...
%     'Display', 'iter-detailed', ...              
%     'Algorithm', 'active-set', ...              
%     'SpecifyObjectiveGradient', true, ...     
%     'SpecifyConstraintGradient', true, ...     
%     'ConstraintTolerance', 1e-10, ...
%     'MaxIterations', 1e4, ...
%     'MaxFunctionEvaluations', 4e4);

% Solve the optimization problem using fmincon
[Y_multiple, dv_multiple, exitflag_multiple,output_multiple,~,~,hessian_multiple] = fmincon( ...
    @obj_function_multiple, Y0_multiple, A, b, [], [], [], [], @constrains_multiple, options);

Iter    = output_multiple.iterations;
FunEval = output_multiple.funcCount;
MaxCon  = output_multiple.constrviolation;
dt = (Y_multiple(end)-Y_multiple(end-1))*K.TU;
dv = dv_multiple*K.VU;

Deltav_multiple = dv*1e3;
Deltat_multiple = dt;

components_out = {'dv [km/s]', 'dt [days]', 'iterations', 'funeval', 'MaxCon'}';
vec_out = [dv; dt; Iter; FunEval; MaxCon];
tab_out = table(components_out, vec_out, 'VariableNames', {'Parameter', 'Value'});
disp(tab_out)

% Define the component names
components = {'x', 'y', 'v_x', 'v_y', 't_i', 't_f'}';

% Evaluate constraints for the optimal solution
% Evaluate the constraints for the optimized solution
[g, c] = constrains_multiple(Y_multiple);
% Create the table
Y_multiple_table = table(components, [Y_multiple(1:4);Y_multiple(17:18)], 'VariableNames', {'Component', 'Value'});

% Display the table
disp(Y_multiple_table);

% Compute and display additional outputs
disp(['Delta-v: ', num2str(dv_multiple)]);
disp(['Exit flag: ', num2str(exitflag_multiple)]);
disp(['Minimum eigenvalue of Hessian: ', num2str(min(eig(hessian_multiple)))]);

% Display the results
disp('Inequality constraints (c):');
disp(g);

disp('Equality constraints (ceq):');
disp(c);

K = constant();
xx1 = Y_multiple(1:4);
xx2 = Y_multiple(5:8);
xx3 = Y_multiple(9:12);
xx4 = Y_multiple(13:16);
t1  = Y_multiple(17);
tN  = Y_multiple(18);

t_grid = linspace(t1, tN, 4).';
xx_nodes = [xx1,xx2,xx3,xx4]';
XX_nodes = rot2ECI(xx_nodes, t_grid, K.mu);

% Solve the system again using ode113 (with optimized states and time spans)
xx0 = Y_multiple(1:4);  % Optimized initial state
tspan = Y_multiple(end-1:end);  % Optimized time span

% Set ODE solver options with high accuracy
options = odeset('RelTol', 2.5e-14, 'AbsTol', 2.5e-14);

% Solve the ODE system for the optimal trajectory
[tt, xx] = ode113(@(t, xx) PBRFBP_rhs(t, xx, K.mu, K.rho, K.oms, K.ms), tspan, xx0, options);

% Convert trajectory from Synodic to ECI frame
XX = rot2ECI(xx, tt, K.mu);
xx_M = zeros(length(tt), 4);
xx_M(:,1) = 1-K.mu;
XX_M = rot2ECI(xx_M, tt, K.mu);
XX_M_0 = XX_M(1,:);

% Plot the optimized trajectory in the synodic frame
figure
hold on
axis equal
grid minor
box on
plot(xx(:, 1), xx(:, 2),'Color', "#0072BD", 'DisplayName', 'S/C Trajectory')  % Plot x-y trajectory in synodic frame
plot_circle(K.ri, -K.mu, 0, 'Parking Orbit', [0.4660, 0.6740, 0.1880])
plot_circle(K.rf, 1-K.mu, 0, 'Arrival Orbit', [0.8500, 0.3250, 0.0980])
for k = 1 : size(xx_nodes, 1)
    plot(xx_nodes(k,1), xx_nodes(k,2),'o')
end
xlabel('x [-]')  % Label x-axis
ylabel('y [-]')  % Label y-axis
legend('S/C Trajectory', 'Parking Orbit', 'Arrival Orbit', 'Nodes')
legend('Location','northwest', 'Orientation','vertical')
axes('position',[.3 .2 .19 .19])
box on % put box around new pair of axes
hold on
plot(xx(1:150, 1), xx(1:150, 2),'Color', "#0072BD", 'DisplayName', 'none')
plot_circle(K.ri, -K.mu, 0, 'Parking Orbit', [0.4660, 0.6740, 0.1880])
axis equal
axes('position',[.6 .2 .19 .19])
box on % put box around new pair of axes
hold on
plot(xx(end-80:end, 1), xx(end-80:end, 2),'Color', "#0072BD", 'DisplayName', 'none')
plot_circle(K.rf, 1-K.mu, 0, 'Arrival Orbit', [0.8500, 0.3250, 0.0980])
axis equal

% Plot the optimized trajectory in the ECI frame
figure
hold on
axis equal
grid minor
box on
plot(XX(:, 1), XX(:, 2),'Color', "#0072BD", 'DisplayName', 'S/C Trajectory') 
plot(XX_M(:, 1), XX_M(:, 2),'--', 'Color', "#D95319", 'DisplayName', 'Moon Orbit')
plot(XX(1, 1), XX(1, 2), '^', 'MarkerSize', 10, 'MarkerEdgeColor', "#0072BD", 'MarkerFaceColor', "#0072BD",'DisplayName', 'Starting Point')
plot(XX(end, 1), XX(end, 2), 'diamond', 'MarkerSize', 10,  'MarkerEdgeColor', "#0072BD", 'MarkerFaceColor', "#0072BD",'DisplayName', 'Arirval Point')
plot(XX_M_0(1), XX_M_0(2), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'DisplayName', 'initial Moon position')
for k = 1 : size(XX_nodes, 1)
    plot(XX_nodes(k,1), XX_nodes(k,2),'o')
end
xlabel('X [-]')  % Label X-axis
ylabel('Y [-]')  % Label Y-axis
legend('S/C Trajectory', 'Moon Orbit', 'Starting Point', 'Arrirval Point', 'Initial Moon position', 'Nodes')
legend('Location','northeast', 'Orientation','vertical')

Im = flip(imread('grafico_pareto_effic.png'));
figure
hold on
image([0,100], [3800,4150], Im)
plot(Deltat_simple, Deltav_simple, 'o', 'MarkerSize', 15, 'MarkerEdgeColor','m', 'MarkerFaceColor','m')
plot(Deltat_multiple, Deltav_multiple, 'o', 'MarkerSize', 15, 'MarkerEdgeColor','g', 'MarkerFaceColor','g')
ax = gca; 
ax.FontSize = 25;
xlim([0,100])
ylim([3800,4150])
yticks([0 3800 3850 3900 3950 4000 4050 4100 4150])
xlabel('$\Delta t$ (days)', 'FontSize', 30)
ylabel('$\Delta v$ (m/s)', 'FontSize', 30)
leg = legend('Simple Shooting', 'Multiple Shooting');
leg.Location = 'northeast';
leg.Orientation = "vertical";
leg.FontSize = 30;
%% Ex 2.4: n-body propagation
clc

% Print the number of loaded kernels for verification
fprintf('Number of LSK  kernels: %d\n', cspice_ktotal('lsk'));
fprintf('Number of SPK  kernels: %d\n', cspice_ktotal('spk'));
fprintf('Number of PCK  kernels: %d\n', cspice_ktotal('pck'));
fprintf('Number of CK   kernels: %d\n', cspice_ktotal('ck'));
fprintf('Number of TEXT kernels: %d\n', cspice_ktotal('TEXT'));
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));

% Define the initial epoch and target orbital angle
ti = Y_simple(end-1);  % Final time from the simple shooting method (Y_simple)
theta_target = wrapTo2Pi(K.oms * ti);  % Target angle for the orbit

% Define the list of celestial bodies to be considered in the n-body simulation
labels = {'Sun'; 'Mercury'; 'Venus'; 'Earth'; 'Moon'; 'Mars Barycenter'; ...
          'Jupiter Barycenter'; 'Saturn Barycenter'; 'Uranus Barycenter'; ...
          'Neptune Barycenter'; 'Pluto Barycenter'};
      
% Initialize the n-body simulation for the given bodies
bodies = nbody_init(labels);

% Define the center and reference frame for the simulation
center = 'Earth';  % Center the simulation on Earth
frame = 'ECLIPJ2000';  % Use the Ecliptic frame for reference (different from J2000)

% Define the initial epoch as a string and convert it to ephemeris time
ref_epoch_str = '2024-Sep-28 00:00:00.0000 TDB';  % Reference epoch in TDB
et0 = cspice_str2et(ref_epoch_str);  % Convert to ephemeris time

% Convert the initial state to the ECI frame and scale to correct units
X0 = rot2ECI(Y_simple(1:4)', ti, K.mu);  % Transform to the ECI frame
X0 = [X0(1:2)' * K.DU; 0; X0(3:4)' * K.VU; 0];  % Scale units to km and km/s
disp('X0_nbody:')
disp(X0)

% Define the final epoch as a string and convert it to ephemeris time
final_epoch_str = '2024-Oct-25 00:00:00.0000 TDB';  % Final epoch
etf = cspice_str2et(final_epoch_str);  % Convert to ephemeris time

% Generate a time vector between initial and final epochs
et_vec = linspace(et0, etf, 1000);  % 1000 time steps between initial and final epochs
theta_vec = nan(size(et_vec));  % Initialize the vector to store angle values

% Loop over each time step to compute the Sun's angle relative to Earth-Moon system
for k = 1 : length(et_vec)
    et = et_vec(k);  % Get current ephemeris time

    % Get the state vectors (position and velocity) of the Moon and Sun
    rv_moon = cspice_spkezr('Moon', et, frame, 'NONE', 'EMB');  % Moon state vector
    rv_sun = cspice_spkezr('Sun', et, frame, 'NONE', 'EMB');    % Sun state vector

    % Normalize the Moon's position and velocity
    rm = rv_moon(1:3) / norm(rv_moon(1:3));
    vm = rv_moon(4:6) / norm(rv_moon(4:6));

    % Compute the unit vectors for the reference frame
    u1 = rm;  % Moon position unit vector
    u3 = cross(rm, vm) / norm(cross(rm, vm));  % Perpendicular unit vector
    u2 = cross(u3, u1) / norm(cross(u3, u1));  % Second perpendicular unit vector
    U = [u1, u2, u3].';  % Rotation matrix for the reference frame

    % Rotate the Sun's position to the new frame
    r_sun_rot = U * rv_sun(1:3);
    
    % Compute the angle of the Sun in the new frame (ecliptic plane)
    theta = wrapTo2Pi(atan2(r_sun_rot(2), r_sun_rot(1)));
    
    % Store the angle in the theta_vec
    theta_vec(k) = theta;
end

threshold = pi / 4;  % Example threshold for a discontinuity
theta_diff = abs(diff(theta_vec));  % Compute differences between consecutive angles
discontinuity_indices = find(theta_diff > threshold);  % Find where differences exceed threshold
theta_vec_fixed = theta_vec;
theta_vec_fixed(discontinuity_indices) = NaN;
% Plot the computed angle (Sun's angle) over time
figure
hold on
plot(et_vec, theta_vec_fixed)  % Plot angle over time
yline(theta_target)  % Plot the target angle as a horizontal line
xlabel('Time (TDB)')  % Time axis
ylabel('Angle (radians)')  % Angle axis

% Now, we solve for the time when the angle reaches the target value using fsolve
fun = @(et) initial_epoch(et) - theta_target;  % Define the function for fsolve
options = optimoptions('fsolve', 'Display', 'iter', 'OptimalityTolerance', 1e-13);  % Solver options
x0 = 7.82e8;  % Initial guess for the target ephemeris time
et_targ = fsolve(fun, x0, options);  % Solve for the target epoch

% Convert the target time to UTC and print it
cspice_et2utc(et_targ, 'C', 10)
t_plot = datetime(cspice_timout(et_vec,'YYYY-MM-DD HR:MN:SC.###'),...
    'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
t_targ = datetime(cspice_timout(et_targ,'YYYY-MM-DD HR:MN:SC.###'),...
    'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
% Plot again with the target time and epoch line
theta_vec = rad2deg(theta_vec);
dth = diff(theta_vec);
theta_vec(abs(dth)>10) = NaN;
figure
hold on
plot(t_plot, theta_vec, "Color",[0, 0.4470, 0.7410], 'DisplayName', '$\theta$')  % Plot angle over time
yline(rad2deg(theta_target), '--', 'LineWidth', 1.5, "Color",[0.8500, 0.3250, 0.0980], 'DisplayName', '$\theta_{target}$')  % Plot target angle
xline(t_targ, '--', 'LineWidth', 1.5, "Color",[0.4660, 0.6740, 0.1880], 'DisplayName', '$t_i$')  % Plot the target time
ylabel('Angle [deg]')
legend

% Set the initial and final times for n-body propagation
ti = et_targ;  % Start time of the n-body propagation
tf = et_targ + (Y_simple(6) - Y_simple(5)) * K.TU * 24 * 3600;  % End time, adjusted to days
cspice_et2utc(tf, 'C', 10)

% Define solver options and perform n-body propagation
options = odeset('RelTol', 2.5e-14, 'AbsTol', 2.5e-14);  % Solver options for accuracy
[tt_nbody, xx_nbody] = ode113(@(t, xx) nbody_shift_rhs(t, xx, bodies, frame, center), ...
    [ti, tf], X0, options);  % Solve the n-body ODE

% Convert the previously computed PBRFBP trajectory to km for comparison
RR_simple = [XX_simple(:, 1) * K.DU, XX_simple(:, 2) * K.DU, zeros(size(XX_simple, 1), 1)];

% Plot the results: compare the n-body propagation with the PBRFBP results
figure
hold on
grid on
plot3(xx_nbody(:, 1), xx_nbody(:, 2), xx_nbody(:, 3),'Color', "#0072BD", 'DisplayName', 'n-body propagation')  % Plot n-body trajectory
plot3(RR_simple(:, 1), RR_simple(:, 2), RR_simple(:, 3),'Color', "#D95319", 'DisplayName', 'PBRFBP')  % Plot PBRFBP trajectory
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
legend('n-body propagation', 'PBRFBP')
view(0, 90)

figure
hold on
grid on
plot3(xx_nbody(:, 1), xx_nbody(:, 2), xx_nbody(:, 3),'Color', "#0072BD", 'DisplayName', 'n-body propagation')  % Plot n-body trajectory
plot3(RR_simple(:, 1), RR_simple(:, 2), RR_simple(:, 3),'Color', "#D95319", 'DisplayName', 'PBRFBP')  % Plot PBRFBP trajectory
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
legend('n-body propagation', 'PBRFBP')
view(15, 15)

%%

function constants = constant()
% Defines and returns a structure containing constants relevant to 
% the Planar-Bicircular Restricted Four-Body Problem (PBRFBP).
%
% Outputs:
%   constants - A structure containing the following fields:
%       Re   - Earth's radius (km)
%       Rm   - Moon's radius (km)
%       hi   - Initial altitude above Earth's surface (km)
%       hf   - Final altitude above Moon's surface (km)
%       DU   - Characteristic distance unit (Earth-Moon distance, km)
%       TU   - Characteristic time unit
%       VU   - Characteristic velocity unit
%       ri   - Initial normalized radial distance from Earth
%       rf   - Final normalized radial distance from Moon
%       mu   - Gravitational parameter ratio of the Earth-Moon system
%       oms  - Rotational velocity of the system (rad/TU)
%       rho  - Radius of the circular orbit of the third body (km)
%       ms   - Mass parameter for the secondary body
%       alpha - Initial angle of the trajectory (radians)
%       beta  - Velocity scaling factor
%       ti    - Initial time
%       delta - Time span of interest
%       T     - Period of the motion
%       vi    - Initial normalized velocity
%       vf    - Final normalized velocity

% Initialize the constants structure
constants = struct();

% Physical and problem-specific parameters
radii_e       = cspice_bodvrd('Earth', 'RADII', 3);
radii_m       = cspice_bodvrd('Moon', 'RADII', 3);
constants.Re  = radii_e(1);  % Earth's radius (km)
constants.Rm  = radii_m(1);   % Moon's radius (km)
constants.hi  = 167;                        % Initial altitude above Earth's surface (km)
constants.hf  = 100;                        % Final altitude above Moon's surface (km)
constants.DU  = 3.84405000e5;               % Characteristic distance unit: Earth-Moon distance (km)
constants.TU  = 4.34256461;                 % Characteristic time unit (Earth-Moon system)
constants.VU  = 1.02454018;                 % Characteristic velocity unit (normalized units)

% Compute normalized radial distances
constants.ri  = (constants.Re + constants.hi) / constants.DU;  % Initial distance (normalized)
constants.rf  = (constants.Rm + constants.hf) / constants.DU;  % Final distance (normalized)

% Gravitational and orbital parameters
mu_e = cspice_bodvrd('Earth', 'GM', 1);
mu_m = cspice_bodvrd('Moon', 'GM', 1);
constants.mu = mu_m/(mu_m+mu_e);
%constants.mu  = 0.012150584269940;  % Gravitational parameter ratio (Moon/Earth+Moon)
constants.oms = -9.25195985e-1;     % Rotational velocity (rad/TU)
constants.rho = 3.88811143e2;       % Radius of circular orbit of the third body (km)
constants.ms  = 3.28900541e5;       % Mass parameter of secondary body

% Trajectory parameters
constants.alpha = 0.2 * pi;         % Initial angle of trajectory (radians)
constants.beta  = 1.41;             % Velocity scaling factor
constants.ti    = 2;                % Initial time (normalized units)
constants.delta = 4;                % Time duration of interest (normalized units)
constants.T     = 23;               % Period of the motion (normalized units)

% Compute initial and final normalized velocities
constants.vi    = sqrt((1 - constants.mu) / constants.ri);  % Initial velocity (normalized)
constants.vf    = sqrt(constants.mu / constants.rf);        % Final velocity (normalized)

end

function Y = par2car(alpha, beta, ti, delta, ri, mu)
% Converts polar orbital parameters into Cartesian state vectors 
% in the rotating reference frame.
%
% Inputs:
%   alpha - Initial angular position (in radians)
%   beta  - Scaling factor for the initial velocity
%   ti    - Initial time
%   delta - Time duration for propagation
%   ri    - Initial radial distance from the primary body
%   mu    - Gravitational parameter ratio for the system
%
% Outputs:
%   Y - Combined state vector and time span:
%       [x0; y0; vx0; vy0; ti; ti + delta]
%       where (x0, y0) are initial positions, 
%       (vx0, vy0) are initial velocities, and 
%       ti and ti+delta are the time bounds.

% Initial radial position
r0 = ri;

% Initial velocity magnitude based on the provided scaling factor
v0  = beta * sqrt((1-mu)/r0);

% Convert polar coordinates to Cartesian in the rotating frame
x0  = r0 * cos(alpha) - mu;      % Initial x-position
y0  = r0 * sin(alpha);           % Initial y-position
vx0 = -(v0 - r0) * sin(alpha);   % Initial x-velocity
vy0 = (v0 - r0)  * cos(alpha);   % Initial y-velocity

% Assemble initial state vector
xx0   = [x0; y0; vx0; vy0];

% Time span for the solution
tspan = [ti; ti + delta];

% Combine state vector and time span
Y = [xx0; tspan];

end

function XX = rot2ECI(xx, t, mu)
% Converts state vectors from the rotating frame to the Earth-Centered Inertial (ECI) frame.
%
% The ECI frame is stationary relative to the primary body, while 
% the rotating frame moves with the system. This function uses a 
% time-dependent rotation matrix for the transformation.
%
% Inputs:
%   xx - State vectors in the rotating frame (N x 4 matrix):
%        [x, y, vx, vy] for each row
%   t  - Time(s) at which to perform the transformation (scalar or vector)
%   mu - Gravitational parameter ratio for the system
%
% Outputs:
%   XX - State vectors in the ECI frame (N x 4 matrix):
%        [X, Y, Vx, Vy] for each row

% Preallocate output matrix for ECI state vectors
XX = nan(size(xx));

% Compute rotation matrix components
Cos = cos(t); % Cosine of the rotation angle
Sin = sin(t); % Sine of the rotation angle

% Apply the rotation transformation for positions
XX(:, 1) = (xx(:, 1) + mu) .* Cos - xx(:, 2) .* Sin; % X = x*cos - y*sin
XX(:, 2) = (xx(:, 1) + mu) .* Sin + xx(:, 2) .* Cos; % Y = x*sin + y*cos

% Apply the rotation transformation for velocities
XX(:, 3) = (xx(:, 3) - xx(:, 2)) .* Cos - (xx(:, 4) + xx(:, 1) + mu) .* Sin; % Vx
XX(:, 4) = (xx(:, 3) - xx(:, 2)) .* Sin + (xx(:, 4) + xx(:, 1) + mu) .* Cos; % Vy

end

function [dxdt] = PBRFBP_rhs(t,xx, mu, rho, oms, ms)
% This function computes the right-hand side (RHS) of the equations of motion 
% for a Planar-Bicircular Restricted Four-Body Problem (PBRFBP).
% Inputs:
%   t   - Current time
%   xx  - State vector [x; y; vx; vy], where:
%         x, y   - Position coordinates in the rotating reference frame
%         vx, vy - Velocity components in the rotating reference frame
%   mu  - Gravitational parameter ratio for the two main bodies
%   rho - Distance of the third body (secondary body) from the origin
%   oms - Angular velocity of the secondary body in the rotating frame
%   ms  - Mass parameter for the secondary body
% Outputs:
%   dxdt - Time derivative of the state vector [dx/dt; dy/dt; dvx/dt; dvy/dt]

% Extract variables from the state vector
x  = xx(1); % x-position
y  = xx(2); % y-position
vx = xx(3); % x-velocity
vy = xx(4); % y-velocity

% Compute distances to the primary and secondary bodies
r1 = sqrt((x + mu)^2 + y^2);                % Distance to the first primary body
r2 = sqrt((x + mu - 1)^2 + y^2);            % Distance to the second primary body
r3 = sqrt((x - rho*cos(oms*t))^2 + ...      % Distance to the secondary body
          (y - rho*sin(oms*t))^2);

% Compute partial derivatives of the effective potential (gravitational forces)
dUdx = x - (1-mu)/r1^3*(mu + x) + ...       % Gravitational influence of body 1
        mu/r2^3*(1 - mu - x) - ...          % Gravitational influence of body 2
        ms*cos(oms*t)/rho^2 - ...           % Influence from secondary body's motion
        ms*(x - rho*cos(oms*t))/r3^3;       % Gravitational effect of the secondary body

dUdy = y - (1-mu)/r1^3*y - ...              % Gravitational influence of body 1
        mu/r2^3*y - ...                     % Gravitational influence of body 2
        ms*sin(oms*t)/rho^2 - ...           % Influence from secondary body's motion
        ms*(y - rho*sin(oms*t))/r3^3;       % Gravitational effect of the secondary body

% Alternative version (if the secondary body's influence is ignored)
% Uncomment these lines and comment out the above `dUdx` and `dUdy` if secondary body effects are not considered.
% dUdx = x - (1-mu)/r1^3*(mu + x) + mu/r2^3*(1 - mu - x);
% dUdy = y - (1-mu)/r1^3*y - mu/r2^3*y;

% Assemble the right-hand side of the equations of motion
dxdt = zeros(4,1);      % Initialize the derivative vector
dxdt(1:2) = xx(3:4);    % dx/dt = vx, dy/dt = vy
dxdt(3)   = dUdx + 2*vy; % dvx/dt = dUdx + Coriolis force term (2*vy)
dxdt(4)   = dUdy - 2*vx; % dvy/dt = dUdy - Coriolis force term (2*vx)

end

function [dXdt]  = PBRFBP_STM(t, X, mu, rho, oms, ms)
% This function computes the time derivative of the state vector and 
% State Transition Matrix (STM) for the Planar-Bicircular Restricted Four-Body Problem (PBRFBP).
%
% The PBRFBP describes the motion of a small particle under the gravitational 
% influence of three primary bodies (two main bodies and a secondary body) 
% in a rotating reference frame, with the secondary body's motion 
% approximated as circular.
%
% Inputs:
%   t   - Current time
%   X   - State vector augmented with the flattened State Transition Matrix (STM):
%         [x; y; vx; vy; PHI(:)], where PHI is the 4x4 STM.
%   mu  - Gravitational parameter ratio for the two main bodies
%   rho - Distance of the secondary body from the origin
%   oms - Angular velocity of the secondary body in the rotating reference frame
%   ms  - Mass parameter for the secondary body
%
% Outputs:
%   dXdt - Time derivative of the augmented state vector [dx/dt; dy/dt; dvx/dt; dvy/dt; dPHIdt(:)]

% Extract state variables from the input vector
x  = X(1); % x-position
y  = X(2); % y-position
vx = X(3); % x-velocity
vy = X(4); % y-velocity

% Reshape the flattened State Transition Matrix (STM) into a 4x4 matrix
PHI = reshape(X(5:end), 4, 4);

% Compute distances to the primary and secondary bodies
r1 = sqrt((x + mu)^2 + y^2);                % Distance to the first primary body
r2 = sqrt((x + mu - 1)^2 + y^2);            % Distance to the second primary body
r3 = sqrt((x - rho*cos(oms*t))^2 + ...      % Distance to the secondary body
          (y - rho*sin(oms*t))^2);

% Compute partial derivatives of the effective potential (gravitational forces)
dUdx = x - (1-mu)/r1^3*(mu+x) + ...         % Contribution from the first primary body
        mu/r2^3*(1-mu-x) - ...              % Contribution from the second primary body
        ms*cos(oms*t)/rho^2 - ...           % Secondary body's influence in the x-direction
        ms*(x - rho*cos(oms*t))/r3^3;       % Gravitational effect of the secondary body

dUdy = y - (1-mu)/r1^3*y - ...              % Contribution from the first primary body
        mu/r2^3*y - ...                     % Contribution from the second primary body
        ms*sin(oms*t)/rho^2 - ...           % Secondary body's influence in the y-direction
        ms*(y - rho*sin(oms*t))/r3^3;       % Gravitational effect of the secondary body

% Initialize the time derivative of the augmented state vector
dXdt = zeros(20,1);

% Derivatives of position (x, y) and velocity (vx, vy)
dXdt(1:2) = X(3:4);           % dx/dt = vx, dy/dt = vy
dXdt(3)   = dUdx + 2*vy;      % dvx/dt = dUdx + Coriolis term (2*vy)
dXdt(4)   = dUdy - 2*vx;      % dvy/dt = dUdy - Coriolis term (2*vx)

% Compute the components of the Jacobian matrix (N)
% N is the Hessian of the effective potential
N1 = (mu-1)/r1^3 - mu/r2^3 - ms/r3^3 - 3*((x + mu)^2)*(mu-1)/r1^5 + ...
     (3*ms*(x - rho*cos(oms*t))^2)/r3^5 + (3*mu*(mu + x - 1)^2)/r2^5 + 1;

N2 = 3*ms*(x - rho*cos(oms*t))*(y - rho*sin(oms*t))/r3^5 + ...
     3*mu*y*(x + mu - 1)/r2^5 - 3*y*(mu - 1)*(x + mu)/r1^5;

N3 = N2; % N is symmetric, so N3 = N2

N4 = (mu-1)/r1^3 - mu/r2^3 - ms/r3^3 - 3*y^2*(mu-1)/r1^5 + ...
     (3*ms*(y - rho*sin(oms*t))^2)/r3^5 + (3*mu*y^2)/r2^5 + 1;

% Assemble the 2x2 Jacobian matrix of the potential
N = [N1, N2; N3, N4];

% Full 4x4 Jacobian matrix of the system dynamics
dfdx = [zeros(2), eye(2); ...                % Velocity derivatives
        N, [0, 2; -2, 0]];                  % Position derivatives and Coriolis terms

% Compute the time derivative of the State Transition Matrix (STM)
dPHIdt = dfdx * PHI;                        % STM evolution equation: dPHI/dt = A * PHI

% Flatten the STM derivative and append it to the state derivatives
dXdt(5:end) = dPHIdt(:);

end

function [dv, G] = obj_function(Y)
% Computes the objective function for the trajectory optimization in the
% Planar-Bicircular Restricted Four-Body Problem (PBRFBP). The function 
% calculates the total delta-v (velocity change) and optionally returns
% the gradient of the objective function.
%
% Inputs:
%   Y - Initial state vector and time span [x0; y0; vx0; vy0; ti; tf]
%
% Outputs:
%   dv - Total delta-v (sum of initial and final velocity changes)
%   G  - Gradient of the objective function with respect to Y (optional)
%
% Uses:
%   constant()    - Function defining problem-specific constants.
%   propagate()   - Function to propagate the trajectory and compute final state.
%   gradient_f()  - Function to compute the gradient of the objective.

% Retrieve constants
K = constant();

% Extract initial state from input vector
xx0 = Y(1:4);  % Initial state: [x0; y0; vx0; vy0]

% Propagate trajectory to compute the final state and state transition matrix
[xxf, PHIf, ~, ~, ~] = propagate(Y);

% Compute the initial velocity magnitude and required initial velocity
v0_actual = sqrt((xx0(3) - xx0(2))^2 + (xx0(4) + xx0(1) + K.mu)^2);
v0_required = sqrt((1 - K.mu) / K.ri);

% Compute the final velocity magnitude and required final velocity
vf_actual = sqrt((xxf(3) - xxf(2))^2 + (xxf(4) + xxf(1) + K.mu - 1)^2);
vf_required = sqrt(K.mu / K.rf);

% Calculate delta-v contributions (initial and final)
dvi = v0_actual - v0_required;  % Initial delta-v
dvf = vf_actual - vf_required;  % Final delta-v

% Total delta-v
dv = dvi + dvf;

% If gradient is requested, compute it
if nargout > 1
    G = gradient_f(Y, xxf, PHIf);
end

end

function [c, ceq, Gc, Gceq] = myconstrains(Y)
% Defines the nonlinear equality and inequality constraints for the
% optimization problem in the Planar-Bicircular Restricted Four-Body
% Problem (PBRFBP). It also computes the gradients of the constraints if
% requested.
%
% Inputs:
%   Y - Initial state vector and time span [x0; y0; vx0; vy0; ti; tf]
%
% Outputs:
%   c    - Inequality constraints (empty in this case, as there are none)
%   ceq  - Equality constraints
%   Gc   - Gradient of inequality constraints (empty in this case)
%   Gceq - Gradient of equality constraints with respect to Y (optional)
%
% Uses:
%   constant()    - Function defining problem-specific constants.
%   propagate()   - Function to propagate the trajectory and compute the final state.
%   gradient_c()  - Function to compute the gradient of the constraints.

% Retrieve constants
K = constant();

% Extract initial state from input vector
xx0 = Y(1:4);  % Initial state: [x0; y0; vx0; vy0]

% Propagate trajectory to compute the final state and state transition matrix
[xxf, PHIf, ~, ~, ~] = propagate(Y);

% Extract initial and final positions and velocities
x0 = xx0(1); y0 = xx0(2); vx0 = xx0(3); vy0 = xx0(4);
xf = xxf(1); yf = xxf(2); vxf = xxf(3); vyf = xxf(4);

% Initial equality constraints (boundary conditions at the initial point):
%   - Distance constraint: Initial position should lie on the initial orbit radius.
%   - Velocity constraint: Initial velocity should align with the orbit.
psi_i = [(x0 + K.mu)^2 + y0^2 - K.ri^2;          % Distance constraint
         (x0 + K.mu) * (vx0 - y0) + y0 * (vy0 + x0 + K.mu)];  % Velocity constraint

% Final equality constraints (boundary conditions at the final point):
%   - Distance constraint: Final position should lie on the final orbit radius.
%   - Velocity constraint: Final velocity should align with the orbit.
psi_f = [(xf + K.mu - 1)^2 + yf^2 - K.rf^2;       % Distance constraint
         (xf + K.mu - 1) * (vxf - yf) + yf * (vyf + xf + K.mu - 1)];  % Velocity constraint

% Inequality constraints (none in this case)
c = [];

% Equality constraints (combine initial and final conditions)
ceq = [psi_i; psi_f];

% If gradients are requested, compute them
if nargout > 2
    Gc = [];  % Gradient of inequality constraints (empty)
    Gceq = gradient_c(Y, xxf, PHIf);  % Gradient of equality constraints
end

end

function [xf, PHIf, tf, X, tt] = propagate(Y)
% Propagates the state and computes the final state and state transition matrix (STM)
% for the Planar-Bicircular Restricted Four-Body Problem (PBRFBP).
%
% Inputs:
%   Y   - State and time span vector: [x0; y0; vx0; vy0; ti; tf]
%
% Outputs:
%   xf    - Final state vector: [x; y; vx; vy] at time tf
%   PHIf  - Final state transition matrix (4x4) at time tf
%   tf    - Final time of propagation
%   X     - Full state trajectory including STM over the integration time span
%   tt    - Time vector corresponding to the trajectory
%
% Uses:
%   constant()    - Function defining problem-specific constants.
%   PBRFBP_STM()  - Function defining the dynamics and STM evolution.

% Retrieve constants
K = constant();

% Extract the initial state and time span
x0 = Y(1:4);  % Initial state: [x0; y0; vx0; vy0]
ti = Y(5);    % Initial time
tf = Y(6);    % Final time
tspan = [ti, tf];  % Time span for propagation

% Initial state transition matrix (STM) is the identity matrix
Phi0 = eye(4);

% Combine the state and the flattened STM into a single vector for integration
X0 = [x0; Phi0(:)];

% Set ODE solver options for high accuracy
options = odeset('RelTol', 2.5e-14, 'AbsTol', 2.5e-14);

% Integrate the dynamics and STM over the time span using ode78
[tt, X] = ode113(@(t, X) PBRFBP_STM(t, X, K.mu, K.rho, K.oms, K.ms), tspan, X0, options);

% Extract the final state from the last row of the integrated solution
xf = X(end, 1:4)';  % Final state: [x; y; vx; vy]

% Extract the final STM by reshaping the relevant part of the last row
PHIf = reshape(X(end, 5:end), 4, 4);  % Final STM (4x4 matrix)

end

function G = gradient_f(Y, xxf, PHI)
% Computes the gradient of the objective function with respect to the decision variables
% in the Planar-Bicircular Restricted Four-Body Problem (PBRFBP).
%
% Inputs:
%   Y    - State and time vector: [x0; y0; vx0; vy0; ti; tf]
%   xxf  - Final state vector: [xf; yf; vxf; vyf]
%   PHI  - State transition matrix (STM) at the final time
%
% Outputs:
%   G    - Gradient of the objective function with respect to Y

% Retrieve problem constants
K = constant();

% Extract initial and final state components
x0 = Y(1);  % Initial x position
y0 = Y(2);  % Initial y position
vx0 = Y(3); % Initial x velocity
vy0 = Y(4); % Initial y velocity
xf = xxf(1); % Final x position
yf = xxf(2); % Final y position
vxf = xxf(3); % Final x velocity
vyf = xxf(4); % Final y velocity

% Compute the distances for the initial and final states
rad1 = sqrt((vx0 - y0)^2 + (vy0 + x0 + K.mu)^2); % Distance at initial time
rad2 = sqrt((vxf - yf)^2 + (vyf + xf + K.mu - 1)^2); % Distance at final time

% Compute the partial derivatives of the objective function with respect to initial and final states
dvidxi = 1 / rad1 * [x0 + vy0 + K.mu; y0 - vx0; vx0 - y0; x0 + vy0 + K.mu];  % Gradient wrt initial state
dvfdxf = 1 / rad2 * [xf + vyf + K.mu - 1; yf - vxf; vxf - yf; xf + vyf + K.mu - 1]; % Gradient wrt final state

% Compute the right-hand side dynamics for the initial and final times
fi = PBRFBP_rhs(Y(5), Y(1:4), K.mu, K.rho, K.oms, K.ms);  % Initial time RHS (dynamics)
ff = PBRFBP_rhs(Y(6), xxf, K.mu, K.rho, K.oms, K.ms);    % Final time RHS (dynamics)

% Initialize the gradients with respect to time as zero
dvidti = 0;  % No dependence on initial time for dvi
dvidtf = 0;  % No dependence on final time for dvi

% Combine gradients wrt initial state and time for the initial state gradient (dvidY)
dvidY = [dvidxi; dvidti; dvidtf];

% Compute the gradients with respect to final state and time
dvfdxi = PHI.' * dvfdxf;  % Gradient wrt initial state (final time)
dvfdti = -dvfdxf.' * PHI * fi;  % Gradient wrt initial time (final time)
dvfdtf = dvfdxf.' * ff;   % Gradient wrt final time (final time)

% Combine gradients with respect to final state and time (dvfdY)
dvfdY = [dvfdxi; dvfdti; dvfdtf];

% Total gradient is the sum of the initial and final state/time gradients
G = dvidY + dvfdY;

end

function Gc = gradient_c(Y, xxf, PHI)
% Computes the gradient of the constraint functions with respect to the decision variables.
%
% Inputs:
%   Y    - State and time vector: [x0; y0; vx0; vy0; ti; tf]
%   xxf  - Final state vector: [xf; yf; vxf; vyf]
%   PHI  - State transition matrix (STM) at the final time
%
% Outputs:
%   Gc   - Gradient of the constraint functions with respect to Y

% Initialize the gradient matrix
Gc = nan(6,4);

% Retrieve problem constants
K = constant();

% Extract initial and final state components from Y and xxf
x0 = Y(1);  % Initial x position
y0 = Y(2);  % Initial y position
vx0 = Y(3); % Initial x velocity
vy0 = Y(4); % Initial y velocity
xf = xxf(1); % Final x position
yf = xxf(2); % Final y position
vxf = xxf(3); % Final x velocity
vyf = xxf(4); % Final y velocity

% Compute the gradients of the constraint functions with respect to the initial state variables

% Constraint 1: (x0 + mu)^2 + y0^2 - r_i^2
dc1dxi = [2*(x0 + K.mu); 2*y0; 0; 0];

% Constraint 2: (x0 + mu)*(vx0 - y0) + y0*(vy0 + x0 + mu)
dc2dxi = [vx0; vy0; x0 + K.mu; y0];

% Compute the gradients of the constraint functions with respect to the final state variables

% Constraint 3: (xf + mu - 1)^2 + yf^2 - r_f^2
dc3dxf = [2*(xf + K.mu - 1); 2*yf; 0; 0];

% Constraint 4: (xf + mu - 1)*(vxf - yf) + yf*(vyf + xf + mu - 1)
dc4dxf = [vxf; vyf; xf + K.mu - 1; yf];

% Compute the right-hand side of the equations of motion for the initial time
xxidot = PBRFBP_rhs(Y(5), Y(1:4), K.mu, K.rho, K.oms, K.ms);

% Compute the gradient of the first and second constraint functions with respect to the decision variables
Gc(:, 1) = [dc1dxi; 0; 0];  % Gradient for constraint 1 wrt initial state (x0, y0, vx0, vy0)
Gc(:, 2) = [dc2dxi; 0; 0];  % Gradient for constraint 2 wrt initial state (x0, y0, vx0, vy0)

% Compute the right-hand side of the equations of motion for the final time
xxfdot = PBRFBP_rhs(Y(6), xxf, K.mu, K.rho, K.oms, K.ms);

% Compute the gradient of the third and fourth constraint functions with respect to the decision variables
Gc(:, 3) = [PHI.' * dc3dxf; -dc3dxf.' * PHI * xxidot; dc3dxf.' * xxfdot];  % Gradient for constraint 3 wrt final state (xf, yf, vxf, vyf)
Gc(:, 4) = [PHI.' * dc4dxf; -dc4dxf.' * PHI * xxidot; dc4dxf.' * xxfdot];  % Gradient for constraint 4 wrt final state (xf, yf, vxf, vyf)

end

function [dv, G] = obj_function_multiple(Y)
% Computes the objective function for the trajectory optimization in the
% Planar-Bicircular Restricted Four-Body Problem (PBRFBP) with multiple shooting.
% The function calculates the total delta-v (velocity change) for the initial 
% and final states and optionally returns the gradient of the objective function.
%
% Inputs:
%   Y - State vector containing the positions, velocities, and times of multiple
%       shooting points (for example, the initial and final states for each segment).
%
% Outputs:
%   dv - Total delta-v (sum of initial and final velocity changes)
%   G  - Gradient of the objective function with respect to Y (optional)
%
% Uses:
%   constant()    - Function defining problem-specific constants.
%   gradient_f_multiple() - Function to compute the gradient of the objective function.

% Retrieve system constants
K = constant();

% Extract initial state (first segment) and final state (last segment) from Y
xx1 = Y(1:4);     % Initial state for the first segment: [x1; y1; vx1; vy1]
xxN = Y(13:16);   % Final state for the last segment: [xN; yN; vxN; vyN]

% Compute the initial velocity magnitude and required initial velocity
v0_actual = sqrt((xx1(3) - xx1(2))^2 + (xx1(4) + xx1(1) + K.mu)^2);
v0_required = K.vi;  % Required initial velocity based on constant

% Compute the final velocity magnitude and required final velocity
vf_actual = sqrt((xxN(3) - xxN(2))^2 + (xxN(4) + xxN(1) + K.mu - 1)^2);
vf_required = K.vf;  % Required final velocity based on constant

% Calculate delta-v contributions (initial and final)
dvi = v0_actual - v0_required;  % Initial delta-v
dvf = vf_actual - vf_required;  % Final delta-v

% Total delta-v (sum of initial and final velocity changes)
dv = dvi + dvf;

% If gradient is requested, compute it
if nargout > 1
    G = gradient_f_multiple(Y);  % Compute gradient of the objective function
end

end

function G = gradient_f_multiple(Y)
% Computes the gradient of the objective function with respect to the state and time variables.
% This function calculates the derivatives of the velocity differences (dv) with respect to the initial 
% and final state variables for a multiple-shooting setup.
%
% Inputs:
%   Y - State vector containing the positions, velocities, and times of multiple
%       shooting points (for example, the initial and final states for each segment).
%
% Outputs:
%   G - The gradient of the objective function with respect to the decision variables (state and time).

% Retrieve system constants
K = constant();

% Extract the initial state for the first segment and final state for the last segment from Y
x1 = Y(1); y1 = Y(2); vx1 = Y(3); vy1 = Y(4);   % Initial state [x1, y1, vx1, vy1]
xN = Y(13); yN = Y(14); vxN = Y(15); vyN = Y(16); % Final state [xN, yN, vxN, vyN]

% Calculate the radial distances for the initial and final velocity differences
rad1 = sqrt((vx1 - y1)^2 + (vy1 + x1 + K.mu)^2);  % Initial point radial distance
rad2 = sqrt((vxN - yN)^2 + (vyN + xN + K.mu - 1)^2);  % Final point radial distance

% Compute the gradients of the objective function with respect to the initial state (xx1)
dv1dx1 = 1 / rad1 * [x1 + vy1 + K.mu; y1 - vx1; vx1 - y1; x1 + vy1 + K.mu];

% Compute the gradients of the objective function with respect to the final state (xxN)
dv2dxN = 1 / rad2 * [xN + vyN + K.mu - 1; yN - vxN; vxN - yN; xN + vyN + K.mu - 1];

% Initialize the gradient vector G as a zero vector of length 18
% This is because the vector Y contains 18 elements (positions, velocities, and time variables)
G = zeros(18, 1);

% Assign the computed gradients for the initial and final state variables to the appropriate indices in G
G(1:4)   = dv1dx1;   % Gradient for the initial state (xx1)
G(13:16) = dv2dxN;   % Gradient for the final state (xxN)

end

function [g, c, Gg, Gc] = constrains_multiple(Y)
% Computes the constraint functions and their gradients for the multiple shooting trajectory optimization problem.
% The function calculates both equality constraints (state matching) and inequality constraints 
% (proximity to the Earth and Moon), and returns their gradients if requested.
%
% Inputs:
%   Y - State vector containing positions, velocities, and times for multiple shooting points.
%
% Outputs:
%   g   - The inequality constraints (proximity to Earth and Moon).
%   c   - The equality constraints (state matching at shooting points).
%   Gg  - Gradient of the inequality constraints with respect to the state and time variables.
%   Gc  - Gradient of the equality constraints with respect to the state and time variables.

% Retrieve system constants
K = constant();

% Extract states and times for the multiple shooting segments
xx1 = Y(1:4); xx2 = Y(5:8); xx3 = Y(9:12); xx4 = Y(13:16);
t1 = Y(17); tN = Y(18);

% Compute time intervals between shooting points
h = (tN - t1) / 3;
t2 = t1 + h; t3 = t2 + h;

% Propagate the trajectory using the provided state vectors and times
[xf1, PHI12, ~, ~, ~]  = propagate([xx1; t1; t2]);
[xf2, PHI23, ~, ~, ~]  = propagate([xx2; t2; t3]);
[xf3, PHI34, ~, ~, ~]  = propagate([xx3; t3; tN]);

% Compute the state matching constraints
z1 = xf1 - xx2;
z2 = xf2 - xx3;
z3 = xf3 - xx4;

% Compute the initial and final state variables and constraints (Earth and Moon proximity)
x1 = xx1(1); y1 = xx1(2); vx1 = xx1(3); vy1 = xx1(4);
xN = xx4(1); yN = xx4(2); vxN = xx4(3); vyN = xx4(4);

% Compute the Earth and Moon proximity constraints
psi_i = [(x1 + K.mu)^2 + y1^2 - K.ri^2; (x1 + K.mu) * (vx1 - y1) + y1 * (vy1 + x1 + K.mu)];
psi_f = [(xN + K.mu - 1)^2 + yN^2 - K.rf^2; (xN + K.mu - 1) * (vxN - yN) + yN * (vyN + xN + K.mu - 1)];

% Combine the equality constraints (state matching and Earth/Moon proximity)
c = [z1; z2; z3; psi_i; psi_f];

% Compute the inequality constraints (proximity to Earth and Moon)
g = [eta(xx1); eta(xx2); eta(xx3); eta(xx4)];

% Compute the gradients of the constraints if requested
if nargout > 2
    Gc = gradient_c_multiple(Y, xf1, xf2, xf3, PHI12, PHI23, PHI34);
    Gg = gradient_g_multiple(Y);
end

end

function Gc = gradient_c_multiple(Y, xf1, xf2, xf3, PHI12, PHI23, PHI34)
% Computes the gradient of the equality constraint functions (state matching) for the multiple shooting
% trajectory optimization problem. The gradients are computed based on the propagated final states (xf1, xf2, xf3)
% and state transition matrices (PHI12, PHI23, PHI34).
%
% Inputs:
%   Y      - State vector containing positions, velocities, and times for multiple shooting points.
%   xf1, xf2, xf3  - Final states of each shooting segment, used for calculating the constraints.
%   PHI12, PHI23, PHI34 - State transition matrices between the shooting segments.
%
% Outputs:
%   Gc     - The gradient of the equality constraint functions with respect to the state and time variables.

% Retrieve system constants
K = constant();

% Extract initial and final state variables and times from the input vector
x1 = Y(1); y1 = Y(2); vx1 = Y(3); vy1 = Y(4);   % Initial state of first segment
xN = Y(13); yN = Y(14); vxN = Y(15); vyN = Y(16); % Final state of last segment
t1 = Y(17); tN = Y(18); h = (tN - t1)/3; % Time interval between shooting points
t2 = t1 + h; t3 = t2 + h;

% Calculate the state derivatives using the system dynamics
x1dot = PBRFBP_rhs(t1, Y(1:4), K.mu, K.rho, K.oms, K.ms);
x2dot = PBRFBP_rhs(t2, Y(5:8), K.mu, K.rho, K.oms, K.ms);
x3dot = PBRFBP_rhs(t3, Y(9:12), K.mu, K.rho, K.oms, K.ms);
xf1dot = PBRFBP_rhs(t2, xf1, K.mu, K.rho, K.oms, K.ms);
xf2dot = PBRFBP_rhs(t3, xf2, K.mu, K.rho, K.oms, K.ms);
xf3dot = PBRFBP_rhs(tN, xf3, K.mu, K.rho, K.oms, K.ms);

% Compute the partial derivatives of the state matching constraints with respect to time
dz1dt1  = -PHI12 * x1dot + xf1dot * 2 / 3;    % Derivative of z1 with respect to t1
dz1dtN  = xf1dot / 3;                        % Derivative of z1 with respect to tN
dz2dt1  = -PHI23 * x2dot * 2 / 3 + xf2dot / 3; % Derivative of z2 with respect to t1
dz2dtN  = -PHI23 * x2dot / 3 + xf2dot * 2 / 3; % Derivative of z2 with respect to tN
dz3dt1  = -PHI34 * x3dot / 3;                % Derivative of z3 with respect to t1
dz3dtN  = -PHI34 * x3dot * 2 / 3 + xf3dot;   % Derivative of z3 with respect to tN

% Compute the gradients of the initial and final state variables (positions and velocities)
dpsiidx1 = [2 * (x1 + K.mu), 2 * y1, 0, 0; vx1, vy1, x1 + K.mu, y1];    % Gradient for initial state
dpsifdxN = [2 * (xN + K.mu - 1), 2 * yN, 0, 0; vxN, vyN, xN + K.mu - 1, yN]; % Gradient for final state

% Initialize the gradient matrix Gc as a zero matrix of size 16x18
Gc = zeros(16, 18);

% Assign the computed gradients to the appropriate blocks in Gc
Gc(1:4, 1:4)    = PHI12;          % State transition matrix for the first segment
Gc(1:4, 5:8)    = -eye(4);        % Gradient of state matching at the first time point
Gc(5:8, 5:8)    = PHI23;          % State transition matrix for the second segment
Gc(5:8, 9:12)   = -eye(4);        % Gradient of state matching at the second time point
Gc(9:12, 9:12)  = PHI34;          % State transition matrix for the third segment
Gc(9:12, 13:16) = -eye(4);        % Gradient of state matching at the third time point

% Assign the computed time derivatives of the state matching constraints to the matrix Gc
Gc(1:4, 17)     = dz1dt1;  % Time derivative for first constraint
Gc(5:8, 17)     = dz2dt1;  % Time derivative for second constraint
Gc(9:12, 17)    = dz3dt1;  % Time derivative for third constraint

Gc(1:4, 18)     = dz1dtN;  % Time derivative for first constraint at final time
Gc(5:8, 18)     = dz2dtN;  % Time derivative for second constraint at final time
Gc(9:12, 18)    = dz3dtN;  % Time derivative for third constraint at final time

% Assign the gradients of the initial and final state variables (Earth and Moon proximity) to Gc
Gc(13:14, 1:4)  = dpsiidx1;   % Gradient for proximity to Earth at initial state
Gc(15:16, 13:16) = dpsifdxN;  % Gradient for proximity to Moon at final state

% Return the transposed gradient matrix
Gc = Gc.';  % Transpose to match the expected dimensions (18x16 matrix)

end

function Gg = gradient_g_multiple(Y)
% Computes the gradient of the inequality constraint functions with respect to the state and time variables.
% The function calculates the gradients of the distance constraints (proximity to the Earth and Moon).
%
% Inputs:
%   Y - State vector containing positions, velocities, and times for multiple shooting points.
%
% Outputs:
%   Gg - The gradient of the inequality constraints with respect to the state and time variables.

% Retrieve system constants
K = constant();

% Extract positions and velocities for the four shooting points
x1 = Y(1); y1 = Y(2);
x2 = Y(5); y2 = Y(6);
x3 = Y(9); y3 = Y(10);
x4 = Y(13); y4 = Y(14);

% Compute the gradients of the distance constraints for each shooting point
S1 = [-2 * (x1 + K.mu), -2 * y1, 0, 0; -2 * (x1 + K.mu - 1), -2 * y1, 0, 0];
S2 = [-2 * (x2 + K.mu), -2 * y2, 0, 0; -2 * (x2 + K.mu - 1), -2 * y2, 0, 0];
S3 = [-2 * (x3 + K.mu), -2 * y3, 0, 0; -2 * (x3 + K.mu - 1), -2 * y3, 0, 0];
S4 = [-2 * (x4 + K.mu), -2 * y4, 0, 0; -2 * (x4 + K.mu - 1), -2 * y4, 0, 0];

% Initialize the gradient matrix Gg as a zero matrix
Gg = zeros(8, 18);

% Assign the computed gradients to the corresponding blocks in Gg
Gg(1:2, 1:4)    = S1;
Gg(3:4, 5:8)    = S2;
Gg(5:6, 9:12)   = S3;
Gg(7:8, 13:16)  = S4;

% Return the transposed gradient matrix
Gg = Gg.';  % Transpose to match the expected dimensions

end

function eta = eta(xx)
% Computes the inequality constraints for proximity to the Earth and Moon based on the state.
%
% Inputs:
%   xx - State vector containing position and velocity.
%
% Outputs:
%   eta - The inequality constraints (proximity to the Earth and Moon).

% Retrieve system constants
K = constant();

% Initialize eta as NaN
eta = nan(2, 1);

% Proximity to Earth constraint
eta(1) = (K.Re / K.DU)^2 - (xx(1) + K.mu)^2 - xx(2)^2;

% Proximity to Moon constraint
eta(2) = (K.Rm / K.DU)^2 - (xx(1) + K.mu - 1)^2 - xx(2)^2;

end

function theta = initial_epoch(et)
% Computes the angle theta representing the Sun's position in the orbital plane
% of the Moon at a specific epoch (et), in the ECLIPJ2000 reference frame.
%
% Inputs:
%   et - Ephemeris time (epoch) for which to compute the angle theta.
%
% Outputs:
%   theta - Angle (in radians) of the Sun's position relative to the Moon's orbital plane.

% Define the center of the coordinate system (Earth-Moon-Barycenter) and the reference frame (ECLIPJ2000).
center = 'EMB';  % Earth-Moon-Barycenter
frame = 'ECLIPJ2000';  % Ecliptic J2000 frame

% Get the state vectors (position and velocity) of the Moon and the Sun at the given epoch (et).
% The function cspice_spkezr returns the state vector of the target body (Moon or Sun) relative to the center (EMB)
% in the specified reference frame (ECLIPJ2000).
rv_moon = cspice_spkezr('Moon', et, frame, 'NONE', center);  % Moon's state vector [position; velocity]
r_sun = cspice_spkpos('Sun', et, frame, 'NONE', center);    % Sun's position vector 

% Normalize the Moon's position and velocity vectors to unit vectors.
% This is done to simplify the calculations and remove the effect of the distance scale.
rm = rv_moon(1:3)/norm(rv_moon(1:3));  % Unit vector for Moon's position
vm = rv_moon(4:6)/norm(rv_moon(4:6));  % Unit vector for Moon's velocity

% Compute the three orthonormal vectors for the orbital plane of the Moon:
% - u1: Unit vector in the direction of the Moon's position vector (rm)
% - u2: Unit vector in the plane perpendicular to both rm and vm (cross product)
% - u3: Unit vector perpendicular to both rm and vm (cross product)
u1 = rm;  % First unit vector (aligned with Moon's position)
u3 = cross(rm, vm)/norm(cross(rm, vm));  % Perpendicular to rm and vm (orbital normal)
u2 = cross(u3, u1)/norm(cross(u3, u1));  % Perpendicular to both u1 and u3 (in orbital plane)

% Combine the three unit vectors into a 3x3 rotation matrix (U)
U = [u1, u2, u3].';  % Rotation matrix from the ECLIPJ2000 frame to the Moon's orbital plane

% Rotate the Sun's position vector (r_sun) into the Moon's orbital plane using the matrix U
% This gives the Sun's position in the orbital plane of the Moon.
r_sun_rot = U * r_sun;  % Rotate the Sun's position into the orbital plane

% Compute the angle theta between the Sun's position and the x-axis in the orbital plane.
% The angle is calculated using the `atan2` function, which returns the angle between the positive x-axis
% and the vector (r_sun_rot(1), r_sun_rot(2)) in the plane, handling full 360-degree angles.
theta = wrapTo2Pi(atan2(r_sun_rot(2), r_sun_rot(1)));  % Angle of the Sun's position in the orbital plane

end

function [bodies] = nbody_init(labels)
%NBODY_INIT Initialize planetary data for n-body propagation
%   Given a set of labels of planets and/or barycentres, returns a
%   cell array populated with structures containing the body label and the
%   associated gravitational constant.
%
%
% Author
%   Name: ALESSANDRO 
%   Surname: MORSELLI
%   Research group: DART
%   Department: DAER
%   University: Politecnico di Milano 
%   Creation: 26/09/2021
%   Contact: alessandro.morselli@polimi.it
%   Copyright: (c) 2021 A. Morselli, Politecnico di Milano. 
%                  All rights reserved.
%
%
% Notes:
%   This material was prepared to support the course 'Satellite Guidance
%   and Navigation', AY 2021/2022.
%
%
% Inputs:
%   labels : [1,n] cell-array with object labels
%
% Outputs:
%   bodies : [1,n] cell-array with struct elements containing the following
%                  fields
%                  |
%                  |--bodies{i}.name -> body label
%                  |--bodies{i}.GM   -> gravitational constant [km**3/s**2]
%
%
% Prerequisites:
%   - MICE (Matlab SPICE)
%   - Populated kernel pool (PCK kernels)
%

% Initialize output
bodies = cell(size(labels));

% Loop over labels
for i = 1:length(labels)
    % Store body label
    bodies{i}.name = labels{i};
    % Store body gravitational constant
    bodies{i}.GM   = cspice_bodvrd(labels{i}, 'GM', 1);
end

end

function [dxdt] = nbody_shift_rhs(t, x, bodies, frame, center)
% Computes the right-hand side of the differential equations for the N-body problem.
% This function calculates the time derivatives of the state vector (position and velocity)
% for a body under the gravitational influence of multiple other bodies in the system.
%
% Inputs:
%   t        - Current time at which to evaluate the right-hand side.
%   x        - State vector of the body being modeled [x; y; z; vx; vy; vz].
%   bodies   - Cell array containing structs for each body, with fields 'name' and 'GM' (gravitational parameter).
%   frame    - Reference frame in which the calculation is performed (e.g., 'ECLIPJ2000').
%   center   - Name of the body being modeled (the center of the system).
%
% Outputs:
%   dxdt     - Time derivatives of the state vector [dx/dt; dy/dt; dz/dt; d(vx)/dt; d(vy)/dt; d(vz)/dt],
%              representing the velocity and acceleration of the body.
%
% Calls:
%   cspice_spkezr  - CSPICE function to get the state vector of a body in a specified frame.
%
% The function calculates the velocity (dx/dt = velocity) and acceleration (d^2x/dt^2 = gravitational force)
% acting on the body from the gravitational pull of all other bodies in the system.

% Initialize the derivative vector dxdt (6x1) for the state vector [position; velocity]
dxdt = zeros(6,1);  

% The first three elements of dxdt correspond to the velocity components (dx/dt = velocity).
dxdt(1:3) = x(4:6);  % Position -> velocity

% Extract the position vector from the state vector
rr = x(1:3);  % Position of the body in the frame

% Initialize an index array for finding the body corresponding to the 'center'
indeces = nan(size(bodies));

% Find the index of the center body in the 'bodies' array
for k = 1:length(bodies)
    indeces(k) = strcmp(bodies{k}.name, center);  % Compare names of bodies
end
index = find(indeces);  % Index of the 'center' body
GM0 = bodies{index}.GM;  % Gravitational parameter of the 'center' body

% Get the position of the center body in the given reference frame (from the Solar System Barycenter 'SSB')
rv_ssb_center = cspice_spkezr(center, t, frame, 'NONE', 'SSB');  % State vector of the center body

% Calculate the magnitude of the position vector rr of the body in question
r = norm(rr);  

% Compute the acceleration of the body due to the gravitational attraction of the center body
aa = -GM0 * rr / r^3;  % Acceleration from the 'center' body

% Loop through all bodies in the system to compute the gravitational influence from each body
for i = 1 : length(bodies)
    if i ~= index  % Skip the center body itself (it has already been accounted for)
        GMi = bodies{i}.GM;  % Gravitational parameter of the i-th body

        % Get the state vector of the i-th body in the specified reference frame
        rv_ssb_body = cspice_spkezr(bodies{i}.name, t, frame, 'NONE', 'SSB');  % State vector of body i

        % Compute the vector from the center body to the i-th body
        ll = rv_ssb_body(1:3) - rv_ssb_center(1:3);  % Relative position of i-th body from the center body
        l = norm(ll);  % Distance between the center body and the i-th body

        % Compute the relative position vector between the body and the i-th body
        dd = rr - ll;  
        d = norm(dd);  % Distance between the body and the i-th body

        % Compute the gravitational acceleration contribution from the i-th body
        aa_i = - GMi * (dd / d^3 + ll / l^3);  % Gravitational acceleration from the i-th body

        % Add this contribution to the total acceleration
        aa = aa + aa_i;
    end
end

% The last three elements of dxdt correspond to the acceleration (d^2x/dt^2)
dxdt(4:6) = aa;  % Acceleration due to the gravitational forces from all bodies

end

function plotStyle
% plotStyle - Applies a custom style to plots in MATLAB.
% This function sets default properties for plot elements such as 
% text, axes, lines, and legends to ensure consistent and aesthetically 
% pleasing formatting throughout your plots.
%
% It applies settings for:
%   - LaTeX interpretation for text, axes, and legends
%   - Grid lines for axes
%   - Line width, marker size, and marker colors
%   - Legend appearance
%   - Font size for axes labels and legends

% Set default interpreter to LaTeX for text, axes, and legends
set(0, 'defaultTextInterpreter', 'Latex');  % Ensure all text is interpreted as LaTeX
set(0, 'defaultAxesTickLabelInterpreter', 'Latex');  % Set LaTeX for tick labels
set(0, 'defaultLegendInterpreter', 'Latex');  % Set LaTeX for legend text

% Enable grid lines on both X and Y axes
set(0, 'defaultAxesXGrid', 'on');  % Enable X-axis grid
set(0, 'defaultAxesYGrid', 'on');  % Enable Y-axis grid

% Set default line properties for plots
set(0, 'defaultLineLineWidth', 1.5);  % Set line width to 1.5 for all lines
set(0, 'defaultLineMarkerSize', 6);   % Set default marker size to 6
set(0, 'defaultLineMarkerEdgeColor', 'k');  % Set marker edge color to black
set(0, 'defaultLineMarkerFaceColor', 'auto');  % Set marker face color to automatically match the plot color

% Set default properties for the legend
set(0, 'defaultLegendLocation', 'northoutside');  % Place legend outside the plot, above
set(0, 'defaultLegendOrientation', 'horizontal');  % Display legend items horizontally
set(0, 'defaultLegendFontSize', 12);  % Set font size of the legend to 12

% Set default font size for axes labels and tick labels
set(0, 'defaultAxesFontSize', 16);  % Set the font size for axes labels and tick labels to 16

end

function plot_circle(r, x_center, y_center, DispName, Color)

   
    % Define the circle
    theta = linspace(0, 2*pi, 500).';
    x = x_center + r*cos(theta);
    y = y_center + r*sin(theta);

    % Plot the circle
    plot(x, y, 'Color', Color, 'DisplayName', DispName);

end

