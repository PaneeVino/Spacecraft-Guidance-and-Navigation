clear; close all; clc; 
format long g
plotStyle;
%% Ex1.1: Lagrange Points
% Define initial guesses for Lagrange points
x0_mat = [ 0,  0;    % Initial guess for L1
           1.5,  0;    % Initial guess for L2
          -1.5,  0;    % Initial guess for L3
           0,  1.5;    % Initial guess for L4
           0, -1.5];   % Initial guess for L5

lag_point = nan(size(x0_mat)); % Initialize solution storage matrix

% Iterate over each initial guess and use fsolve to find Lagrange points
options = optimoptions("fsolve", "Display","none");  % Suppress display
options.OptimalityTolerance = 1e-10; % High precision for solution
for i = 1:size(x0_mat, 1)
    % Call fsolve for each initial guess
    solution = fsolve(@gradient_of_potential, x0_mat(i, :), options);
    lag_point(i, :) = solution; % Store the solution
end

% Define gravitational parameter
mu = 0.012150; 

% Initialize storage for Lagrange points and Jacobi constants
lagrange_points = zeros(5, 6); 
C_vec = nan(5, 1);

% Compute Jacobi constant for each Lagrange point
for k = 1 : length(C_vec)
    lagrange_points(k, 1:2) = lag_point(k, :); % Store x and y coordinates
    xx = lagrange_points(k, :); % Extract Lagrange point position
    [C, ~] = jacobi_constant(xx, mu); % Compute Jacobi constant
    C_vec(k) = C; % Store Jacobi constant
end

% Extract x and y coordinates of Lagrange points
x_vec = lagrange_points(:, 1); 
y_vec = lagrange_points(:, 2);

% Names for Lagrange points
names = ["L1"; "L2"; "L3"; "L4"; "L5"]; 

% Create table to display results
lagr_tab = table(names, x_vec, y_vec, C_vec, ...
    'VariableNames', {'Point', 'x', 'y', 'C'});

% Display the table
disp(lagr_tab);

% Plotting the results
x_E = -mu; y_E = 0; % Earth's position
x_M = 1 - mu; y_M = 0; % Moon's position

% Generate grid for contour plot
val = 1.5;
x_vals = linspace(-val, val, 1000); % Range of x-values
y_vals = linspace(-val, val, 1000); % Range of y-values
[X, Y] = meshgrid(x_vals, y_vals); % Create meshgrid

% Evaluate norm of gradient of potential
Z = norm_gradient(X, Y); 

% Plot the contour of the norm of the gradient
cmap = parula; % Colormap
lmin = 0;      % Minimum contour level
lmax = 2;    % Maximum contour level
figure;
hold on;
contour(X, Y, Z, lmin:0.05:lmax); % Contour plot
hcb = colorbar('eastoutside', 'Ticks', linspace(lmin, lmax, 5)); % Add colorbar
hcb.FontSize = 15; % Colorbar font size
colormap(cmap); % Set colormap
title(hcb, 'U [-]', 'interpreter', 'latex', 'fontsize', 20); % Colorbar title
set(hcb, 'TickLabelInterpreter', 'latex'); % Use LaTeX for tick labels
clim([lmin lmax]); % Set colorbar limits
xlabel('x [-]'); % X-axis label
ylabel('y [-]'); % Y-axis label
grid minor; % Add minor grid lines

% Plot Earth's position
plot(x_E, y_E, 'o', "MarkerSize", 20, "MarkerFaceColor", "g");

% Plot Moon's position
plot(x_M, y_M, 'o', "MarkerSize", 12, "MarkerFaceColor", [0.5 0.5 0.5]);

% Plot and annotate Lagrange points
for kk = 1 : length(C_vec)
    xl = lagrange_points(kk, 1); % Lagrange point x-coordinate
    yl = lagrange_points(kk, 2); % Lagrange point y-coordinate
    txt = names(kk); % Lagrange point name
    plot(xl, yl, 'o', "MarkerSize", 8, "MarkerFaceColor", "r"); % Plot point
    text(xl, yl, txt, "FontSize", 15, 'color', 'k', "HorizontalAlignment","left", 'VerticalAlignment','bottom'); % Annotate point
end

% Add legend
legend('', 'Earth', 'Moon', 'Lagrange Points'); % Legend
%% Ex 1.2: Halo orbit with J=3.09
% Define the initial conditions for the halo orbit
x0  = 1.068792441776;
y0  = 0;
z0  = 0.071093328515;
vx0 = 0;
vy0 = 0.319422926485;
vz0 = 0;

% Combine initial conditions into a state vector
xx0  = [x0; y0; z0; vx0; vy0; vz0];
r_M  = [1-mu; 0; 0];
r_l2 = [lag_point(2, 1); lag_point(2, 2); 0];

% Gravitational parameter (mass ratio of the two primary bodies)
mu  = 0.012150;

% Target Jacobi constant
C_target = 3.09;
tf_0 = 3.0;
% Generate initial guess for tf
[~, ~, tf_guess] = propagate_STM(0, xx0, tf_0, mu, true);

% Initialize tf
tf = tf_guess;

% Find the initial conditions for a halo orbit with the target Jacobi constant
[x0, tf, iter] = find_halo(xx0, tf, C_target, mu);
disp('X0=')
disp(x0)
% Set options for the ODE solver with specified tolerances
options = odeset('reltol', 1e-12, 'abstol', 1e-12);

% Integrate the equations of motion from 0 to 2*tf with the found initial conditions
[~, xx] = ode78(@(t,x) CR3BP_rhs(t,x,mu), [0 2*tf], x0, options);

% Transpose the solution matrix for easier indexing
xx = xx.';

% Extract position vectors from the state vector
rr = xx(1:3, :);

% Plot the 3D trajectory of the halo orbit
figure
grid on
hold on
plot3(r_M(1), r_M(2), r_M(3), 'o', "MarkerFaceColor", [0.5 0.5 0.5], ...
    "MarkerSize", 17, "DisplayName", "Moon")
plot3(r_l2(1), r_l2(2), r_l2(3), 'o', "MarkerFaceColor", "r", ...
    "MarkerSize", 8, "DisplayName", "L2 Point")
plot3(rr(1, :), rr(2, :), rr(3, :),'Color',"#0072BD", ...
    "DisplayName", "Halo Orbit")
axis equal
view(15, 15)
xlabel('x [-]')
ylabel('y [-]')
zlabel('z [-]')
xlim([0.95 1.2])
ylim([-0.15 0.15])
zlim([-0.15 0.1]);
%title('Halo Orbit with J = 3.09')
legend

%% Ex 3.3: numerical continuation

% Gravitational parameter (mass ratio of the two primary bodies)
mu  = 0.012150;

% Create a vector for continuation values and corresponding target Jacobi constants
vec = 0.01 * linspace(0, 5, 6);
C_vec = 3.09 - vec;
% Define the custom colors as RGB values
Color = flip([
    0, 0.4470, 0.7410;  % Blue
    0.8500, 0.3250, 0.0980;  % Orange
    0.9290, 0.6940, 0.1250;  % Yellow
    0.4940, 0.1840, 0.5560;  % Purple
    0.4660, 0.6740, 0.1880;  % Green
    0.3010, 0.7450, 0.9330  % Light Blue
]);

% Initialize arrays to store the results of the continuation process
X0   = nan(6, length(C_vec));
Tf   = nan(length(C_vec), 1);
Iter = nan(length(C_vec), 1);

X0(:, 1) = x0;
Tf(1) = tf;

% Perform numerical continuation to find initial conditions for varying Jacobi constants
for i = 2 : length(C_vec)

    % Set the current target Jacobi constant
    C_target = C_vec(i);
    % Find the initial conditions for the current target Jacobi constant
    [x0, tf, iter] = find_halo(xx0, tf,C_target, mu);

    % Store the results
    X0(:, i) = x0;
    Tf(i)    = tf;
    Iter(i)  = iter;
    xx0      = x0;

end

% Set options for the ODE solver with specified tolerances
options = odeset('reltol', 1e-12, 'abstol', 1e-12);

% Plot the results of the numerical continuation
figure
hold on
grid on
plot3(r_M(1), r_M(2), r_M(3), 'o', "MarkerFaceColor", [0.5 0.5 0.5], ...
    "MarkerSize", 17, "DisplayName", "Moon")
plot3(r_l2(1), r_l2(2), r_l2(3), 'o', "MarkerFaceColor", "r", ...
    "MarkerSize", 8, "DisplayName", "L2 Point")
for i =  1 : length(C_vec)

    % Retrieve the initial conditions and final time for the current target Jacobi constant
    x0 = X0(:, i);
    tf = Tf(i);
    % Integrate the equations of motion from 0 to 2*tf with the current initial conditions
    [~, xx] = ode78(@(t,x) CR3BP_rhs(t,x,mu), [0 2*tf], x0, options);
    xx = xx.';  % Transpose the solution matrix for easier indexing
    rr = xx(1:3, :);  % Extract position vectors from the state vector

    % Plot the 3D trajectory of the halo orbit
    plot3(rr(1, :), rr(2, :), rr(3, :), "Color", Color(length(C_vec) - i + 1, :)) 
end
disp('X0=')
disp(X0(:, end))
axis equal
view(15, 15)
xlabel('x [-]')
ylabel('y [-]')
zlabel('z [-]')
xlim([0.95 1.2])
ylim([-0.15 0.15])
zlim([-0.2 0.1]);
legend("Moon", "L2")
colormap(Color);
c = colorbar;  % Add the colorbar
c.Ticks = linspace(min(C_vec), max(C_vec), length(C_vec));  % Set ticks to match C_vec values
%c.TickLabels = C_vec;  % Use rounded values of C_vec as labels
c.FontSize = 15;
clim([min(C_vec)-0.005,max(C_vec)+0.005]);  % Adjust colorbar limits
title(c, 'C [-]', 'interpreter', 'latex', 'fontsize', 20); % Colorbar title
set(c, 'TickLabelInterpreter', 'latex'); % Use LaTeX for tick labels
%c.Position = [0.76, 0.2, 0.03, 0.65];  % Example: Move to a custom position


xlabel('x [-]')
ylabel('y [-]')
zlabel('z [-]')
%title('Numerical Continuation of Halo Orbits')
%%
function U = norm_gradient(x, y)
    % Constants
    mu = 0.012150;
    
    % Distances to the two foci
    r1 = sqrt((x + mu).^2 + y.^2);
    r2 = sqrt((x - (1 - mu)).^2 + y.^2); % Corrected distance to second focus
    
    % Gradients of the potential
    dUdx = x - (1 - mu) ./ r1.^3 .* (x + mu) + mu ./ r2.^3 .* (1-mu-x);
    dUdy = y - (1 - mu) ./ r1.^3 .* y - mu ./ r2.^3 .* y;
    
    % Squared norm of the gradient (squared potential gradient)
    U2 = dUdx.^2 + dUdy.^2;
    U = sqrt(U2);
end

function gradU = gradient_of_potential(xy)
    % Constants
    mu = 0.012150;
    
    % Extract x and y from the input vector
    x = xy(1);
    y = xy(2);
    
    % Distances to the two foci
    r1 = sqrt((x + mu).^2 + y.^2);
    r2 = sqrt((x - (1 - mu)).^2 + y.^2);
    
    % Gradients of the potential
    dUdx = x - (1 - mu) ./ r1.^3 .* (x + mu) + mu ./ r2.^3 .* (1-mu-x);
    dUdy = y - (1 - mu) ./ r1.^3 .* y - mu ./ r2.^3 .* y;
    
    % Return the gradient as a vector (for fsolve to find the zeros)
    gradU = [dUdx; dUdy];
end

function dxdt = CR3BP_rhs(~, xx, mu)
    % CR3BP_rhs: Computes the derivatives for the Circular Restricted Three-Body Problem
    % Inputs:
    %   ~   : Placeholder for time (unused)
    %   xx  : State vector [x, y, z, vx, vy, vz]
    %   mu  : Gravitational parameter (mass ratio of the two primary bodies)
    % Output:
    %   dxdt: Derivative of the state vector

    % Initialize right-hand-side
    dxdt = zeros(6, 1);

    % Extract positions from state vector
    x = xx(1);
    y = xx(2);
    z = xx(3);

    % Extract velocities from state vector
    vx = xx(4);
    vy = xx(5);
    vz = xx(6);

    % Derivative of position is velocity
    dxdt(1:3) = [vx; vy; vz];

    % Compute distances to the two primary bodies
    r1 = sqrt((x + mu)^2 + y^2 + z^2);
    r2 = sqrt((x + mu - 1)^2 + y^2 + z^2);

    % Compute partial derivatives of the effective potential
    dUdx = x - (1 - mu) / r1^3 * (x + mu) + mu / r2^3 * (1 - mu - x);
    dUdy = y - (1 - mu) / r1^3 * y - mu / r2^3 * y;
    dUdz = - (1 - mu) / r1^3 * z - mu / r2^3 * z;

    % Compute accelerations due to gravitational and centrifugal forces
    dxdt(4:6) = [2 * vy + dUdx; -2 * vx + dUdy; dUdz];
end

function PHI = STM_numerical(t0, xx0, tf, xf, mu, evtFlag)
    % STM_numerical: Computes the State Transition Matrix (STM) for the CR3BP
    % Inputs:
    %   t0      : Initial time
    %   xx0     : Initial state vector [x, y, z, vx, vy, vz]
    %   tf      : Final time
    %   xf      : Final state vector [x, y, z, vx, vy, vz]
    %   mu      : Gravitational parameter (mass ratio of the two primary bodies)
    %   evtFlag : Event flag for custom event handling
    % Output:
    %   PHI     : State Transition Matrix (6x6)

    % Initialize PHI with NaN values
    PHI = nan(6, 6);

    % Time of flight
    tof = tf - t0;

    % Set options for the ODE solver with specified tolerances and event function
    options = odeset('reltol', 1e-12, 'abstol', 1e-12, 'Events', @(x,y) xz_plane_crossing(x,y,evtFlag));

    % Loop through each state variable to compute perturbed trajectories
    for i = 1 : 6

        % Initialize perturbation vector
        ee = zeros(6, 1);
        % Determine the perturbation size based on machine precision
        e = max(abs(xx0(i)) * sqrt(eps), sqrt(eps));
        ee(i) = e; 

        % Perturb the initial state
        xx0_pert = xx0 + ee;

        % Integrate the equations of motion from t0 to tf with the perturbed initial state
        [~, xx] = ode78(@(t,x) CR3BP_rhs(t,x,mu), [0 tof], xx0_pert, options);
        xx = xx.';  % Transpose the solution matrix
        xf_pert = xx(1:6, end);  % Extract the perturbed final state

        % Calculate the column of the State Transition Matrix
        PHI(:, i) = (xf_pert - xf) / e;
    end
end

function [STM, xxf, ttf] = propagate_STM(t0, xx0, tf, mu, varargin)
    % propagate_STM: Propagates the State Transition Matrix (STM) and the state vector for the CR3BP
    % Inputs:
    %   t0      : Initial time
    %   xx0     : Initial state vector [x, y, z, vx, vy, vz]
    %   tf      : Final time
    %   mu      : Gravitational parameter (mass ratio of the two primary bodies)
    %   varargin: Optional argument to set event flag
    % Outputs:
    %   STM     : State Transition Matrix (6x6)
    %   xxf     : Final state vector at time tf
    %   ttf     : Final time after integration

    % Check if event flag is provided, otherwise set default to true
    if nargin > 4
        evtFlag = varargin{1};
    else
        evtFlag = true;
    end

    % Calculate the time of flight
    tof = tf - t0;

    % Set options for the ODE solver with specified tolerances and event function
    options = odeset('reltol', 1e-12, 'abstol', 1e-12, 'Events', @(x,y) xz_plane_crossing(x,y,evtFlag));

    % Integrate the equations of motion from t0 to tf
    [tt, xx] = ode78(@(t,x) CR3BP_rhs(t,x,mu), [0 tof], xx0, options);
    xx = xx.';  % Transpose the solution matrix
    xxf = xx(1:6, end);  % Extract the final state vector
    ttf = tt(end);  % Extract the final time

    % Compute the State Transition Matrix using numerical differentiation
    STM = STM_numerical(t0, xx0, tf, xxf, mu, evtFlag);
end

function [C, dJdx] = jacobi_constant(xx, mu)
    % jacobi_constant: Computes the Jacobi constant and its gradient for the CR3BP
    % Inputs:
    %   xx  : State vector [x, y, z, vx, vy, vz]
    %   mu  : Gravitational parameter (mass ratio of the two primary bodies)
    % Outputs:
    %   C   : Jacobi constant
    %   dJdx: Gradient of the Jacobi constant with respect to state vector

    % Extract positions from state vector
    x = xx(1); 
    y = xx(2); 
    z = xx(3);

    % Extract velocities from state vector
    vx = xx(4); 
    vy = xx(5); 
    vz = xx(6);

    % Compute Jacobi constant C
    C = (mu - 1.0) * 1.0 / sqrt((mu + x)^2 + y^2 + z^2) * -2.0 - mu * (mu - 1.0) - vx^2 - vy^2 - vz^2 + x^2 + y^2 + mu * 1.0 / sqrt((mu + x - 1.0)^2 + y^2 + z^2) * 2.0;

    % Compute the gradient of the Jacobi constant (dJdx)
    dJdx = [...
        x * 2.0 + (mu * 2.0 + x * 2.0) * (mu - 1.0) * 1.0 / ((mu + x)^2 + y^2 + z^2)^(3.0 / 2.0) - mu * (mu * 2.0 + x * 2.0 - 2.0) * 1.0 / ((mu + x - 1.0)^2 + y^2 + z^2)^(3.0 / 2.0), ...
        y * 2.0 - mu * y * 1.0 / ((mu + x - 1.0)^2 + y^2 + z^2)^(3.0 / 2.0) * 2.0 + y * (mu - 1.0) * 1.0 / ((mu + x)^2 + y^2 + z^2)^(3.0 / 2.0) * 2.0, ...
        mu * z * 1.0 / ((mu + x - 1.0)^2 + y^2 + z^2)^(3.0 / 2.0) * -2.0 + z * (mu - 1.0) * 1.0 / ((mu + x)^2 + y^2 + z^2)^(3.0 / 2.0) * 2.0, ...
        vx * -2.0, ...
        vy * -2.0, ...
        vz * -2.0 ...
    ];

end

function [xx0, tf, iter] = find_halo(x00, tf, C_target, mu)
    % find_halo: Finds initial conditions for a halo orbit in the CR3BP
    % Inputs:
    %   x00      : Initial state vector guess [x, y, z, vx, vy, vz]
    %   C_target : Target Jacobi constant
    %   tf       : initial final time guess
    %   mu       : Gravitational parameter (mass ratio of the two primary bodies)
    % Outputs:
    %   xx0      : Adjusted initial state vector
    %   tf       : Final time
    %   iter     : Number of iterations performed

    % Extract initial positions and velocities from the guess state vector
    x0 = x00(1);
    y0 = x00(2);
    z0 = x00(3);
    vx0 = x00(4);
    vy0 = x00(5);
    vz0 = x00(6);

    % Initialize errors and tolerance
    err_y  = 1;
    err_vx = 1;
    err_vz = 1;  
    err_C  = 1;

    % Set maximum number of iterations and the desired tolerance
    Nmax    = 50;
    iter    = 0;
    tol     = 1e-13;

    % Initialize perturbation corrections
    dx  = 0;
    dz  = 0;
    dvy = 0;
    dt  = 0;

    % Iterative process to adjust initial conditions
    while (abs(err_y) > tol || abs(err_vx) > tol || abs(err_vz) > tol || abs(err_C) > tol) && (iter < Nmax)

        % Update initial conditions with the previous corrections
        x0 = x0 + dx;
        z0 = z0 + dz;
        vy0 = vy0 + dvy;
        xx0 = [x0; y0; z0; vx0; vy0; vz0];
        tf  = tf + dt;

        % Propagate the state and STM
        [STM, xf, tf] = propagate_STM(0, xx0, tf, mu, false);
        f = CR3BP_rhs(0, xf, mu);
        [C, dJdx] = jacobi_constant(xx0, mu);

        % Construct the matrix A and the vector b for the linear system
        A  =  [STM, f; dJdx, 0];
        A(:, [2, 4, 6]) = [];  % Remove columns for velocities and time
        A([1, 3, 5], :) = [];  % Remove rows for velocities

        % Calculate errors 
        err_y  = 0-xf(2);
        err_vx = 0-xf(4);
        err_vz = 0-xf(6);
        err_C  = C_target-C;
        % Assembly vector
        b = [err_y; err_vx; err_vz; err_C];

        % Solve the linear system for the perturbations
        vet = A \ b;

        % Extract the corrections from the solution vector
        dx  = vet(1);
        dz  = vet(2);
        dvy = vet(3);
        dt  = vet(4);

        % Calculate errors in velocities and Jacobi constant
        err_y   = abs(xf(2));
        err_vx  = abs(xf(4));
        err_vz  = abs(xf(6));
        err_C   = abs(C - C_target);

        % Increment the iteration counter
        iter = iter + 1;
    end
end

function [value, isterminal, direction] = xz_plane_crossing(~, xx, isTerminal)
    % xz_plane_crossing: Event function to detect crossing of the xz-plane
    % Inputs:
    %   ~         : Placeholder for time (unused)
    %   xx        : State vector [x, y, z, vx, vy, vz]
    %   isTerminal: Logical flag to stop integration when the event is detected
    % Outputs:
    %   value     : The value that is checked to determine if an event occurs
    %   isterminal: A flag to determine if the integration should stop when the event occurs
    %   direction : Direction of the zero crossing (0 means all directions)

    % The event occurs when the y-coordinate is zero (crossing the xz-plane)
    value = xx(2) - 0;

    % Terminal event flag (stop integration if true)
    isterminal = isTerminal;

    % Detect crossings in both directions (y going positive or negative)
    direction = 0;
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
