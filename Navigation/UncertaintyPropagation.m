clear; close all; clc;
plotStyle;

%% Point 1: LinCov and UT Methods
% Initial position (rri) and velocity (vvi) vectors
rri = [-0.011965533749906; -0.017025663128129];
vvi = [10.718855256727338; 0.116502348513671];

% Initial and final times
ti  = 1.282800225339865;
tf  = 9.595124551366348;

% Initial covariance matrix
P0  = [1.041*1e-15, 6.026*1e-17, 5.647*1e-16, 4.577*1e-15;
       6.026*1e-17, 4.287*1e-18, 4.312*1e-17, 1.855*1e-16; 
       5.647*1e-16, 4.312*1e-17, 4.432*1e-16, 1.455*1e-15;
       4.577*1e-15, 1.855*1e-16, 1.455*1e-15, 2.822*1e-14];

% Combine initial position and velocity into the initial state vector
x0 = [rri; vvi];

% Define a time grid with 5 equally spaced points between ti and tf
t_grid = linspace(ti, tf, 5);

% Compute the propagated states and covariances using LinCov method
[YY_lin, PP_lin] = LinCov_method(x0, P0, t_grid);

% Tuning parameters for Unscented Transform (UT)
alpha = 1; beta = 2;

% Compute the propagated states and covariances using UT method
[YY_ut, PP_ut] = UT_method(x0, P0, alpha, beta, t_grid);

% LinCov results
mur_lin = YY_lin(1:2, :);    % Position mean
muv_lin = YY_lin(3:4, :);    % Velocity mean
Pr_lin = PP_lin(1:2, 1:2, :); % Position covariance
Pv_lin = PP_lin(3:4, 3:4, :); % Velocity covariance

% UT results
mur_ut = YY_ut(1:2, :);      % Position mean
muv_ut = YY_ut(3:4, :);      % Velocity mean
Pr_ut = PP_ut(1:2, 1:2, :);  % Position covariance
Pv_ut = PP_ut(3:4, 3:4, :);  % Velocity covariance

% Mean and covariance at the final time step
mu_lin = mur_lin(:, end);    % LinCov final mean
mu_ut  = mur_ut(:, end);     % UT final mean
P_lin  = Pr_lin(:, :, end);  % LinCov final covariance
P_ut   = Pr_ut(:, :, end);   % UT final covariance

% Compute principal standard deviations (square root of eigenvalues)
sigma_princ_lin = sqrt(eig(P_lin));
disp('sigma_princ_lincov=')
disp(sigma_princ_lin);

sigma_princ_ut = sqrt(eig(P_ut));
disp('sigma_princ_ut=')
disp(sigma_princ_ut);

figure; 
hold on;
axis equal;

% Plot the confidence ellipses at 3σ level
plot_conf_ellipse(P_lin, mu_lin, 3, [0, 0.4470, 0.7410], 'Ellipse \ LinCov');  % LinCov
plot_conf_ellipse(P_ut, mu_ut, 3, [0.8500, 0.3250, 0.0980], 'Ellipse \ UT');    % UT

% Plot the means as points
plot(mu_lin(1), mu_lin(2), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0, 0.4470, 0.7410], 'MarkerFaceColor', [0, 0.4470, 0.7410], 'DisplayName', '$\hat{\mathbf{r}} \quad LinCov$')
plot(mu_ut(1), mu_ut(2), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.8500, 0.3250, 0.0980], 'MarkerFaceColor', [0.8500, 0.3250, 0.0980], 'DisplayName', '$\hat{\mathbf{r}} \quad UT$')

% Add legend and labels
legend('Location', 'northeast', 'Orientation', 'vertical');
xlabel('x [-]');
ylabel('y [-]');

% Zoomed inset plot
axes('position', [.2 .2 .19 .19]);
box on; 
hold on;

% Plot the zoomed-in ellipses and means
plot(mu_lin(1), mu_lin(2), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0, 0.4470, 0.7410], 'MarkerFaceColor', [0, 0.4470, 0.7410], 'DisplayName', '$\hat{\mathbf{r}} \quad LinCov$')
plot(mu_ut(1), mu_ut(2), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.8500, 0.3250, 0.0980], 'MarkerFaceColor', [0.8500, 0.3250, 0.0980], 'DisplayName', '$\hat{\mathbf{r}} \quad UT$')
plot_conf_ellipse(P_lin, mu_lin, 3, [0, 0.4470, 0.7410], 'Ellipse \ LinCov');  % LinCov
plot_conf_ellipse(P_ut, mu_ut, 3, [0.8500, 0.3250, 0.0980], 'Ellipse \ UT');    % UT

% Set zoomed-in plot limits
xlim([0.98685, 0.98695]);
ylim([-4.75, -4.6]*1e-3);
axis equal;

%% Point 2: Monte Carlo simulation and Comparison with other methods
n_sam = 1000; % Number of Monte Carlo samples
[YY_mc, PP_mc, Y_tens] = MC_method(x0, P0, t_grid, n_sam);

% Extract position and velocity means from MC simulation
mur_mc = YY_mc(1:2, :);  % Mean position over time
muv_mc = YY_mc(3:4, :);  % Mean velocity over time
Pr_mc = PP_mc(1:2, 1:2, :);  % Position covariance over time
Pv_mc = PP_mc(3:4, 3:4, :);  % Velocity covariance over time

% Final mean and covariance of position
mu_mc = mur_mc(:, end);        % Final mean position
P_mc = Pr_mc(:, :, end);       % Final covariance
sigma_princ_mc = sqrt(eig(P_mc)); % Principal standard deviations
disp('sigma_princ_mc=')
disp(sigma_princ_mc);          % Display standard deviations

% Extract final sampled positions
r_sam = Y_tens(1:2, :, end);  % Final sampled positions

figure;
hold on;

% Plot the confidence ellipses at 3σ level
plot_conf_ellipse(P_lin, mu_lin, 3, [0, 0.4470, 0.7410], 'Ellipse \ LinCov');  % LinCov
plot_conf_ellipse(P_ut, mu_ut, 3, [0.8500, 0.3250, 0.0980], 'Ellipse \ UT');    % UT
plot_conf_ellipse(P_mc, mu_mc, 3, [0.4660, 0.6740, 0.1880], 'Ellipse \ MC');    % MC

% Plot the means as points
plot(mu_lin(1), mu_lin(2), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0, 0.4470, 0.7410], ...
    'MarkerFaceColor', [0, 0.4470, 0.7410], 'DisplayName', '$\hat{\mathbf{r}} \quad LinCov$');
plot(mu_ut(1), mu_ut(2), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor', [0.8500, 0.3250, 0.0980], 'DisplayName', '$\hat{\mathbf{r}} \quad UT$');
plot(mu_mc(1), mu_mc(2), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880], ...
    'MarkerFaceColor', [0.4660, 0.6740, 0.1880], 'DisplayName', '$\hat{\mathbf{r}} \quad MC$');

% Plot sampled positions from Monte Carlo
plot(r_sam(1, :), r_sam(2, :), '.', 'Color', 'k', 'DisplayName', 'MC \ samples');

% Legend and labels
legend('Location', 'northeast', 'Orientation', 'vertical');
xlabel('x [-]');
ylabel('y [-]');
axis equal;
%%
m = length(t_grid);  % Number of time steps
lr_lin = zeros(m, 1); lv_lin = lr_lin;  % LinCov uncertainties
lr_ut = lr_lin; lv_ut = lr_lin;         % UT uncertainties
lr_mc = lr_lin; lv_mc = lr_lin;         % MC uncertainties

for k = 1:m
    lr_lin(k) = 3 * sqrt(max(eig(Pr_lin(:, :, k)))); % LinCov position
    lv_lin(k) = 3 * sqrt(max(eig(Pv_lin(:, :, k)))); % LinCov velocity
    
    lr_ut(k) = 3 * sqrt(max(eig(Pr_ut(:, :, k))));   % UT position
    lv_ut(k) = 3 * sqrt(max(eig(Pv_ut(:, :, k))));   % UT velocity
    
    lr_mc(k) = 3 * sqrt(max(eig(Pr_mc(:, :, k))));   % MC position
    lv_mc(k) = 3 * sqrt(max(eig(Pv_mc(:, :, k))));   % MC velocity
end

% Plot Time History of Position Uncertainty
figure;
subplot(3, 1, 1);
semilogy(t_grid, lr_lin, 'o--', 'Color', [0, 0.4470, 0.7410]);
xlabel('Time [-]');
ylabel('$3\sqrt{\max(\lambda_{i}(P_r))}$', 'Interpreter', 'latex');
title('LinCov');
grid on;

subplot(3, 1, 2);
semilogy(t_grid, lr_ut, 'o--', 'Color', [0.8500, 0.3250, 0.0980]);
xlabel('Time [-]');
ylabel('$3\sqrt{\max(\lambda_{i}(P_r))}$', 'Interpreter', 'latex');
title('UT');
grid on;

subplot(3, 1, 3);
semilogy(t_grid, lr_mc, 'o--', 'Color', [0.4660, 0.6740, 0.1880]);
xlabel('Time [-]');
ylabel('$3\sqrt{\max(\lambda_{i}(P_r))}$', 'Interpreter', 'latex');
title('MC');
grid on;

% Plot Time History of Velocity Uncertainty
figure;
subplot(3, 1, 1)
semilogy(t_grid, lv_lin, 'o--', 'Color', [0, 0.4470, 0.7410]);
xlabel('Time [-]');
ylabel('$3\sqrt{\max(\lambda_{i}(P_v))}$', 'Interpreter', 'latex');
title('LinCov');
grid on;

subplot(3, 1, 2)
semilogy(t_grid, lv_ut, 'o--', 'Color', [0.8500, 0.3250, 0.0980]);
xlabel('Time [-]');
ylabel('$3\sqrt{\max(\lambda_{i}(P_v))}$', 'Interpreter', 'latex');
title('UT');
grid on;

subplot(3, 1, 3)
semilogy(t_grid, lv_mc, 'o--', 'Color', [0, 0.4470, 0.7410]);
xlabel('Time [-]');
ylabel('$3\sqrt{\max(\lambda_{i}(P_v))}$', 'Interpreter', 'latex');
title('MC');
grid on;

%%
% Extract the position and velocity components at the final time step
x = Y_tens(1, :, end);   % Position x at the final time step
y = Y_tens(2, :, end);   % Position y at the final time step
vx = Y_tens(3, :, end);  % Velocity vx at the final time step
vy = Y_tens(4, :, end);  % Velocity vy at the final time step

% Standardize the extracted data
zx = standardize(x);   % Standardized position x
zy = standardize(y);   % Standardized position y
zvx = standardize(vx); % Standardized velocity vx
zvy = standardize(vy); % Standardized velocity vy
%%
% Create a 2x2 grid of QQ plots
figure;

% Plot the QQ plot for standardized position x
subplot(2, 2, 1)
qqplot(zx); % Use standardized data for x
ylabel('Quantiles of $\overline{x}$', 'Interpreter','latex')
title('')

% Plot the QQ plot for standardized position y
subplot(2, 2, 2)
qqplot(zy)  % Use standardized data for y
ylabel('Quantiles of $\overline{y}$', 'Interpreter','latex')
title('')

% Plot the QQ plot for standardized velocity vx
subplot(2, 2, 3)
qqplot(zvx)  % Use standardized data for vx
ylabel('Quantiles of $\overline{v_x}$', 'Interpreter','latex')
title('')

% Plot the QQ plot for standardized velocity vy
subplot(2, 2, 4)
qqplot(zvy)  % Use standardized data for vy
ylabel('Quantiles of $\overline{v_y}$', 'Interpreter','latex')
title('')

nbins = round(sqrt(n_sam));
figure
subplot(2, 2, 1)
hold on
histogram(zx, nbins)
xlabel('$\overline{x}$ [-]')
ylabel('N. samples')

subplot(2, 2, 2)
hold on
histogram(zy, nbins)
xlabel('$\overline{y}$ [-]')
ylabel('N. samples')

subplot(2, 2, 3)
hold on
histogram(zvx, nbins)
xlabel('$\overline{v}_x$ [-]')
ylabel('N. samples')

subplot(2, 2, 4)
hold on
histogram(zvy, nbins)
xlabel('$\overline{v}_y$ [-]')
ylabel('N. samples')
%%

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
constants.Re  = 6.378136600000000e+03;  % Earth's radius (km)
constants.Rm  = 1.737400000000000e+03;  % Moon's radius (km)
constants.hi  = 167;                    % Initial altitude above Earth's surface (km)
constants.hf  = 100;                    % Final altitude above Moon's surface (km)
constants.DU  = 3.84405000e5;           % Characteristic distance unit: Earth-Moon distance (km)
constants.TU  = 4.34256461;             % Characteristic time unit (Earth-Moon system)
constants.VU  = 1.02454018;             % Characteristic velocity unit (normalized units)

% Compute normalized radial distances
constants.ri  = (constants.Re + constants.hi) / constants.DU;  % Initial distance (normalized)
constants.rf  = (constants.Rm + constants.hf) / constants.DU;  % Final distance (normalized)

% Gravitational and orbital parameters
constants.mu  = 0.012150584269940;  % Gravitational parameter ratio (Moon/Earth+Moon)
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

function plot_conf_ellipse(P, mu, n_sigma, Color, DisplayName)
    % Plot a confidence ellipse for a 2D Gaussian distribution
    % P: Covariance matrix (2x2)
    % mu: Mean vector (2x1)
    % n_sigma: Confidence level (number of standard deviations)

    % Check that P is a 2x2 matrix and mu is a 2x1 vector
    if size(P, 1) ~= 2 || size(P, 2) ~= 2
        error('Covariance matrix P must be 2x2');
    end
    if length(mu) ~= 2
        error('Mean vector mu must be 2x1');
    end
    
    % Create a unit circle parametrized by angles alpha
    alpha = linspace(0, 2*pi, 100);
    circle = [cos(alpha); sin(alpha)];

    % Perform Singular Value Decomposition (SVD) of the covariance matrix P
    [R, D] = svd(P);  % R is the rotation matrix, D is the diagonal matrix with eigenvalues

    % Extract the square root of eigenvalues from D (standard deviations along each axis)
    d = sqrt(diag(D));  % d will be a 2x1 vector of standard deviations

    % Apply the transformation: R * d scales the unit circle and rotates it
    ellipse = n_sigma * R * diag(d) * circle;

    % Compute the x and y coordinates of the ellipse
    x = mu(1) + ellipse(1, :);  % x-coordinates of the ellipse
    y = mu(2) + ellipse(2, :);  % y-coordinates of the ellipse

    % Plot the ellipse
    plot(x, y, 'Color', Color, 'DisplayName', DisplayName);
end

function [yhat_tens, Py_tens] = LinCov_method(x0, P0, t_grid)
% Linearized covariance propagation method
%
% Inputs:
%   x0     - Initial state vector
%   P0     - Initial covariance matrix
%   t_grid - Time grid for propagation
%
% Outputs:
%   yhat_tens - Mean state at each time step (n x m)
%   Py_tens   - Covariance matrices at each time step (n x n x m)

    % Constants
    K = constant();
    n = length(x0);  % State dimension
    m = length(t_grid); % Number of time steps

    % Preallocate results
    yhat_tens = zeros(n, m); % Mean state
    Py_tens = zeros(n, n, m); % Covariance matrices

    % Initial state and STM
    PHI0 = eye(n); % Initial state transition matrix
    X0 = [x0; PHI0(:)]; % Combine state and flattened STM

    % ODE solver options
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

    % Propagate state and STM using ode78
    [~, X] = ode78(@(t, X) PBRFBP_STM(t, X, K.mu, K.rho, K.oms, K.ms), t_grid, X0, options);

    % Process results
    for k = 1:m
        % Extract state at time step k
        yhat_tens(:, k) = X(k, 1:n).'; % Extract mean state

        % Extract and reshape STM
        PHI = reshape(X(k, n+1:end), n, n); % Reshape flattened STM

        % Propagate covariance
        Py_tens(:, :, k) = PHI * P0 * PHI.'; % Covariance propagation
    end
end

function [yhat_tens, Py_tens] = UT_method(x0, P0, alpha, beta, t_grid)
% Unscented Transform method for uncertainty propagation
%
% Inputs:
%   x0     - Initial state vector
%   P0     - Initial covariance matrix
%   alpha  - UT scaling parameter
%   beta   - UT tuning parameter
%   t_grid - Time grid for propagation
%
% Outputs:
%   yhat_tens - Mean state at each time step (n x m)
%   Py_tens   - Covariance matrices at each time step (n x n x m)

    % Problem constants
    n = length(x0);      % State dimension
    m = length(t_grid);  % Number of time steps
    K = constant();      % Constants for the dynamical system

    % UT scaling parameters
    lambda = n * (alpha^2 - 1);
    gamma = sqrt(n + lambda);

    % Compute sigma points and weights
    try
        L = chol(P0, 'lower'); % Cholesky decomposition (lower triangular)
    catch
        error('Initial covariance matrix P0 is not positive definite.');
    end

    B = gamma * L; % Scaled square root of P0
    Wm = [lambda / (n + lambda), repmat(1 / (2 * (n + lambda)), 1, 2 * n)];
    Wc = [lambda / (n + lambda) + (1 - alpha^2 + beta), repmat(1 / (2 * (n + lambda)), 1, 2 * n)];

    % Generate sigma points
    Xi = [x0, x0 + B, x0 - B]; % Shape: n x (2n+1)

    % Preallocate tensors
    y = zeros(n, m, 2 * n + 1); % Propagated states
    yhat_tens = zeros(n, m);    % Mean state
    Py_tens = zeros(n, n, m);   % Covariance matrices

    % ODE solver options
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

    % Propagate each sigma point
    for k = 1:2 * n + 1
        % Initial state for sigma point
        x0_sigma = Xi(:, k);

        % Propagate using ODE solver
        [~, x] = ode78(@(t, x) PBRFBP_rhs(t, x, K.mu, K.rho, K.oms, K.ms), t_grid, x0_sigma, options);

        % Store propagated state
        y(:, :, k) = x.'; % Transpose to match dimensions
    end

    % Compute mean and covariance at each time step
    for i = 1:m
        % Extract propagated states at time t_grid(i)
        Y_i = squeeze(y(:, i, :)); % Shape: (n x 2n+1)

        % Compute mean
        yhat = Y_i * Wm'; % Weighted sum of sigma points

        % Compute covariance
        Py = zeros(n, n);
        for j = 1:2 * n + 1
            diff = Y_i(:, j) - yhat;
            Py = Py + Wc(j) * (diff * diff');
        end

        % Store results
        yhat_tens(:, i) = yhat;
        Py_tens(:, :, i) = Py;
    end
end

function [yhat_tens, Py_tens, Y_tens] = MC_method(x0, P0, t_grid, n_sam)
% Monte Carlo method for uncertainty propagation
%
% Inputs:
%   x0     - Initial state vector
%   P0     - Initial covariance matrix
%   t_grid - Time grid for propagation
%   n_sam  - Number of Monte Carlo samples
%
% Outputs:
%   yhat_tens - Mean state at each time step (n x m)
%   Py_tens   - Covariance matrix at each time step (n x n x m)
%   Y_tens    - Propagated samples (n x n_sam x m)

    % Constants
    K = constant();
    n = length(x0);  % State dimension
    m = length(t_grid);  % Number of time steps
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

    % Preallocate tensors
    Y_tens = nan(n, n_sam, m);
    yhat_tens = zeros(n, m);
    Py_tens = zeros(n, n, m);

    % Monte Carlo sampling and propagation
    for k = 1:n_sam
        % Sample initial state from the distribution
        x0_sample = mvnrnd(x0, P0)';
        
        % Propagate using ode78
        [~, x] = ode78(@(t, x) PBRFBP_rhs(t, x, K.mu, K.rho, K.oms, K.ms), t_grid, x0_sample, options);

        % Store results in Y_tens
        Y_tens(:, k, :) = permute(x, [2, 1]); % Transpose and permute to match dimensions
    end

    % Compute mean and covariance at each time step
    for i = 1:m
        % Extract samples at the i-th time step
        Y_i = squeeze(Y_tens(:, :, i)); % Shape: (n x n_sam)

        % Compute sample mean
        yhat = mean(Y_i, 2);

        % Compute sample covariance
        Py = cov(Y_i');

        % Store results
        yhat_tens(:, i) = yhat;
        Py_tens(:, :, i) = Py;
    end
end

function z = standardize(x)
% STANDARDIZE Standardizes a vector to have zero mean and unit variance.
%
% Inputs:
%   x - A vector of samples (1 x n or n x 1).
%
% Outputs:
%   z - Standardized vector (same size as x),
%       where z = (x - mean(x)) / std(x).

    mu = mean(x);       % Compute mean of the vector
    sigma = std(x);     % Compute standard deviation of the vector

    % Check for zero standard deviation to avoid division by zero
    if sigma == 0
        error('Standard deviation is zero. Cannot standardize the vector.');
    end

    z = (x - mu) / sigma; % Standardize the vector
end
