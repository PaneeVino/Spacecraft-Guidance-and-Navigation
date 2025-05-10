clear; cspice_kclear; close all; clc; 
plotStyle; % Apply custom plot styles

% Add paths for required libraries and kernels
% addpath("kernels/")
% addpath("../mice/src/mice/")
% addpath("../mice/lib/")
% addpath("sgp4/")
% addpath("tle/")
cspice_kclear;
cspice_furnsh('assignment02.tm') % Load SPICE kernels
lander = struct('name', 'MOONLANDER', 'El_min', 0, 'lat', 78, 'lon', 15, 'alt', 0);
%% Visibility Windows
rr0 = [4307.844185282820; -1317.980749248651; 2109.210101634011]; % Initial position (km)
vv0 = [-0.110997301537882; -0.509392750828585; 0.815198807994189]; % Initial velocity (km/s)
epoch_0 = '2024-11-18T16:30:00.000';
epoch_f = '2024-11-18T20:30:00.000';

sigma_rho = 100; % Measurement noise (m)

et0 = cspice_str2et(epoch_0); % Convert start epoch to ephemeris time
etf = cspice_str2et(epoch_f); % Convert end epoch to ephemeris time
freq = 30; % Sampling frequency (s)
npoints = round((etf - et0) / freq) + 1; % Number of time points
et_vec = linspace(et0, etf, npoints); % Generate time vector

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); % ODE solver settings
mu = cspice_bodvrd('Moon', 'GM', 1); % Moon's gravitational parameter
odefun = @(t, x) two_body_rhs(t, x, mu); % Define two-body dynamics function
xx0 = [rr0; vv0]; % Initial state vector

% Propagate orbiter state using ODE78
[~, xx_sc_m] = ode78(@(t, x) two_body_rhs(t, x, mu), et_vec, xx0, options);
xx_sc_m = xx_sc_m.'; % Transpose state matrix for consistency

% Transform orbiter position to Moon-centered, Moon-fixed (MCMF) frame
r_sc_mci = xx_sc_m(1:3, :); % Extract positions from state vector
r_sc_mcmf = nan(size(r_sc_mci)); % Initialize Moon-fixed positions
R_IN2MF = cspice_pxform('J2000', 'IAU_MOON', et_vec); % Transformation matrix

for k = 1:length(et_vec)
    r_sc_mcmf(:, k) = R_IN2MF(:, :, k) * r_sc_mci(:, k); % Apply transformation
end

% Compute lander visibility
Y_pred = pointing_lander(lander, et_vec, r_sc_mci, true);
[Y_pred, et_vec_vis] = visibility(Y_pred, et_vec, lander);

% Check if lander is always visible
flag = length(et_vec) == length(et_vec_vis);
if flag
    disp('Lander always visible')
else
    disp('Lander not always visible')
end


t  = datetime(cspice_timout(et_vec,'YYYY-MM-DD HR:MN:SC.###'),...
    'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
figure

plot(t,Y_pred(2,:), 'DisplayName', 'Elevation')
yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName','Minimum Elevation')
ylabel('El [deg]')
ylim([-10 90])
legend

[X, Y, Z] = sphere(30);
radii = cspice_bodvrd('Moon', 'RADII', 3);
re = radii(1); rp = radii(3);
Rm = (re+rp)/2;
X = X * Rm;
Y = Y * Rm;
Z = Z * Rm;

lat_rad = lander.lat*cspice_rpd;
lon_rad = lander.lon*cspice_rpd;
alt = lander.alt;
flat = (re - rp) / re;
r_land_mcmf = cspice_pgrrec('Moon', lon_rad, lat_rad, alt, re, flat);
figure
plot3(r_sc_mcmf(1, :), r_sc_mcmf(2, :), r_sc_mcmf(3, :), 'DisplayName', 'S/C trajectory')
hold on
plot3(r_land_mcmf(1), r_land_mcmf(2), r_land_mcmf(3), 'o', 'MarkerSize', 10, 'MarkerEdgeColor','b', 'MarkerFaceColor','b', 'DisplayName', 'Lander')
text(r_land_mcmf(1), r_land_mcmf(2), r_land_mcmf(3)+400, 'Lander', 'FontSize', 20, 'Color', 'b', 'VerticalAlignment','bottom', 'HorizontalAlignment','left')
axis equal
surf(X, Y, Z, 'FaceColor','w')
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
%% Simulate measurements
Y_pred = pointing_lander(lander, et_vec, r_sc_mci, false);
[Y_pred, et_vec_vis] = visibility(Y_pred, et_vec, lander);

flag = length(et_vec)==length(et_vec_vis);
if flag
    disp('Lander always visible')
else
    disp('Lander not always visible')
end

range = Y_pred(3, :); % Extract ranges from predictions
mu_pos = [r_sc_mci; range].'; % Combine position and range for measurement
sigma_pos = 100 * 1e-3; % Measurement noise standard deviation
P_pos = sigma_pos^2 * eye(4); % Measurement noise covariance
pos_meas = mvnrnd(mu_pos, P_pos).'; % Generate noisy measurements

err = mu_pos.'-pos_meas;
z_err = reshape(err, 4*size(err, 2), 1);
mu_z = mean(z_err);
sigma_err = std(z_err);
z_err = (z_err-mu_z)/sigma_err;
nbins = ceil(sqrt(length(z_err)));
xx = linspace(-3, 3, 1e3);
yy = normpdf(xx);
figure
hold on
histogram(z_err, nbins, 'Normalization','pdf')
plot(xx, yy, 'r', 'LineWidth',1.5)
xline(-3, 'k--', 'LineWidth', 1.5)
xline(3, 'k--', 'LineWidth', 1.5)
xlabel('Values [$z_i$]')
ylabel('Normalized N. Samples')

boolean = (z_err < -3) | (z_err > 3);
outliers = sum(boolean);
perc = outliers / length(z_err) * 100
%% Estimate Orbiter state

r_meas = pos_meas(1:3, :);
P0 = diag([10, 1, 1, 0.001, 0.001, 0.001]); % Initial covariance
x0_guess = mvnrnd(xx0, P0).'; % Generate initial guess
R = sigma_pos^2 * eye(3); % Measurement noise covariance
[xx_hat, P] = UKF(r_meas, x0_guess, P0, et_vec, R, mu); % Run UKF

P_r = P(1:3, 1:3, :);
P_v = P(4:6, 4:6, :);
sigma_r = nan(length(et_vec), 1);
sigma_v = sigma_r;

for k = 1 : length(et_vec)

    Pr = P_r(:, :, k); Pv = P_v(:, :, k);
    sr = sqrt(max(eig(Pr)));
    sv = sqrt(max(eig(Pv)));

    sigma_r(k) = sr;
    sigma_v(k) = sv;

end

eps_est = xx_hat - xx_sc_m;
eps_est_rr = eps_est(1:3, :);
eps_est_vv = eps_est(4:6, :);
eps_est_r  = vecnorm(eps_est_rr);
eps_est_v  = vecnorm(eps_est_vv);

figure
subplot(1, 2, 1)
semilogy(t, eps_est_r)
hold on
semilogy(t, 3*sigma_r)
ylabel('$\varepsilon_r$, 3$\sigma_r$ [km]')
legend('$\varepsilon_r$', '3$\sigma_r$')

subplot(1, 2, 2)
semilogy(t, eps_est_v)
hold on
semilogy(t, 3*sigma_v)
ylabel('$\varepsilon_v$, 3$\sigma_v$ [km/s]')
legend('$\varepsilon_v$', '3$\sigma_v$')
%% Estimate Lander coordinates
xx0_land = [xx0; 78 * cspice_rpd; 15 * cspice_rpd]; % Include lander in initial state
P0 = diag([10, 1, 1, 0.001, 0.001, 0.001, 0.00001, 0.00001]); % Initial covariance
x0_guess = mvnrnd(xx0_land, P0).'; % Generate initial guess
R = sigma_pos^2 * eye(4); % Measurement noise covariance
[xx_hat, P] = UKF_land(pos_meas, x0_guess, P0, et_vec, R, mu); % Run UKF for orbiter + lander

P_r = P(1:3, 1:3, :);
P_v = P(4:6, 4:6, :);
P_lat = P(7, 7, :);
P_lon = P(8, 8, :);
sigma_r = nan(length(et_vec), 1);
sigma_v = sigma_r;
sigma_lat = sigma_r;
sigma_lon = sigma_r;

for k = 1 : length(et_vec)

    Pr = P_r(:, :, k); Pv = P_v(:, :, k);
    Plat = P_lat(:, :, k); Plon = P_lon(:, :, k);

    sr   = sqrt(max(eig(Pr)));
    sv   = sqrt(max(eig(Pv)));
    slat = sqrt(max(eig(Plat)));
    slon = sqrt(max(eig(Plon)));

    sigma_r(k)   = sr;
    sigma_v(k)   = sv;
    sigma_lat(k) = slat*cspice_dpr;
    sigma_lon(k) = slon*cspice_dpr;

end

nt = length(et_vec);
lat_ex = 78.229772;
lon_ex = 15.407786;
xx_hat(7:8, :) = xx_hat(7:8, :) * cspice_dpr;
lat_ref = lat_ex * ones(1, nt);
lon_ref = lon_ex * ones(1, nt);

eps_est = xx_hat - [xx_sc_m; lat_ref; lon_ref];
eps_est_rr = eps_est(1:3, :);
eps_est_vv = eps_est(4:6, :);
eps_est_Lat = eps_est(7, :);
eps_est_Lon = eps_est(8, :);

eps_est_r    = vecnorm(eps_est_rr);
eps_est_v    = vecnorm(eps_est_vv);
eps_est_lat  = abs(eps_est_Lat);
eps_est_lon  = abs(eps_est_Lon);

figure
subplot(2, 2, 1)
semilogy(t, eps_est_r)
hold on
semilogy(t, 3*sigma_r)
ylabel('$\varepsilon_r$, 3$\sigma_r$ [km]')
legend('$\varepsilon_r$', '3$\sigma_r$')

subplot(2, 2, 2)
semilogy(t, eps_est_v)
hold on
semilogy(t, 3*sigma_v)
ylabel('$\varepsilon_v$, 3$\sigma_v$ [km/s]')
legend('$\varepsilon_v$', '3$\sigma_v$')

subplot(2, 2, 3)
semilogy(t, eps_est_lat)
hold on
semilogy(t, 3*sigma_lat)
ylabel('$\varepsilon_{lat}$, 3$\sigma_{lat}$ [deg]')
legend('$\varepsilon_{lat}$', '3$\sigma_{lat}$')

subplot(2, 2, 4)
semilogy(t, eps_est_lon)
hold on
semilogy(t, 3*sigma_lon)
ylabel('$\varepsilon_{lon}$, 3$\sigma_{lon}$ [deg]')
legend('$\varepsilon_{lon}$', '3$\sigma_{lon}$')
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

function dy = two_body_rhs(~, y, mu)
% This function computes the time derivative of the state vector for the
% two-body problem.
%
% Inputs:
%   ~   - Unused time variable (required by MATLAB ODE solvers)
%   y   - State vector [x; y; z; vx; vy; vz], where:
%         x, y, z     - Position components (in km or a consistent unit)
%         vx, vy, vz  - Velocity components (in km/s or a consistent unit)
%   mu  - Standard gravitational parameter (G * M, where G is the gravitational 
%         constant and M is the mass of the central body, in km^3/s^2)
%
% Outputs:
%   dy  - Time derivative of the state vector [vx; vy; vz; ax; ay; az], where:
%         vx, vy, vz  - Velocity components
%         ax, ay, az  - Acceleration components due to gravity
    % Compute the distance from the central body
    r = sqrt(y(1)^2 + y(2)^2 + y(3)^2); % Magnitude of the position vector

    % Compute the time derivative of the state vector
    dy = [y(4);                      % dx/dt = vx
          y(5);                      % dy/dt = vy
          y(6);                      % dz/dt = vz
          -(mu / r^3) * y(1);        % ax = -mu * x / r^3
          -(mu / r^3) * y(2);        % ay = -mu * y / r^3
          -(mu / r^3) * y(3)];       % az = -mu * z / r^3
end

function Y_pred = pointing_lander(lander, et_vec, r_sc_m, flag)
% This function computes the pointing angles (Azimuth, Elevation) and range
% of a spacecraft relative to a lander on the Moon
%
% Inputs:
%   lander   - Struct containing the lander's parameters:
%              .lat   - Latitude of the lander (degrees)
%              .lon   - Longitude of the lander (degrees)
%              .alt   - Altitude of the lander (meters)
%              .name  - (Optional) Name of the lander for dynamic frame mode
%   et_vec   - Vector of ephemeris times (ET) at which computations are performed
%   r_sc_m   - Spacecraft position vectors in the inertial frame 'J2000'
%              [3 x n] array, where n is the number of time steps
%   flag     - Boolean flag to determine the mode of computation:
%              true  - topocentric frame based on estimated latitude and longitude
%              false - topocentric frame using SPICE kernels
%
% Outputs:
%   Y_pred   - Matrix containing computed pointing data:
%              [Azimuth; Elevation; Range], where:
%              Azimuth (degrees), Elevation (degrees), and Range (km)

frame1 = 'J2000'; % Define the inertial reference frame

if flag
    
    lat = lander.lat; % Lander latitude (degrees)
    lon = lander.lon; % Lander longitude (degrees)
    alt = lander.alt; % Lander altitude (meters)

    % Retrieve Moon's radii and compute flattening
    radii = cspice_bodvrd('Moon', 'RADII', 3); % [Equatorial radius, Polar radius]
    re = radii(1); rp = radii(3); % Equatorial and polar radii (km)
    flat = (re - rp) / re; % Flattening of the Moon

    % Convert lander coordinates to radians and altitude to km
    lat_rad = deg2rad(lat); lon_rad = deg2rad(lon); alt_km = alt / 1000;

    % Compute the lander's position in the Moon's body-fixed frame
    r_lander = cspice_pgrrec('Moon', lon_rad, lat_rad, alt_km, re, flat);

    % Define the Moon-fixed frame and compute transformations
    frame2 = 'IAU_Moon'; % Moon body-fixed frame
    R_IN2MF = cspice_pxform(frame1, frame2, et_vec); % Transform J2000 -> IAU_Moon
    R_MF2TOPO = cspice_eul2m(lat_rad - pi, pi - lon_rad, pi / 2, 2, 1, 2); % Transform IAU_Moon -> Topocentric

    % Initialize array for spacecraft position in the topocentric frame
    r_sc_TOPO = nan(size(r_sc_m));

    % Loop through ephemeris times to compute topocentric positions
    for k = 1:length(et_vec)
        r_sc_mf = R_IN2MF(:, :, k) * r_sc_m(:, k); % Spacecraft in Moon-fixed frame
        r_rel = r_sc_mf - r_lander; % Relative position to the lander
        r_sc_TOPO(:, k) = R_MF2TOPO * r_rel; % Transform to topocentric frame
    end

else
    
    name = lander.name; % Lander's name for SPICE kernel
    frame2 = strcat(name, '_TOPO'); % Topocentric frame defined by SPICE kernel

    % Compute the lander's position in the inertial frame
    r_lander = cspice_spkpos(name, et_vec, frame1, 'NONE', 'Moon');

    % Compute the relative position of the spacecraft to the lander
    r_rel = r_sc_m - r_lander;

    % Transform spacecraft position to the topocentric frame
    R = cspice_pxform(frame1, frame2, et_vec); % Transform J2000 -> Topocentric
    r_sc_TOPO = nan(size(r_sc_m));
    for k = 1:length(et_vec)
        r_sc_TOPO(:, k) = R(:, :, k) * r_rel(:, k); % Apply topocentric transformation
    end

end

% Compute the range, azimuth, and elevation from topocentric positions
[range, Az, El] = cspice_reclat(r_sc_TOPO); % Convert to spherical coordinates
Az = rad2deg(Az); % Convert azimuth to degrees
El = rad2deg(El); % Convert elevation to degrees

% Output the computed pointing data
Y_pred = [Az; El; range];
end

function [Y_pred, et] = visibility(Y_pred, et, station)
% This function filters the visibility data based on the elevation angle 
% of a ground station. It removes data points where the elevation angle 
% is below the minimum allowable value.
%
% Inputs:
%   Y_pred   - Matrix containing predicted visibility data:
%              [Azimuth; Elevation; Range] for each time step
%              Azimuth (degrees), Elevation (degrees), Range (km)
%   et       - Vector of ephemeris times corresponding to each prediction
%   station  - Struct containing the station's parameters:
%              .El_min - Minimum elevation angle for visibility (degrees)
%
% Outputs:
%   Y_pred   - Filtered visibility data containing only valid points where 
%              Elevation > El_min
%   et       - Filtered vector of ephemeris times corresponding to valid
%              visibility points

% Extract the minimum allowable elevation angle for the station
El_min = station.El_min;

% Extract the elevation angles from the visibility data
El = Y_pred(2, :); % Elevation (degrees)

% Create a logical array for points meeting the visibility condition
bool = double(El > El_min); % 1 for valid points, 0 for invalid points

% Identify indices where the visibility condition is satisfied
indeces = bool ~= 0; % Logical index of valid data points

% Filter the ephemeris times and visibility data using valid indices
et = et(indeces);             % Keep only valid ephemeris times
Y_pred = Y_pred(:, indeces);  % Keep only valid visibility data
end

function [X, Wm, Wc] = sigma_points(x0, P, alpha, beta)
% This function computes the sigma points and associated weights for the 
% Unscented Transform (UT) based on the input state and covariance.
%
% Inputs:
%   x0    - Initial state vector [n x 1]
%   P     - Initial covariance matrix [n x n]
%   alpha - Scaling parameter for controlling the spread of sigma points
%   beta  - Parameter
%
% Outputs:
%   X     - Sigma points matrix [n x (2n + 1)], where each column is a sigma point
%   Wm    - Weight vector for computing the mean [1 x (2n + 1)]
%   Wc    - Weight vector for computing the covariance [1 x (2n + 1)]

% Number of state dimensions
n = length(x0);

% Compute scaling parameter lambda
lambda = n * (alpha^2 - 1);

% Compute the square root scaling factor for sigma points
gamma = sqrt(n + lambda);

% Perform Cholesky decomposition of the covariance matrix
try
    L = chol(P, 'lower'); % Lower triangular matrix such that P = L * L'
catch
    error('Initial covariance matrix P is not positive definite.');
end

% Compute the scaled sigma point spread matrix
B = gamma * L; % Scale the Cholesky factor by gamma

% Define the mean and covariance weights for the sigma points
Wm = [lambda / (n + lambda), repmat(1 / (2 * (n + lambda)), 1, 2 * n)];
Wc = [lambda / (n + lambda) + (1 - alpha^2 + beta), repmat(1 / (2 * (n + lambda)), 1, 2 * n)];

% Generate sigma points
X = [x0, x0 + B, x0 - B];

end

function [xk_minus, P_minus, yk_hat, Xprop, Y] = prediction_step(Xprec, Wm, Wc, tspan, mu)
% This function performs the prediction step of the Unscented Kalman Filter (UKF).
% It propagates the sigma points through the dynamics, computes the predicted mean 
% and covariance, and maps sigma points to the measurement space.
%
% Inputs:
%   Xprec  - Matrix of sigma points from the previous step [n x (2n + 1)]
%   Wm     - Weight vector for mean computation [1 x (2n + 1)]
%   Wc     - Weight vector for covariance computation [1 x (2n + 1)]
%   tspan  - Time span for propagation [1 x 2] (start and end time)
%   mu     - Gravitational parameter (km^3/s^2)
%
% Outputs:
%   xk_minus - Predicted mean state vector [n x 1]
%   P_minus  - Predicted covariance matrix [n x n]
%   yk_hat   - Predicted measurement mean (position components) [3 x 1]
%   Xprop    - Propagated sigma points [n x (2n + 1)]
%   Y        - Sigma points mapped to the measurement space [3 x (2n + 1)]

% Define ODE solver options
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Define the two-body dynamics function
odefun = @(t, x) two_body_rhs(t, x, mu);

% Dimensions of the state and sigma points
n = size(Xprec, 1); % Number of state variables
m = size(Xprec, 2); % Number of sigma points

% Initialize predicted covariance matrix and propagated sigma points
P_minus = zeros(n);          % Predicted covariance [n x n]
Xprop = nan(size(Xprec));    % Propagated sigma points [n x (2n + 1)]

% Propagate each sigma point through the dynamics
for k = 1:m
    % Extract the current sigma point
    x0 = Xprec(:, k);
    
    % Numerically integrate the two-body equations of motion
    [~, xx] = ode78(odefun, tspan, x0, options);
    
    % Store the propagated sigma point (final state from integration)
    Xprop(:, k) = xx(end, :).';
end

% Compute the predicted mean state vector
xk_minus = Xprop * Wm';

% Compute the predicted covariance matrix
for k = 1:m
    % Compute the deviation of each sigma point from the mean
    diff = Xprop(:, k) - xk_minus;
    
    % Update the covariance using the weighted outer product of deviations
    P_minus = P_minus + Wc(k) * (diff * diff.');
end

% Map the propagated sigma points to the measurement space (position only)
Y = Xprop(1:3, :); % Extract position components [3 x (2n + 1)]

% Compute the predicted measurement mean
yk_hat = Y * Wm';

end

function [xk_plus, P_plus] = update_step(Wc, R, Y, yk_hat, Xprop, xk_minus, Pk_minus, yk)
% This function performs the update step of the Unscented Kalman Filter (UKF).
% It incorporates new measurements to refine the predicted state and covariance.
%
% Inputs:
%   Wc       - Covariance weights for the sigma points [1 x (2n + 1)]
%   R        - Measurement noise covariance matrix [ny x ny]
%   Y        - Sigma points in the measurement space [ny x (2n + 1)]
%   yk_hat   - Predicted measurement mean vector [ny x 1]
%   Xprop    - Propagated sigma points in the state space [nx x (2n + 1)]
%   xk_minus - Predicted mean state vector [nx x 1]
%   Pk_minus - Predicted covariance matrix [nx x nx]
%   yk       - Measurement vector [ny x 1]
%
% Outputs:
%   xk_plus  - Updated mean state vector [nx x 1]
%   P_plus   - Updated covariance matrix [nx x nx]

% Dimensions of state and measurement vectors
nx = length(xk_minus); % Number of state variables
ny = length(yk_hat);   % Number of measurement variables

% Initialize measurement covariance (Pyy) and cross-covariance (Pxy)
Pyy = R;                     % Start with measurement noise covariance
Pxy = zeros(nx, ny);         % Initialize cross-covariance matrix

% Compute Pyy (measurement covariance) and Pxy (cross-covariance)
for k = 1:size(Xprop, 2)
    % Deviation of the sigma points in state space
    diffx = Xprop(:, k) - xk_minus;
    
    % Deviation of the sigma points in measurement space
    diffy = Y(:, k) - yk_hat;
    
    % Update measurement covariance
    Pyy = Pyy + Wc(k) * (diffy * diffy.');
    
    % Update cross-covariance
    Pxy = Pxy + Wc(k) * (diffx * diffy.');
end

% Compute the Kalman gain
Kk = Pxy / Pyy;

% Update the mean state vector
xk_plus = xk_minus + Kk * (yk - yk_hat);

% Update the state covariance matrix
P_plus = Pk_minus - Kk * Pyy * Kk.';

end

function [xx_hat, P] = UKF(meas, xx0, P0, tvec, R, mu)
% This function implements the Unscented Kalman Filter (UKF) to estimate
% the state of a system based on measurements over time.
%
% Inputs:
%   meas   - Measurement matrix [ny x nt], where ny is the number of
%            measurements per time step and nt is the number of time steps
%   xx0    - Initial state vector [nx x 1]
%   P0     - Initial covariance matrix [nx x nx]
%   tvec   - Time vector [1 x nt] with the time steps for propagation
%   R      - Measurement noise covariance matrix [ny x ny]
%   mu     - Gravitational parameter (km^3/s^2)
%
% Outputs:
%   xx_hat - State estimates at each time step [nx x nt]
%   P      - State covariance matrices at each time step [nx x nx x nt]

% Parameters for the Unscented Transform
alpha = 0.01;  % Scaling parameter to control the spread of sigma points
beta = 2;      % Optimal for Gaussian distributions

% Initialize the state and covariance for the first time step
xprec = xx0;          % Prior state estimate
Pprec = P0;           % Prior covariance estimate
tprec = tvec(1);      % Initial time

% Dimensions
nx = length(xx0);     % Number of state variables
nt = length(tvec);    % Number of time steps

% Initialize output arrays
xx_hat = nan(nx, nt);     % State estimates
P = nan(nx, nx, nt);      % Covariance matrices

% Set initial state and covariance
xx_hat(:, 1) = xx0;
P(:, :, 1) = P0;

% Iterate through time steps
for k = 2:nt
    % Current time and time span for propagation
    tk = tvec(k);          % Current time
    tspan = [tprec, tk];   % Time interval for propagation
    
    % Current measurement
    yk = meas(:, k);       % Measurement at time tk

    % Generate sigma points for the prior state and covariance
    [Xprec, Wm, Wc] = sigma_points(xprec, Pprec, alpha, beta);

    % Prediction step: propagate sigma points through dynamics
    [xk_minus, P_minus, yk_hat, Xprop, Y] = prediction_step(Xprec, Wm, Wc, tspan, mu);

    % Update step: incorporate the measurement to refine the state estimate
    [xk_plus, P_plus] = update_step(Wc, R, Y, yk_hat, Xprop, xk_minus, P_minus, yk);

    % Ensure covariance matrix symmetry (numerical stability)
    P_plus = (P_plus + P_plus.') / 2;

    % Store the updated state and covariance
    xx_hat(:, k) = xk_plus;     % Store state estimate
    P(:, :, k) = P_plus;        % Store covariance estimate

    % Prepare for the next iteration
    tprec = tk;                 % Update time
    xprec = xk_plus;            % Update state
    Pprec = P_plus;             % Update covariance
end

end

function [xk_minus, P_minus, yk_hat, Xprop, Y] = prediction_step_land(Xprec, Wm, Wc, tspan, mu)
% This function performs the prediction step of the Unscented Kalman Filter (UKF) for a
% system involving both an orbiter and a lander on the Moon's surface.
%
% Inputs:
%   Xprec   - Matrix of sigma points from the previous step [n x (2n + 1)], where:
%             - First 6 rows correspond to the orbiter state (position and velocity)
%             - Last 2 rows correspond to the lander's latitude and longitude
%   Wm      - Weight vector for mean computation [1 x (2n + 1)]
%   Wc      - Weight vector for covariance computation [1 x (2n + 1)]
%   tspan   - Time span for propagation [1 x 2] (start and end time)
%   mu      - Gravitational parameter of the Moon (km^3/s^2)
%
% Outputs:
%   xk_minus - Predicted mean state vector [n x 1]
%   P_minus  - Predicted covariance matrix [n x n]
%   yk_hat   - Predicted measurement mean (orbiter position and range) [4 x 1]
%   Xprop    - Propagated sigma points [n x (2n + 1)]
%   Y        - Sigma points mapped to the measurement space [4 x (2n + 1)]

% ODE solver options
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

% Define the dynamics function for the two-body problem
odefun = @(t, x) two_body_rhs(t, x, mu);

% Retrieve Moon's radii and compute flattening
radii = cspice_bodvrd('Moon', 'RADII', 3); % [Equatorial radius, Polar radius]
re = radii(1); % Equatorial radius (km)
rp = radii(3); % Polar radius (km)
flat = (re - rp) / re; % Flattening factor

% Dimensions of the state and sigma points
n = size(Xprec, 1); % Number of state variables
m = size(Xprec, 2); % Number of sigma points

% Initialize predicted covariance matrix and propagated sigma points
P_minus = zeros(n);         % Predicted covariance [n x n]
Xprop = nan(size(Xprec));   % Propagated sigma points [n x (2n + 1)]
Y = nan(4, m);              % Measurement sigma points [4 x (2n + 1)]

% Separate sigma points for orbiter and lander states
Xprec_orb = Xprec(1:6, :);  % Orbiter state (position and velocity)
Xprec_land = Xprec(7:8, :); % Lander state (latitude and longitude)

% Propagate each sigma point through the dynamics
for k = 1:m
    % Orbiter propagation
    x0 = Xprec_orb(:, k);        % Initial state of the orbiter
    Lat = Xprec_land(1, k);      % Lander latitude (radians)
    Lon = Xprec_land(2, k);      % Lander longitude (radians)
    
    % Numerically integrate the orbiter's motion
    [~, xx] = ode78(odefun, tspan, x0, options);
    xf = xx(end, :).';           % Final state of the orbiter

    % Transform orbiter position to the Moon-fixed frame
    r_orb = xf(1:3);             % Orbiter position in J2000
    R_IN2MF = cspice_pxform('J2000', 'IAU_Moon', tspan(end)); % J2000 -> Moon-fixed
    r_orb_mf = R_IN2MF * r_orb;  % Orbiter position in Moon-fixed frame

    % Compute the lander's position in the Moon-fixed frame
    r_land = cspice_georec(Lon, Lat, 0, re, flat); % Lander position in Moon-fixed

    % Compute the relative position and range
    r_rel = r_orb_mf - r_land;   % Relative position
    range = norm(r_rel);         % Range (distance)

    % Store sigma points in the measurement space
    Y(:, k) = [r_orb; range];    % Orbiter position and range

    % Propagate sigma points in the state space
    Xprop(:, k) = [xf; Lat; Lon]; % Append latitude and longitude to propagated state
end

% Compute the predicted mean state vector
xk_minus = Xprop * Wm'; % Weighted mean of propagated sigma points

% Compute the predicted measurement mean
yk_hat = Y * Wm'; % Weighted mean of measurements

% Compute the predicted covariance matrix
for k = 1:m
    % Compute the deviation of each sigma point from the mean
    diff = Xprop(:, k) - xk_minus;

    % Update the covariance using the weighted outer product of deviations
    P_minus = P_minus + Wc(k) * (diff * diff.');
end

end

function [xx_hat, P] = UKF_land(meas, xx0, P0, tvec, R, mu)
% This function implements the Unscented Kalman Filter (UKF) for a system 
% involving both an orbiter and a lander. It estimates the state of the 
% system based on measurements over time.
%
% Inputs:
%   meas   - Measurement matrix [ny x nt], where:
%            ny is the number of measurements per time step
%            nt is the number of time steps
%   xx0    - Initial state vector [nx x 1], including:
%            - Orbiter state (position and velocity)
%            - Lander state (latitude and longitude)
%   P0     - Initial covariance matrix [nx x nx]
%   tvec   - Time vector [1 x nt] containing time steps for propagation
%   R      - Measurement noise covariance matrix [ny x ny]
%   mu     - Gravitational parameter of the Moon (km^3/s^2)
%
% Outputs:
%   xx_hat - State estimates at each time step [nx x nt]
%   P      - State covariance matrices at each time step [nx x nx x nt]

% Parameters for the Unscented Transform
alpha = 0.01; % Scaling parameter to control sigma point spread
beta = 2;     % Parameter for prior knowledge of Gaussian distributions

% Initialization
xprec = xx0;           % Prior state estimate
Pprec = P0;            % Prior covariance estimate
tprec = tvec(1);       % Initial time

% Dimensions
nx = length(xx0);      % Number of state variables
nt = length(tvec);     % Number of time steps

% Allocate memory for outputs
xx_hat = nan(nx, nt);  % State estimates
P = nan(nx, nx, nt);   % Covariance matrices

% Set initial state and covariance
xx_hat(:, 1) = xx0;
P(:, :, 1) = P0;

% Loop through time steps
for k = 2:nt
    % Current time and time span for propagation
    tk = tvec(k);         % Current time
    tspan = [tprec, tk];  % Time interval for propagation
    
    % Current measurement
    yk = meas(:, k);      % Measurement at time tk

    % Generate sigma points for the prior state and covariance
    [Xprec, Wm, Wc] = sigma_points(xprec, Pprec, alpha, beta);

    % Prediction step: propagate sigma points through dynamics
    [xk_minus, P_minus, yk_hat, Xprop, Y] = prediction_step_land(Xprec, Wm, Wc, tspan, mu);

    % Update step: incorporate the measurement to refine the state estimate
    [xk_plus, P_plus] = update_step(Wc, R, Y, yk_hat, Xprop, xk_minus, P_minus, yk);

    % Ensure covariance symmetry for numerical stability
    P_plus = (P_plus + P_plus.') / 2;

    % Store the updated state and covariance
    xx_hat(:, k) = xk_plus; % Store state estimate
    P(:, :, k) = P_plus;    % Store covariance estimate

    % Update prior values for the next iteration
    tprec = tk;             % Update time
    xprec = xk_plus;        % Update state
    Pprec = P_plus;         % Update covariance
end

end
