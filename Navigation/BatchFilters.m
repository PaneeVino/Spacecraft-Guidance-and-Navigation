clear; cspice_kclear; close all; clc;
plotStyle;

% addpath("kernels/")
% addpath("../mice/src/mice/")
% addpath("../mice/lib/")
% addpath("sgp4/")
% addpath("tle/")
% % 
cspice_kclear
cspice_furnsh('assignment02.tm')

sigma_Az = 125*1e-3;
sigma_El = 125*1e-3;
sigma_range = 1e-2;
d = [1/sigma_Az; 1/sigma_El; 1/sigma_range];
Wm = diag(d);
kourou   = struct('name', 'Kourou', 'frequency', 60, 'El_min', 6, ...
    'Wm', Wm,  'cost', 30000);
troll    = struct('name', 'Troll', 'frequency', 30, 'El_min', 0, ...
    'Wm', Wm, 'cost', 35000);
svalbard = struct('name', 'Svalbard', 'frequency', 60,'El_min', 8, ...
    'Wm', Wm, 'cost', 35000);

%% Visibility Windows
% Constants and Inputs
SAT_ID = 36036;                    % Satellite ID (SMOS)
filename = "36036.3le";            % TLE file
whichconst = 72;                   % Gravitational constant set (WGS-72)
[ satrec, ~, ~ ] = read_3LE(SAT_ID, filename, whichconst); % Read TLE data
mu = satrec.mu;                    % Gravitational parameter

% Extract and display TLE epoch information
[year, mon, day, hr, min, sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year, mon, day, hr, min, sec]);
etref = cspice_str2et(sat_epoch_str); % Convert to ephemeris time (ET)
fprintf('Satellite num ID: %d\n', satrec.satnum);
fprintf('TLE reference epoch: UTC %s\n', sat_epoch_str);

% Compute satellite state at epoch using SGP4
[satrec, rteme, vteme] = sgp4(satrec, 0.0); % TEME position and velocity

% Convert TEME state to osculating elements
elts = cspice_oscelt([rteme; vteme], etref, mu); % Osculating elements
fprintf('*** Osculating orbital elements ***\n');
fprintf('SMA   [km]:  %.5f\n', elts(1)/(1 - elts(2))); % Semi-major axis
fprintf('ECC        :  %.8f\n', elts(2));              % Eccentricity
fprintf('INC  [deg]: %.5f\n', elts(3) * cspice_dpr()); % Inclination
fprintf('RAAN [deg]: %.5f\n', elts(4) * cspice_dpr()); % RAAN
fprintf('ARGP [deg]: %.5f\n', elts(5) * cspice_dpr()); % Argument of perigee
fprintf('M.AN [deg]: %.5f\n', elts(6) * cspice_dpr()); % Mean anomaly

% Nutation and precession parameters
arcsec2rad = pi / (180 * 3600);
ddpsi = -0.114761 * arcsec2rad; % Nutation in longitude
ddeps = -0.007531 * arcsec2rad; % Nutation in obliquity
ttt = cspice_unitim(0.0, 'ET', 'TDT') / cspice_jyear() / 100; % Centuries since J2000

% Convert TEME to ECI state
[reci, veci] = teme2eci(rteme, vteme, [0.0; 0.0; 0.0], ttt, ddpsi, ddeps);
xxref = [reci; veci]; % Reference state in ECI

% Time bounds for visibility computation
epoch_0 = '2024-11-18T20:30:00.000';
epoch_f = '2024-11-18T22:15:00.000';
et0 = cspice_str2et(epoch_0);
etf = cspice_str2et(epoch_f);

% Propagate orbit and compute visibility for each ground station
[reci_k, ~, et_vec_k] = two_body_prop(xxref, et0, etf, kourou, mu, etref);
[reci_t, ~, et_vec_t] = two_body_prop(xxref, et0, etf, troll, mu, etref);
[reci_s, ~, et_vec_s] = two_body_prop(xxref, et0, etf, svalbard, mu, etref);

% Compute antenna pointing data for each station
[Y_k_id, et_vec_k] = pointing_antenna(kourou, et_vec_k, reci_k);
[Y_t_id, et_vec_t] = pointing_antenna(troll, et_vec_t, reci_t);
[Y_s_id, et_vec_s] = pointing_antenna(svalbard, et_vec_s, reci_s);

% Compute visibility windows
[Y_k_id, et_vec_k] = visibility(Y_k_id, et_vec_k, kourou);
[Y_t_id, et_vec_t] = visibility(Y_t_id, et_vec_t, troll);
[Y_s_id, et_vec_s] = visibility(Y_s_id, et_vec_s, svalbard);

% Extract Azimuth and Elevation for plotting
Az_k_id = Y_k_id(1, :);
Az_t_id = Y_t_id(1, :);
Az_s_id = Y_s_id(1, :);
El_k_id = Y_k_id(2, :);
El_t_id = Y_t_id(2, :);
El_s_id = Y_s_id(2, :);

% Convert times for plotting
t_k = datetime(cspice_timout(et_vec_k, 'YYYY-MM-DD HR:MN:SC.###'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
t_t = datetime(cspice_timout(et_vec_t, 'YYYY-MM-DD HR:MN:SC.###'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
t_s = datetime(cspice_timout(et_vec_s, 'YYYY-MM-DD HR:MN:SC.###'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');

% Plot Azimuth and Elevation for all stations
figure
subplot(2, 3, 1), plot(t_k, Az_k_id, 'o', 'Markersize', 4, 'MarkerEdgeColor', [0, 0.4470, 0.7410], 'MarkerFaceColor', [0, 0.4470, 0.7410]);
ylabel('Az [deg]'), title('Kourou');
subplot(2, 3, 2), plot(t_t, Az_t_id, 'o', 'Markersize', 4, 'MarkerEdgeColor', [0.8500, 0.3250, 0.0980], 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
ylabel('Az [deg]'), title('Troll');
subplot(2, 3, 3), plot(t_s, Az_s_id, 'o', 'Markersize', 4, 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880], 'MarkerFaceColor', [0.4660, 0.6740, 0.1880]);
ylabel('Az [deg]'), title('Svalbard');
subplot(2, 3, 4), plot(t_k, El_k_id, 'o', 'Markersize', 4, 'MarkerEdgeColor', [0, 0.4470, 0.7410], 'MarkerFaceColor', [0, 0.4470, 0.7410]);
ylabel('El [deg]'), title('Kourou');
subplot(2, 3, 5), plot(t_t, El_t_id, 'o', 'Markersize', 4, 'MarkerEdgeColor', [0.8500, 0.3250, 0.0980], 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
ylabel('El [deg]'), title('Troll');
subplot(2, 3, 6), plot(t_s, El_s_id, 'o', 'Markersize', 4, 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880], 'MarkerFaceColor', [0.4660, 0.6740, 0.1880]);
ylabel('El [deg]'), title('Svalbard');

% Display visibility windows
vis_k = [t_k(1) t_k(end)];
vis_t = [t_t(1) t_t(end)];
vis_s = [t_s(1) t_s(end)];
disp(vis_k);
disp(vis_t);
disp(vis_s);


%% Simulate Measurements

% Constants and Inputs
SAT_ID = 36036;                    % Satellite ID (SMOS)
filename = "36036.3le";            % TLE file
whichconst = 72;                   % Gravitational constant set (WGS-72)

% Read TLE and extract satellite parameters
[ satrec, longstr1, longstr2 ] = read_3LE(SAT_ID, filename, whichconst); % Read TLE data
mu = satrec.mu;                    % Gravitational parameter
satnum = satrec.satnum;            % Satellite number

% Extract and display TLE epoch information
[year, mon, day, hr, min, sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year, mon, day, hr, min, sec]);
sat_epoch_et = cspice_str2et(sat_epoch_str); % Convert epoch to ephemeris time (ET)
fprintf('Satellite num ID: %d\n', satnum);
fprintf('TLE reference epoch: UTC %s\n', sat_epoch_str);

% Compute satellite state at epoch using SGP4
[satrec, rteme, vteme] = sgp4(satrec, 0.0); % TEME position and velocity

% Convert TEME state to osculating elements for reference
elts = cspice_oscelt([rteme; vteme], sat_epoch_et, mu);

% Propagate satellite state using SGP4 for each ground station's visibility window
[reci_k, veci_k, et_vec_k] = sgp4_prop(et_vec_k(1), et_vec_k(end), kourou, sat_epoch_et, satrec);
[reci_t, ~, et_vec_t] = sgp4_prop(et_vec_t(1), et_vec_t(end), troll, sat_epoch_et, satrec);
[reci_s, ~, et_vec_s] = sgp4_prop(et_vec_s(1), et_vec_s(end), svalbard, sat_epoch_et, satrec);

% Compute pointing data (Azimuth, Elevation) for each ground station
[Y_k_id, et_vec_k] = pointing_antenna(kourou, et_vec_k, reci_k);
[Y_t_id, et_vec_t] = pointing_antenna(troll, et_vec_t, reci_t);
[Y_s_id, et_vec_s] = pointing_antenna(svalbard, et_vec_s, reci_s);

% Simulate noisy measurements for each station (Add Gaussian noise)
Y_k = simulate_measurements(Y_k_id, kourou);
Y_t = simulate_measurements(Y_t_id, troll);
Y_s = simulate_measurements(Y_s_id, svalbard);

% Compute visibility windows for each station (filter data below elevation threshold)
[Y_k, et_vec_k] = visibility(Y_k, et_vec_k, kourou);
[Y_t, et_vec_t] = visibility(Y_t, et_vec_t, troll);
[Y_s, et_vec_s] = visibility(Y_s, et_vec_s, svalbard);

% Extract Azimuth, Elevation, and Range for plotting
Az_k = Y_k(1, :); Az_t = Y_t(1, :); Az_s = Y_s(1, :);
El_k = Y_k(2, :); El_t = Y_t(2, :); El_s = Y_s(2, :);
Rng_k = Y_k(3, :); Rng_t = Y_t(3, :); Rng_s = Y_s(3, :);

% Convert ephemeris times to datetime for plotting
t_k = datetime(cspice_timout(et_vec_k, 'YYYY-MM-DD HR:MN:SC.###'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
t_t = datetime(cspice_timout(et_vec_t, 'YYYY-MM-DD HR:MN:SC.###'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
t_s = datetime(cspice_timout(et_vec_s, 'YYYY-MM-DD HR:MN:SC.###'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');

% Plot Azimuth, Elevation, and Range for all stations
figure
% Azimuth plots
subplot(3, 3, 1), plot(t_k, Az_k, 'o', 'Markersize', 4, 'MarkerEdgeColor', [0, 0.4470, 0.7410], 'MarkerFaceColor', [0, 0.4470, 0.7410]);
ylabel('Az [deg]'), title('Kourou');
subplot(3, 3, 2), plot(t_t, Az_t, 'o', 'Markersize', 4, 'MarkerEdgeColor', [0.8500, 0.3250, 0.0980], 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
ylabel('Az [deg]'), title('Troll');
subplot(3, 3, 3), plot(t_s, Az_s, 'o', 'Markersize', 4, 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880], 'MarkerFaceColor', [0.4660, 0.6740, 0.1880]);
ylabel('Az [deg]'), title('Svalbard');

% Elevation plots
subplot(3, 3, 4), plot(t_k, El_k, 'o', 'Markersize', 4, 'MarkerEdgeColor', [0, 0.4470, 0.7410], 'MarkerFaceColor', [0, 0.4470, 0.7410]);
ylabel('El [deg]'), title('Kourou');
subplot(3, 3, 5), plot(t_t, El_t, 'o', 'Markersize', 4, 'MarkerEdgeColor', [0.8500, 0.3250, 0.0980], 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
ylabel('El [deg]'), title('Troll');
subplot(3, 3, 6), plot(t_s, El_s, 'o', 'Markersize', 4, 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880], 'MarkerFaceColor', [0.4660, 0.6740, 0.1880]);
ylabel('El [deg]'), title('Svalbard');

% Range plots
subplot(3, 3, 7), plot(t_k, Rng_k, 'o', 'Markersize', 4, 'MarkerEdgeColor', [0, 0.4470, 0.7410], 'MarkerFaceColor', [0, 0.4470, 0.7410]);
ylabel('Range [km]'), title('Kourou');
subplot(3, 3, 8), plot(t_t, Rng_t, 'o', 'Markersize', 4, 'MarkerEdgeColor', [0.8500, 0.3250, 0.0980], 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
ylabel('Range [km]'), title('Troll');
subplot(3, 3, 9), plot(t_s, Rng_s, 'o', 'Markersize', 4, 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880], 'MarkerFaceColor', [0.4660, 0.6740, 0.1880]);
ylabel('Range [km]'), title('Svalbard');

% Save results in structured arrays for further analysis
kourou.Y = Y_k; kourou.tspan = [et0, et_vec_k];
troll.Y = Y_t; troll.tspan = [et0, et_vec_t];
svalbard.Y = Y_s; svalbard.tspan = [et0, et_vec_s];
%%
% Define ODE solver options
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-20);

% Time span for integration (from reference time to target time)
tspan = [etref, et0];

% Numerical integration for two-body problem
[~, xx] = ode78(@(t, x) two_body_rhs(t, x, mu), tspan, xxref, options);
x0_guess = xx(end, :).'; % Extract final state as initial guess

% Earth's mean radius (used for J2 perturbation calculations)
Re = cspice_bodvrd('Earth', 'RADII', 3); 
Re = (Re(3) + Re(2)) / 2; % Average radius (polar and equatorial)

% Numerical integration for two-body problem with J2 perturbation
[~, xx] = ode78(@(t, x) two_body_J2_rhs(t, x, mu, Re), tspan, xxref, options);
x0_guess_J2 = xx(end, :).'; % Extract final state for J2 model

% SGP4 propagation for the same time span
[reci_sgp4, veci_sgp4, ~] = sgp4_prop(etref, et0, kourou, sat_epoch_et, satrec);

%% Orbit Determination Using Kourou Station Only
% Combine position and velocity from SGP4 as the reference state
x0_sgp4 = [reci_sgp4(:, end); veci_sgp4(:, end)];

% Define the cost function for optimization (Kourou station only)
stations = kourou; % Station information
fun = @(x) costfunction(x, stations, false); % Cost function handle

% Set optimization options
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');

% Perform least-squares estimation
[x_only_k, resnorm, residual, exitflag, ~, ~, Jac] = lsqnonlin(fun, x0_guess, [], [], options);

% Compute covariance matrix from the Jacobian
B = (Jac.' * Jac) \ eye(size(Jac, 2));

% Estimate covariance and transform to orbital elements
P_only_k = resnorm / (length(residual) - length(x_only_k)) * B;
sigma_vec_only_k = linear_mapping(x_only_k, P_only_k, mu);
sigma_a_only_k = sigma_vec_only_k(1); % Semi-major axis uncertainty
sigma_i_only_k = rad2deg(sigma_vec_only_k(3)); % Inclination uncertainty in degrees

% Compute estimation errors in position and velocity
err = x_only_k - x0_sgp4; % Error between estimated and reference states
err_r = norm(err(1:3)); % Position error (norm)
err_v = norm(err(4:6)); % Velocity error (norm)

% Extract position and velocity covariance matrices
Pr = P_only_k(1:3, 1:3); % Position covariance matrix
Pv = P_only_k(4:6, 4:6); % Velocity covariance matrix

% Compute square root of trace of covariance matrices
sr_tr_Pr_only_k = sqrt(trace(Pr)); % sqrt(trace(Pr))
sr_tr_Pv_only_k = sqrt(trace(Pv)); % sqrt(trace(Pv))

disp('Only Kourou');
fprintf('ExitFlag: %d\n', exitflag); % Display optimization exit flag
disp('x - x0= ');
disp(err); % Display estimation errors
disp('err_r= ');
disp(err_r); % Display position error magnitude
disp('err_v= ');
disp(err_v); % Display velocity error magnitude
disp('P:');
disp(P_only_k); % Display covariance matrix
disp('sqrt(tr(Pr)):');
disp(sr_tr_Pr_only_k); % Display sqrt(trace(Pr))
disp('sqrt(tr(Pv)):');
disp(sr_tr_Pv_only_k); % Display sqrt(trace(Pv))
disp('sigma_a [km]:');
disp(sigma_a_only_k); % Display semi-major axis uncertainty
disp('sigma_i [deg]:');
disp(sigma_i_only_k); % Display inclination uncertainty

%% Orbit Determination Using All Stations

% Combine all ground stations
stations = [kourou; troll; svalbard];

% Define the cost function for optimization
fun = @(x) costfunction(x, stations, false);

% Set optimization options
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');

% Perform least-squares estimation using data from all stations
[x_all_stat, resnorm, residual, exitflag, ~, ~, Jac] = lsqnonlin(fun, x0_guess, [], [], options);

% Compute covariance matrix from the Jacobian
B = (Jac.' * Jac) \ eye(size(Jac, 2));

% Estimate covariance and transform to orbital elements
P_all_stat = resnorm / (length(residual) - length(x_all_stat)) * B;
sigma_vec_all_stat = linear_mapping(x_all_stat, P_all_stat, mu);
sigma_a_all_stat = sigma_vec_all_stat(1); % Semi-major axis uncertainty [km]
sigma_i_all_stat = rad2deg(sigma_vec_all_stat(3)); % Inclination uncertainty [deg]

% Compute estimation errors in position and velocity
err = x_all_stat - x0_sgp4; % Error between estimated and reference states
err_r = norm(err(1:3)); % Position error (norm) [km]
err_v = norm(err(4:6)); % Velocity error (norm) [km/s]

% Extract position and velocity covariance matrices
Pr = P_all_stat(1:3, 1:3); % Position covariance matrix
Pv = P_all_stat(4:6, 4:6); % Velocity covariance matrix

% Compute square root of trace of covariance matrices
sr_tr_Pr_all_stat = sqrt(trace(Pr)); % sqrt(trace(Pr)) for position uncertainty
sr_tr_Pv_all_stat = sqrt(trace(Pv)); % sqrt(trace(Pv)) for velocity uncertainty

disp('All stations');
fprintf('ExitFlag: %d\n', exitflag); % Display optimization exit flag
disp('x - x0= ');
disp(err); % Display estimation errors
disp('err_r= ');
disp(err_r); % Display position error magnitude
disp('err_v= ');
disp(err_v); % Display velocity error magnitude
disp('P:');
disp(P_all_stat); % Display covariance matrix
disp('sqrt(tr(Pr)):');
disp(sr_tr_Pr_all_stat); % Display sqrt(trace(Pr)) for position uncertainty
disp('sqrt(tr(Pv)):');
disp(sr_tr_Pv_all_stat); % Display sqrt(trace(Pv)) for velocity uncertainty
disp('sigma_a [km]:');
disp(sigma_a_all_stat); % Display semi-major axis uncertainty
disp('sigma_i [deg]:');
disp(sigma_i_all_stat); % Display inclination uncertainty

%% Orbit Determination Using All Stations with J2 Perturbation

% Combine all ground stations
stations = [kourou; troll; svalbard];

% Define the cost function for optimization with J2 perturbation
fun = @(x) costfunction(x, stations, true); % Third argument "true" enables J2 dynamics

% Set optimization options
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter');

% Perform least-squares estimation using data from all stations with J2 perturbation
[x_also_J2, resnorm, residual, exitflag, ~, ~, Jac] = lsqnonlin(fun, x0_guess_J2, [], [], options);

% Compute covariance matrix from the Jacobian
B = (Jac.' * Jac) \ eye(size(Jac, 2));

% Estimate covariance and transform to orbital elements
P_also_J2 = resnorm / (length(residual) - length(x_also_J2)) * B;
sigma_vec_also_J2 = linear_mapping(x_also_J2, P_also_J2, mu);
sigma_a_also_J2 = sigma_vec_also_J2(1); % Semi-major axis uncertainty [km]
sigma_i_also_J2 = rad2deg(sigma_vec_also_J2(3)); % Inclination uncertainty [deg]

% Compute estimation errors in position and velocity
err = x_also_J2 - x0_sgp4; % Error between estimated and reference states
err_r = norm(err(1:3)); % Position error magnitude [km]
err_v = norm(err(4:6)); % Velocity error magnitude [km/s]

% Extract position and velocity covariance matrices
Pr = P_also_J2(1:3, 1:3); % Position covariance matrix
Pv = P_also_J2(4:6, 4:6); % Velocity covariance matrix

% Compute square root of trace of covariance matrices
sr_tr_Pr_also_J2 = sqrt(trace(Pr)); % sqrt(trace(Pr)) for position uncertainty
sr_tr_Pv_also_J2 = sqrt(trace(Pv)); % sqrt(trace(Pv)) for velocity uncertainty

% Extract estimated position vector for visualization
mu_r = x_also_J2(1:3);

disp('Also J2');
fprintf('ExitFlag: %d\n', exitflag); % Display optimization exit flag
disp('x - x0= ');
disp(err); % Display estimation errors
disp('err_r= ');
disp(err_r); % Display position error magnitude
disp('err_v= ');
disp(err_v); % Display velocity error magnitude
disp('P:');
disp(P_also_J2); % Display covariance matrix
disp('sqrt(tr(Pr)):');
disp(sr_tr_Pr_also_J2); % Display sqrt(trace(Pr)) for position uncertainty
disp('sqrt(tr(Pv)):');
disp(sr_tr_Pv_also_J2); % Display sqrt(trace(Pv)) for velocity uncertainty
disp('sigma_a [km]:');
disp(sigma_a_also_J2); % Display semi-major axis uncertainty
disp('sigma_i [deg]:');
disp(sigma_i_also_J2); % Display inclination uncertainty

% Generate the 3σ confidence ellipsoid for position uncertainty
[X, Y, Z] = confidence_ellipsoid(Pr, mu_r, 3);

% Plot the ellipsoid and reference points
figure;
surf(X, Y, Z, 'FaceColor', 'cyan', 'EdgeColor', 'k', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.17); % Ellipsoid
hold on;

% Plot the estimated position and reference position (SGP4)
plot3(mu_r(1), mu_r(2), mu_r(3), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r'); % Estimated position
plot3(x0_sgp4(1), x0_sgp4(2), x0_sgp4(3), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b'); % Reference position

text(mu_r(1), mu_r(2), mu_r(3), '$\mathbf{\hat{r}}_0$', 'FontSize', 15, 'Color', 'k', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');
text(x0_sgp4(1), x0_sgp4(2), x0_sgp4(3), '$\mathbf{r}_{0,SGP4}$', 'FontSize', 15, 'Color', 'k', ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Interpreter', 'latex');

xlabel('X [km]', 'Interpreter', 'latex');
ylabel('Y [km]', 'Interpreter', 'latex');
zlabel('Z [km]', 'Interpreter', 'latex');

axis equal;

%% Trade-Off Analysis
% This script evaluates the trade-off between performance and cost for orbit determination
% using different combinations of ground stations. It computes the normalized uncertainties
% and performance metrics for semi-major axis and inclination.

% Extract orbital elements
SMA = elts(1) / (1 - elts(2)); % Semi-major axis [km]
ECC = elts(2);                % Eccentricity
INC = elts(3) * cspice_dpr(); % Inclination [deg]
RAAN = elts(4) * cspice_dpr(); % RAAN [deg]
ARGP = elts(5) * cspice_dpr(); % Argument of perigee [deg]
MAN = elts(6) * cspice_dpr();  % Mean anomaly [deg]

% Perform trade-off analysis for each station and combination
stations = kourou;
[cost_k, sigma_a_k, sigma_i_k] = trade_off_analysis(stations, x0_guess_J2);

stations = troll;
[cost_t, sigma_a_t, sigma_i_t] = trade_off_analysis(stations, x0_guess_J2);

stations = svalbard;
[cost_s, sigma_a_s, sigma_i_s] = trade_off_analysis(stations, x0_guess_J2);

stations = [kourou; troll];
[cost_kt, sigma_a_kt, sigma_i_kt] = trade_off_analysis(stations, x0_guess_J2);

stations = [kourou; svalbard];
[cost_ks, sigma_a_ks, sigma_i_ks] = trade_off_analysis(stations, x0_guess_J2);

stations = [troll; svalbard];
[cost_ts, sigma_a_ts, sigma_i_ts] = trade_off_analysis(stations, x0_guess_J2);

stations = [kourou; troll; svalbard];
[cost_kts, sigma_a_kts, sigma_i_kts] = trade_off_analysis(stations, x0_guess_J2);

% Normalize uncertainties for semi-major axis and inclination
sigma_a_vec_ad = [sigma_a_k; sigma_a_t; sigma_a_s; sigma_a_kt; sigma_a_ks; sigma_a_ts] / SMA;
sigma_i_vec_ad = [sigma_i_k; sigma_i_t; sigma_i_s; sigma_i_kt; sigma_i_ks; sigma_i_ts] / INC;

% Store costs for each configuration
cost_vec = [cost_k; cost_t; cost_s; cost_kt; cost_ks; cost_ts];

% Define colors for each station combination
Colors = [0, 0.4470, 0.7410;    % Blue
          0.8500, 0.3250, 0.0980; % Orange
          0.4660, 0.6740, 0.1880; % Green
          0.3010, 0.7450, 0.9330; % Light Blue
          0.4940, 0.1840, 0.5560; % Purple
          0.9290, 0.6940, 0.1250]; % Yellow

% Create a log-log plot of normalized uncertainties
figure;
loglog(sigma_a_vec_ad(1), sigma_i_vec_ad(1), 'o', 'MarkerSize', 10, ...
       'MarkerEdgeColor', Colors(1, :), 'MarkerFaceColor', Colors(1, :));
hold on;
for k = 2:length(sigma_a_vec_ad)
    color = Colors(k, :);
    loglog(sigma_a_vec_ad(k), sigma_i_vec_ad(k), 'o', 'MarkerSize', 10, ...
           'MarkerEdgeColor', color, 'MarkerFaceColor', color);
end

% Add legend and labels
legend('Only Kourou', 'Only Troll', 'Only Svalbard', 'Kourou + Troll', ...
       'Kourou + Svalbard', 'Troll + Svalbard', 'Location', 'northwest', ...
       'Orientation', 'vertical');
xlabel('$\overline{\sigma}_{aa}$', 'Interpreter', 'latex'); % Normalized semi-major axis uncertainty
ylabel('$\overline{\sigma}_{ii}$', 'Interpreter', 'latex'); % Normalized inclination uncertainty


% Compute a performance metric based on the distance from the origin
dist = sqrt(sigma_i_vec_ad.^2 + sigma_a_vec_ad.^2); % Euclidean distance
perf = 1 ./ dist; % Higher performance corresponds to smaller uncertainty

% Define station combinations and costs
station_combinations = ["Kourou", "Troll", "Svalbard", ...
                        "Kourou + Troll", "Kourou + Svalbard", "Troll + Svalbard"];
costs = [cost_k, cost_t, cost_s, cost_kt, cost_ks, cost_ts];

% Create a table summarizing performance and costs
results_table = table(station_combinations', perf, costs', ...
    'VariableNames', {'Stations', 'Perf', 'Cost'});

% Display the table in the command window
disp(results_table);
%%
% Define nutation parameters
ddpsi = -0.114761 * arcsec2rad; 
ddeps = -0.007531 * arcsec2rad;  

% Define the time interval for ground track propagation
epoch_0 = '2024-11-18T0:0:00.000';  % Start time in UTC
epoch_f = '2024-11-18T23:59:00.000'; % End time in UTC

% Convert time strings to ephemeris time (ET)
et0 = cspice_str2et(epoch_0);  % Start time in ET
etf = cspice_str2et(epoch_f);  % End time in ET

% Define the number of time points for propagation
npoints = 5000;  % Number of points to generate in the time vector
et_vec = linspace(et0, etf, npoints);  % Time vector from start to end

% Preallocate arrays for satellite position and velocity in ECI frame
reci = nan(3, npoints);  % Position in ECI frame
veci = nan(3, npoints);  % Velocity in ECI frame

% Loop through each time step to compute position and velocity
for i = 1:npoints
    % Compute time since TLE epoch in minutes
    tsince = (et_vec(i) - sat_epoch_et) / 60;  % Time since epoch in minutes

    % Propagate the satellite state using SGP4
    [~, rteme, vteme] = sgp4(satrec, tsince);

    % Compute centuries from J2000 epoch in TDT
    ttt = cspice_unitim(et_vec(i), 'ET', 'TDT') / cspice_jyear() / 100;

    % Convert position and velocity from TEME to ECI
    [reci(:, i), veci(:, i)] = ...
        teme2eci(rteme, vteme, [0.0; 0.0; 0.0], ttt, ddpsi, ddeps);
end

% Combine position and velocity vectors in ECI
xxeci = [reci; veci];

% Precompute the transformation matrix from J2000 to IAU_EARTH
Rot = cspice_pxform('J2000', 'IAU_EARTH', et_vec);

% Preallocate array for positions in IAU_EARTH frame
refix = nan(size(reci));

% Transform positions from ECI to IAU_EARTH frame
for k = 1:length(et_vec)
    refix(:, k) = Rot(:, :, k) * reci(:, k);
end

% Retrieve Earth radii and flattening factor from SPICE kernels
radii = cspice_bodvrd('EARTH', 'RADII', 3);  % Radii in km
flat = (radii(1) - radii(3)) / radii(1);     % Flattening factor

% Convert Cartesian coordinates to geodetic coordinates (lon, lat, alt)
[lon, lat, alt] = cspice_recgeo(refix, radii(1), flat);

% Convert longitude and latitude from radians to degrees
lon = lon * cspice_dpr;  % Convert to degrees
lat = lat * cspice_dpr;  % Convert to degrees

% Detect discontinuities in longitude due to wrapping
lon_diff = abs(diff(lon));                       % Compute differences between consecutive points
discontinuity_indices = find(lon_diff > 180);    % Find large jumps (>180°)

% Insert NaN at discontinuities to break the plot
lon(discontinuity_indices) = NaN;
lat(discontinuity_indices) = NaN;

% Load Earth map image
Im = flip(imread('Earth_NASA.jpg'));

% Plot the ground track on the Earth map
figure;
hold on;

% Plot the Earth map as a background
image([-180, 180], [-90, 90], Im);

% Plot the satellite ground track
plot(lon, lat, 'LineWidth', 1.5, 'Color', 'g');

% Plot specific locations with markers and labels
% Kourou
plot(-52.80466, 5.25144, 'o', 'MarkerSize', 7, ...
    'MarkerEdgeColor', [0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
text(-52.80466, 5.25144, 'Kourou', 'FontSize', 20, ...
    'Color', [0.8500, 0.3250, 0.0980], ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Troll
plot(2.536103, -72.011977, 'o', 'MarkerSize', 7, ...
    'MarkerEdgeColor', [0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
text(2.536103, -72.011977, 'Troll', 'FontSize', 20, ...
    'Color', [0.8500, 0.3250, 0.0980], ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Svalbard
plot(15.407786, 78.229772, 'o', 'MarkerSize', 7, ...
    'MarkerEdgeColor', [0.8500, 0.3250, 0.0980], ...
    'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
text(15.407786, 78.229772, 'Svalbard', 'FontSize', 20, ...
    'Color', [0.8500, 0.3250, 0.0980], ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

% Configure axis limits and labels
xlim([-170, 170]);
ylim([-90, 90]);
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');

%%
% Define constants and parameters
mu = cspice_bodvrd('Earth', 'GM', 1); % Earth's gravitational parameter [km^3/s^2]
Re = cspice_bodvrd('Earth', 'RADII', 3); % Earth's radii [km]
Re = (Re(3) + Re(2)) / 2; % Average equatorial and polar radius [km]
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); % ODE solver tolerances

% Define start and end epochs
epoch_0 = '2024-11-18T20:30:00.000'; % Start time
epoch_f = '2024-11-20T20:30:00.000'; % End time
et0 = cspice_str2et(epoch_0); % Convert start time to ephemeris time (ET)
etf = cspice_str2et(epoch_f); % Convert end time to ET

% Propagate satellite motion for three ground stations (Kourou, Troll, Svalbard)
[rr_k, vv_k, et_vec_k] = two_body_prop_J2(xxref, et0, etf, kourou, mu, Re, etref);
[rr_t, vv_t, et_vec_t] = two_body_prop_J2(xxref, et0, etf, troll, mu, Re, etref);
[rr_s, vv_s, et_vec_s] = two_body_prop_J2(xxref, et0, etf, svalbard, mu, Re, etref);

% Calculate pointing information for each ground station
[Y_k_id, et_vec_k] = pointing_antenna(kourou, et_vec_k, rr_k);
[Y_t_id, et_vec_t] = pointing_antenna(troll, et_vec_t, rr_t);
[Y_s_id, et_vec_s] = pointing_antenna(svalbard, et_vec_s, rr_s);

% Initialize visibility flags
flag_k = nan(size(et_vec_k));
flag_t = nan(size(et_vec_t));
flag_s = nan(size(et_vec_s));

coord_k = 1; % Kourou
coord_t = 2; % Troll
coord_s = 3; % Svalbard

% Check visibility for Kourou
for k = 1:length(et_vec_k)
    if Y_k_id(2, k) >= kourou.El_min % Check if elevation meets minimum requirement
        flag_k(k) = coord_k; % Assign visibility flag for Kourou
    end
end

% Check visibility for Troll
for k = 1:length(et_vec_t)
    if Y_t_id(2, k) >= troll.El_min % Check if elevation meets minimum requirement
        flag_t(k) = coord_t; % Assign visibility flag for Troll
    end
end

% Check visibility for Svalbard
for k = 1:length(et_vec_s)
    if Y_s_id(2, k) >= svalbard.El_min % Check if elevation meets minimum requirement
        flag_s(k) = coord_s; % Assign visibility flag for Svalbard
    end
end

% Convert ephemeris times to datetime for plotting
t_k = datetime(cspice_timout(et_vec_k, 'YYYY-MM-DD HR:MN:SC.###'), ...
    'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
t_t = datetime(cspice_timout(et_vec_t, 'YYYY-MM-DD HR:MN:SC.###'), ...
    'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
t_s = datetime(cspice_timout(et_vec_s, 'YYYY-MM-DD HR:MN:SC.###'), ...
    'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');

% Define colors for each station
Colors = [0, 0.4470, 0.7410; % Blue for Kourou
          0.8500, 0.3250, 0.0980; % Orange for Troll
          0.4660, 0.6740, 0.1880]; % Green for Svalbard

% Plot visibility flags over time
Markersize = 7;
figure;
hold on;
plot(t_k, flag_k, 'o', 'MarkerSize', Markersize, 'MarkerEdgeColor', Colors(1, :), ...
    'MarkerFaceColor', Colors(1, :)); % Kourou visibility
plot(t_t, flag_t, 'o', 'MarkerSize', Markersize, 'MarkerEdgeColor', Colors(2, :), ...
    'MarkerFaceColor', Colors(2, :)); % Troll visibility
plot(t_s, flag_s, 'o', 'MarkerSize', Markersize, 'MarkerEdgeColor', Colors(3, :), ...
    'MarkerFaceColor', Colors(3, :)); % Svalbard visibility

% Configure plot
ylim([0 4]); % Set y-axis limits
yticks(1:3); % Set y-axis ticks for ground stations
yticklabels({'Kourou', 'Troll', 'Svalbard'}); % Label ground stations
ytickangle(45); % Rotate y-axis tick labels for clarity
grid on;

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

function dy = two_body_J2_rhs(et, y, mu, Re)
% This function computes the time derivative of the state vector for the 
% two-body problem, including the effect of the J2 perturbation.
%
% Inputs:
%   et  - Ephemeris time (seconds past J2000; used for frame transformations)
%   y   - State vector [x; y; z; vx; vy; vz], where:
%         x, y, z     - Position components in the J2000 inertial frame
%         vx, vy, vz  - Velocity components in the J2000 inertial frame
%   mu  - Gravitational parameter (km^3/s^2) for the central body
%   Re  - Equatorial radius of the central body (in km)
%
% Outputs:
%   dy  - Time derivative of the state vector:
%         [dx/dt; dy/dt; dz/dt; dvx/dt; dvy/dt; dvz/dt], including J2 acceleration.
    % Gravitational harmonic coefficient for J2
    J2 = 0.0010826269;

    % Calculate the magnitude of the position vector
    r = sqrt(y(1)^2 + y(2)^2 + y(3)^2);

    % Initialize the derivative of the state vector
    dy = [y(4);                      % dx/dt = vx
          y(5);                      % dy/dt = vy
          y(6);                      % dz/dt = vz
          -(mu / r^3) * y(1);        % ax = -mu * x / r^3
          -(mu / r^3) * y(2);        % ay = -mu * y / r^3
          -(mu / r^3) * y(3)];       % az = -mu * z / r^3

    % Extract position components in the J2000 frame
    rr = y(1:3);

    % Transform position vector to the ITRF93 frame
    % (using ephemeris time for accurate rotation matrices)
    rotm = cspice_pxform('J2000', 'ITRF93', et);  % Rotation matrix from J2000 to ITRF93
    rr_itrf = rotm * rr;                          % Position vector in the ITRF93 frame

    % Recompute the norm in the ITRF93 frame
    r_itrf = norm(rr_itrf);

    % Compute the J2 perturbation acceleration in the ITRF93 frame
    aj2_itrf = 3/2 * mu * J2 * (Re / r_itrf)^2 * (rr_itrf / r_itrf^3) .* ...
               (5 * (rr_itrf(3) / r_itrf)^2 - [1; 1; 3]);

    % Transform aj2 acceleration back to the J2000 inertial frame
    rotm_back = cspice_pxform('ITRF93', 'J2000', et);  % Rotation matrix from ITRF93 to J2000
    aj2 = rotm_back * aj2_itrf;

    % Add the J2 acceleration to the velocity derivatives
    dy(4:6) = dy(4:6) + aj2;
end

function [rr, vv, et_vec] = two_body_prop(xxref, et0, etf, station, mu, etref)
% This function propagates the state vector [position; velocity] using the two-body 
% problem dynamics. It outputs the propagated position and velocity vectors at 
% specified time intervals and the corresponding ephemeris times.
%
% Inputs:
%   xxref   - Reference state vector [x; y; z; vx; vy; vz] at the reference time (in J2000 frame)
%   et0     - Initial ephemeris time (start of propagation, in seconds past J2000)
%   etf     - Final ephemeris time (end of propagation, in seconds past J2000)
%   station - Struct containing station parameters, including:
%             station.frequency: Measurement frequency (seconds)
%   mu      - Gravitational parameter (km^3/s^2) for the central body
%   etref   - Reference ephemeris time corresponding to `xxref` (in seconds past J2000)
%
% Outputs:
%   rr      - Propagated position vectors [x; y; z] (in km), each column corresponds to a time step
%   vv      - Propagated velocity vectors [vx; vy; vz] (in km/s), each column corresponds to a time step
%   et_vec  - Vector of ephemeris times (in seconds past J2000) for the propagation steps

    % Extract the measurement frequency (in seconds)
    freq = station.frequency;

    % Calculate the number of propagation points based on the measurement frequency
    npoints = round((etf - et0) / freq) + 1;

    % Initial propagation from reference time `etref` to the start time `et0`
    % Using high-precision numerical integration with ode78
    tspan = [etref, et0]; % Time span for propagation
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

    % Propagate the state from reference time to initial time
    [~, xx] = ode78(@(t, x) two_body_rhs(t, x, mu), tspan, xxref, options);
    xx0 = xx(end, :).'; % Extract the final state at `et0`

    % Generate the time vector for propagation between `et0` and `etf`
    et_vec = linspace(et0, etf, npoints);

    % Propagate the state from `et0` to `etf` at the specified intervals
    [~, xx] = ode78(@(t, x) two_body_rhs(t, x, mu), et_vec, xx0, options);
    xx = xx.';

    % Extract the propagated position and velocity components
    rr = xx(1:3, :); % Position vectors [x; y; z]
    vv = xx(4:6, :); % Velocity vectors [vx; vy; vz]
end

function [rr, vv, et_vec] = two_body_prop_J2(xxref, et0, etf, station, mu, Re, etref)
% This function propagates the state vector [position; velocity] using the two-body 
% problem dynamics with J2 perturbation. It outputs the propagated position 
% and velocity vectors at specified time intervals and the corresponding ephemeris times.
%
% Inputs:
%   xxref   - Reference state vector [x; y; z; vx; vy; vz] at the reference time (in J2000 frame)
%   et0     - Initial ephemeris time (start of propagation, in seconds past J2000)
%   etf     - Final ephemeris time (end of propagation, in seconds past J2000)
%   station - Struct containing station parameters, including:
%             station.frequency: Measurement frequency (seconds)
%   mu      - Gravitational parameter (km^3/s^2) for the central body
%   Re      - Equatorial radius of the central body (in km)
%   etref   - Reference ephemeris time corresponding to `xxref` (in seconds past J2000)
%
% Outputs:
%   rr      - Propagated position vectors [x; y; z] (in km), each column corresponds to a time step
%   vv      - Propagated velocity vectors [vx; vy; vz] (in km/s), each column corresponds to a time step
%   et_vec  - Vector of ephemeris times (in seconds past J2000) for the propagation steps

    % Extract the measurement frequency (in seconds)
    freq = station.frequency;

    % Calculate the number of propagation points based on the measurement frequency
    npoints = round((etf - et0) / freq) + 1;

    % Initial propagation from reference time `etref` to the start time `et0`
    % Using high-precision numerical integration with ode78
    tspan = [etref, et0]; % Time span for propagation
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); 

    % Propagate the state from reference time to initial time
    [~, xx] = ode78(@(t, x) two_body_J2_rhs(t, x, mu, Re), tspan, xxref, options);
    xx0 = xx(end, :).'; % Extract the final state at `et0`

    % Generate the time vector for propagation between `et0` and `etf`
    et_vec = linspace(et0, etf, npoints);

    % Propagate the state from `et0` to `etf` at the specified intervals
    [~, xx] = ode78(@(t, x) two_body_J2_rhs(t, x, mu, Re), et_vec, xx0, options);

    xx = xx.';

    % Extract the propagated position and velocity components
    rr = xx(1:3, :); % Position vectors [x; y; z]
    vv = xx(4:6, :); % Velocity vectors [vx; vy; vz]
end

function [reci, veci, et_vec] = sgp4_prop(et0, etf, station, sat_epoch_et, satrec)
% This function propagates the position and velocity vectors of the S/C
% using the SGP4 model and converts them from the TEME frame to the ECI frame.
%
% Inputs:
%   et0        - Initial ephemeris time (start of propagation, in seconds past J2000)
%   etf        - Final ephemeris time (end of propagation, in seconds past J2000)
%   station    - Struct containing station parameters
%   sat_epoch_et - Epoch time of the satellite in ephemeris time (ET) (in seconds past J2000)
%   satrec     - Structure containing the TLE parameters of the satellite
%
% Outputs:
%   reci      - Propagated position vectors in the ECI frame [x; y; z] (in km), each column corresponds to a time step
%   veci      - Propagated velocity vectors in the ECI frame [vx; vy; vz] (in km/s), each column corresponds to a time step
%   et_vec    - Vector of ephemeris times (in seconds past J2000) for the propagation steps

    % Conversion factor from arcseconds to radians
    arcsec2rad = pi / (180*3600);

    % Nutation parameters for converting from TEME to ECI at 18/11/2024
    ddpsi = -0.114752 * arcsec2rad;  % Nutation in longitude (rad)
    ddeps = -0.007529 * arcsec2rad; % Nutation in obliquity (rad)

    % Extract the measurement frequency from the station structure (in seconds)
    freq = station.frequency;

    % Calculate the number of propagation points based on the measurement frequency
    % This defines how many steps are taken from the initial to the final ephemeris time
    npoints = round((etf - et0) / freq) + 1;

    % Generate the vector of ephemeris times between `et0` and `etf` at the specified frequency
    et_vec = linspace(et0, etf, npoints);

    % Initialize the position and velocity matrices (3xN) to store the results
    reci = nan(3, npoints);  % Position vectors in ECI frame (in km)
    veci = nan(3, npoints);  % Velocity vectors in ECI frame (in km/s)

    % Propagate the satellite state for each ephemeris time step
    for i = 1:npoints
        % Calculate the time since the TLE epoch in minutes
        tsince = (et_vec(i) - sat_epoch_et) / 60;

        % Use the sgp4 function to compute the satellite state at the current time
        % `rteme` and `vteme` are the position and velocity in the TEME frame
        [~, rteme, vteme] = sgp4(satrec, tsince);

        % Convert the current ephemeris time to Julian centuries since 2000-01-01
        ttt = cspice_unitim(et_vec(i), 'ET', 'TDT') / cspice_jyear() / 100;

        % Convert the position and velocity from TEME to ECI frame using nutation corrections
        % `teme2eci` uses the conversion parameters and the Julian centuries from 2000-01-01
        [reci(:, i), veci(:, i)] = teme2eci(rteme, vteme, [0.0; 0.0; 0.0], ttt, ddpsi, ddeps);
    end
end

function [Y_pred, et_vec] = pointing_antenna(station, et_vec, reci_vec)
% This function calculates the predicted pointing of the antenna, converting
% satellite position vectors to the local topocentric coordinate system and 
% obtaining the range, azimuth, and elevation angles for the antenna's orientation.
%
% Inputs:
%   station   - Struct containing station parameters, including:
%               station.name: Name of the station
%   et_vec    - Vector of ephemeris times (in seconds past J2000)
%   reci_vec  - Matrix of satellite position vectors in the ECI frame [x; y; z] (in km),
%               each column corresponds to a time step
%
% Outputs:
%   Y_pred    - Predicted pointing parameters, including:
%               Y_pred(1, :) - Azimuth angles (in degrees)
%               Y_pred(2, :) - Elevation angles (in degrees)
%               Y_pred(3, :) - Range (in km)
%   et_vec    - The input vector of ephemeris times (unchanged)

    % Define the reference frames
    frame1 = 'J2000';  % ECI frame
    station_name = station.name;  % Get station name from input structure
    frame2 = strcat(station_name, '_TOPO');  % Topocentric frame for the station

    % Get the position of the station in the ECI frame at each time step
    r_station = cspice_spkpos(station_name, et_vec, frame1, 'NONE', 'Earth');

    % Calculate the rotation matrix from the ECI frame to the station's topocentric frame
    R = cspice_pxform(frame1, frame2, et_vec);

    % Calculate the relative position vector of the satellite with respect to the station
    r_rel = reci_vec - r_station;

    % Initialize the topocentric position vectors (in the station's frame)
    r_topo = nan(size(r_rel));  % Allocate space for the topocentric position vectors

    % Loop over each time step to convert from ECI to topocentric coordinates
    for k = 1:length(et_vec)
        % Apply the rotation matrix to convert each relative position vector to the topocentric frame
        r_topo(:, k) = R(:, :, k) * r_rel(:, k);
    end

    % Convert the topocentric position vectors to range, azimuth, and elevation
    [range, Az, El] = cspice_reclat(r_topo);

    % Convert azimuth and elevation from radians to degrees
    Az = Az * cspice_dpr;  % Azimuth angle in degrees
    El = El * cspice_dpr;  % Elevation angle in degrees

    % Return the predicted pointing parameters: azimuth, elevation, and range
    Y_pred = [Az; El; range];
end

function [Y_pred, et] = visibility(Y_pred, et, station)
% This function filters the predicted pointing data based on he minimum elevation angle 
% for visibility. It retains only the data points where the elevation angle is greater than
% the minimum elevation angle of the station.
%
% Inputs:
%   Y_pred    - Predicted pointing parameters, including:
%               Y_pred(1, :) - Azimuth angles (in degrees)
%               Y_pred(2, :) - Elevation angles (in degrees)
%               Y_pred(3, :) - Range (in km)
%   et        - Vector of ephemeris times (in seconds past J2000)
%   station   - Struct containing station parameters
%
% Outputs:
%   Y_pred    - Filtered predicted pointing parameters after visibility check:
%               Y_pred(1, :) - Azimuth angles (in degrees)
%               Y_pred(2, :) - Elevation angles (in degrees)
%               Y_pred(3, :) - Range (in km)
%   et        - Filtered ephemeris times corresponding to the visible times

    % Retrieve the minimum elevation angle from the station parameters (in degrees)
    El_min = station.El_min;

    % Extract the elevation angles from the predicted pointing parameters
    El = Y_pred(2, :);

    % Create a boolean array indicating whether each elevation angle exceeds the minimum threshold
    bool = double(El > El_min);

    % Find the indices where the elevation angle is greater than the minimum threshold
    indeces = bool ~= 0;

    % Filter the ephemeris times and predicted pointing data based on the visibility condition
    et = et(indeces);         % Retain only the visible times
    Y_pred = Y_pred(:, indeces);  % Retain only the corresponding visibility data
end

function Y = simulate_measurements(Y, station)
% This function simulates measurement noise on the predicted pointing parameters (Azimuth, Elevation, Range).
% The function adds random noise, modeled as Gaussian noise, to the predicted values based on the station's
% measurement noise covariance matrix.
%
% Inputs:
%   Y        - Matrix of predicted pointing parameters, including:
%              Y(1, :) - Azimuth angles (in degrees)
%              Y(2, :) - Elevation angles (in degrees)
%              Y(3, :) - Range (in km)
%   station  - Struct containing station parameters
%
% Outputs:
%   Y        - Simulated measurements, with added Gaussian noise

    % Retrieve the measurement noise covariance matrix from the station structure
    Wm = station.Wm;

    % Compute the measurement noise covariance matrix for each measurement by inverting Wm and squaring it
    R = (Wm \ eye(size(Wm))) .^ 2;

    % Convert the predicted values (Y) into a row vector (mu) for use in the multivariate normal distribution
    mu = Y.';

    % Generate simulated measurements by adding Gaussian noise to the predicted values
    Y = mvnrnd(mu, R).';
end

function residuals = costfunction(x, stations, J2)
% This function computes the residuals between the simulated and measured pointing data (Azimuth, Elevation, Range).
% It integrates the satellite's trajectory and compares the predicted pointing (from the propagation) 
% with the actual measurements (with some noise). The cost function residuals are weighted by the 
% measurement noise covariance matrix for optimization purposes.
%
% Inputs:
%   x         - Initial state vector [position; velocity] of the satellite
%   stations  - Array of station structs, each containing:
%               station.Wm: Measurement noise covariance matrix
%               station.Y: Measured pointing data (Azimuth, Elevation, Range)
%               station.tspan: Time span for measurements
%   J2        - Boolean flag indicating whether to include J2 perturbation in the orbital propagation
%
% Outputs:
%   residuals - Vector of residuals for the cost function, which will be minimized during optimization.
%               Residuals correspond to the difference between simulated and measured data, weighted by the measurement covariance.

    % Retrieve the gravitational parameter for Earth
    mu = cspice_bodvrd('Earth', 'GM', 1);

    % Set ODE solver options for high precision
    options = odeset('Reltol', 1.e-13, 'Abstol', 1.e-20);

    % Define the ODE function for propagation depending on whether J2 perturbation is included
    if J2
        % Earth’s equatorial radius (mean value for J2 perturbation)
        Re = cspice_bodvrd('Earth', 'RADII', 3);
        Re = (Re(3) + Re(2)) / 2;
        % Include J2 perturbation in the ODE function
        odefun = @(t, x) two_body_J2_rhs(t, x, mu, Re);
    else
        % Without J2 perturbation, use simple two-body dynamics
        odefun = @(t, x) two_body_rhs(t, x, mu);
    end

    % Initialize the number of residuals (sum of residuals for all stations)
    n_stations = length(stations);
    n_res = 0;
    for kk = 1:n_stations
        n_res = n_res + size(stations(kk).Y, 2);  % Add the number of measurements for each station
    end

    % Allocate space for the residuals matrix
    residuals = nan(3, n_res);
    nprec = 1;

    % Loop over each station
    for kk = 1:n_stations

        % Extract station-specific information
        station = stations(kk);
        Wm = station.Wm;  % Measurement noise covariance matrix
        Y = station.Y;    % Measured pointing data (Azimuth, Elevation, Range)
        tspan = station.tspan;  % Time span of the measurements
        n = size(Y, 2);  % Number of measurements for this station
        res = zeros(size(Y));  % Initialize the residual vector for this station

        % Propagate the satellite's state using the specified ODE function
        [~, x_prop] = ode78(odefun, tspan, x, options);  % ODE integration
        reci_vec = x_prop(:, 1:3).';  % Extract the position vectors (ECI)

        % Compute the predicted antenna pointing from the satellite's position
        [y_pred, ~] = pointing_antenna(station, tspan, reci_vec);
        Y_pred = y_pred(:, 2:end);  % Predicted values (excluding the first time step)

        % Loop over each time step to compute the residuals
        for k = 1:length(tspan)-1
            % Compute the difference between predicted and measured values
            dAz = angdiff(deg2rad(Y_pred(1, k)), deg2rad(Y(1, k)));  % Azimuth difference
            dEl = angdiff(deg2rad(Y_pred(2, k)), deg2rad(Y(2, k)));  % Elevation difference
            dr  = Y_pred(3, k) - Y(3, k);  % Range difference

            % Weight the residuals by the measurement noise covariance matrix
            diff_meas_weighted = Wm * [dAz; dEl; dr];
            res(:, k) = diff_meas_weighted;  % Store the weighted residuals
        end

        % Store the residuals for this station
        residuals(:, nprec:nprec+n-1) = res;
        nprec = nprec + n;  % Update the pointer to the next position in the residuals array

    end

end

function kep = car2kep(xx, mu)
% This function converts a satellite's state vector from Cartesian coordinates
% (position and velocity) to orbital elements.
%
% Inputs:
%   xx  - State vector in Cartesian coordinates:
%         xx(1:3) - Position vector [x; y; z] (in km)
%         xx(4:6) - Velocity vector [vx; vy; vz] (in km/s)
%   mu  - Gravitational parameter of the central body (in km^3/s^2)
%
% Outputs:
%   kep - Orbital elements corresponding to the input state vector:
%         kep(1)  - Semi-major axis (a) (in km)
%         kep(2)  - Eccentricity (e)
%         kep(3)  - Inclination (i) (in radians)
%         kep(4)  - right ascension of the ascending node (Ω) (in radians)
%         kep(5)  - Argument of periapsis (ω) (in radians)
%         kep(6)  - True anomaly (θ) (in radians)

    % Extract the position and velocity vectors from the input state vector `xx`
    rr = xx(1:3);  % Position vector [x; y; z] (in km)
    vv = xx(4:6);  % Velocity vector [vx; vy; vz] (in km/s)

    % Ensure that the position and velocity vectors are column vectors
    if size(rr,1) == 1
        rr = rr.';  % Transpose position vector if it's a row vector
    end

    if size(vv,1) == 1
        vv = vv.';  % Transpose velocity vector if it's a row vector
    end

    % Compute the orbital parameters
    r = norm(rr);  % Magnitude of the position vector (in km)
    v = norm(vv);  % Magnitude of the velocity vector (in km/s)
    
    % Define the z-axis unit vector
    k = [0 0 1]';

    % Semi-major axis (a)
    a = -mu / (v^2 - 2*mu/r);

    % Calculate the angular momentum vector `hh` and its magnitude `h`
    hh = cross(rr, vv);  % Angular momentum vector (in km^2/s)
    h = norm(hh);  % Magnitude of the angular momentum

    % Calculate the eccentricity vector `ee` and its magnitude `e`
    ee = cross(vv, hh) / mu - rr / r;  % Eccentricity vector
    e = norm(ee);  % Eccentricity (dimensionless)

    % If eccentricity is zero (circular orbit), set the eccentricity vector to [1; 0; 0]
    if e == 0
        ee = [1; 0; 0];
    end

    % Inclination angle (i) from the angular momentum vector
    i = acos(hh(3) / h);  % Inclination (in radians)

    % Calculate the unit vector for the ascending node (N)
    N = cross(k, hh) / norm(cross(hh, k));  % Node vector

    % If inclination is zero (equatorial orbit), set N to [1; 0; 0]
    if i == 0
        i = eps;  % Set small inclination to avoid division by zero
        N = [1; 0; 0];
    end

    % Longitude of the ascending node (Ω) from the node vector
    Om = acos(N(1));  % Ascending node longitude (in radians)
    if N(2) < 0  % Adjust for the correct quadrant
        Om = 2 * pi - Om;
    end

    % Argument of periapsis (ω) from the eccentricity vector
    om = acos(dot(N, ee) / e);  % Argument of periapsis (in radians)
    if ee(3) < 0  % Adjust for the correct quadrant
        om = 2 * pi - om;
    end

    % True anomaly (θ) from the position vector and eccentricity vector
    theta = acos(dot(rr, ee) / (r * e));  % True anomaly (in radians)
    if dot(rr, vv) < 0  % If the dot product is negative, the satellite is past periapsis
        theta = 2 * pi - theta;  % Adjust the true anomaly
    end

    % Return the orbital elements as a column vector
    kep = [a; e; i; Om; om; theta];  % Orbital elements [a, e, i, Om, om, theta]
end

function J = jacobian_car2kep(xx, mu)
% This function computes the Jacobian matrix of the car2kep function, which
% converts Cartesian state vectors (position and velocity) to orbital elements.
%
% Inputs:
%   xx  - State vector in Cartesian coordinates:
%         xx(1:3) - Position vector [x; y; z] (in km)
%         xx(4:6) - Velocity vector [vx; vy; vz] (in km/s)
%   mu  - Gravitational parameter of the central body (in km^3/s^2)
%
% Outputs:
%   J   - Jacobian matrix (n x n) where n is the number of elements in `xx`.
%         Each column contains the partial derivatives of the orbital elements
%         with respect to the corresponding component of the Cartesian state vector.

    % Get the number of elements in the input state vector `xx`
    n = length(xx);
    
    % Initialize the Jacobian matrix with NaN values
    J = nan(n);

    % Loop through each component of the input state vector `xx`
    for k = 1:n
        % Perturb the k-th component of `xx` by a small amount `ek`
        pert = zeros(size(xx));  % Initialize perturbation vector
        ek = sqrt(eps) * max(1, xx(k));  % Perturbation magnitude (small)

        % Apply the perturbation to create perturbed state vectors
        pert(k) = ek;  % Set the k-th component to the perturbation magnitude
        xx_plus  = xx + pert;  % State vector with positive perturbation
        xx_minus = xx - pert;  % State vector with negative perturbation

        % Convert the perturbed Cartesian vectors to Keplerian elements
        K_plus  = car2kep(xx_plus, mu);  % Keplerian elements for perturbed state
        K_minus = car2kep(xx_minus, mu); % Keplerian elements for perturbed state

        % Compute the finite difference for the partial derivative (central difference)
        J(:, k) = (K_plus - K_minus) / (2 * ek);  % Central difference approximation
    end
end

function sigma_vec = linear_mapping(xx_car, P_car, mu)
% This function performs a linear transformation of the state vector covariance
% from Cartesian coordinates to Keplerian orbital elements. It computes the 
% standard deviations (sigmas) for the Keplerian elements using the Jacobian
% matrix of the transformation.
%
% Inputs:
%   xx_car - State vector in Cartesian coordinates [position; velocity], 
%            where xx_car(1:3) is the position vector and xx_car(4:6) is the 
%            velocity vector (in km and km/s).
%   P_car  - Covariance matrix in Cartesian coordinates (6x6 matrix).
%            This matrix describes the uncertainties (covariances) in the 
%            Cartesian state vector components.
%   mu     - Gravitational parameter of the central body (in km^3/s^2).
%
% Outputs:
%   sigma_vec - A vector of standard deviations for the orbital elements
%               (semi-major axis, eccentricity, inclination, tight ascension of 
%               ascending node, argument of periapsis, and true anomaly).
%               The output is a 1x6 vector, where each element corresponds
%               to the standard deviation of a specific orbital element.

    % Compute the Jacobian matrix of the transformation from Cartesian to Keplerian
    % orbital elements using the car2kep function.
    J = jacobian_car2kep(xx_car, mu);  % Jacobian matrix (6x6)

    % Perform the linear transformation of the covariance matrix from Cartesian
    % coordinates to Keplerian orbital elements.
    P_kep = J * P_car * J.';  % Covariance matrix in Keplerian coordinates

    % Extract the standard deviations (square roots of the diagonal elements) 
    % from the covariance matrix of orbital elements.
    sigma_vec = sqrt(diag(P_kep));  % Standard deviations (1x6 vector)
end

function [cost, sigma_a, sigma_i] = trade_off_analysis(stations, x0_guess)
% This function computes the standard deviations of orbital 
% elements (semi-major axis and inclination) and the total cost.
%
% Inputs:
%   stations   - A structure array containing the details of the observation 
%                stations, including the cost function for each station.
%   x0_guess   - An initial guess for the orbital parameters (6x1 vector).
%                This guess is used as the starting point for the optimization.
%
% Outputs:
%   cost       - The total cost, summing up the individual costs 
%                from all stations.
%   sigma_a    - The standard deviation of the semi-major axis (in km).
%   sigma_i    - The standard deviation of the inclination (in degrees).
    
    % Initialize cost to zero
    cost = 0;
    
    % Get the gravitational parameter of Earth (in km^3/s^2)
    mu = cspice_bodvrd('Earth', 'GM', 1);
    
    % Define the cost function that will be minimized during optimization
    fun = @(x) costfunction(x, stations, true);  % Cost function with J2 effect included
    
    % Set options
    options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'none');
    
    % Solve the navigation problem
    [x, resnorm, residual, ~, ~, ~, Jac] = lsqnonlin(fun, x0_guess, [], [], options);
    
    % Compute the matrix B
    B = (Jac.' * Jac) \ eye(size(Jac, 2));
    
    % Compute the covariance matrix P
    P = resnorm / (length(residual) - length(x)) * B;
    
    % Compute the standard deviations of the orbital elements
    sigma_vec = linear_mapping(x, P, mu);
    
    % Extract the standard deviations for the semi-major axis and inclination
    sigma_a = sigma_vec(1);  % Standard deviation of semi-major axis (in km)
    sigma_i = rad2deg(sigma_vec(3));  % Standard deviation of inclination (in degrees)
    
    % Sum up the costs for each station to get the total cost
    for k = 1:length(stations)
        station = stations(k);
        cost = cost + station.cost;  % Add the individual station cost
    end
end

function [ellipsoid_x, ellipsoid_y, ellipsoid_z] = confidence_ellipsoid(P, mu, n_sigma)
    % Confidence ellipsoid for a 3D Gaussian distribution
    % P: Covariance matrix (3x3)
    % mu: Mean vector (3x1)
    % n_sigma: Confidence level (number of standard deviations)
    [xe, ye, ze] = sphere(100);
    points = [xe(:), ye(:), ze(:)]';

    % Perform Singular Value Decomposition of P
    [R, D] = svd(P);  % R is the rotation matrix, D is the diagonal matrix with eigenvalues

    % Extract the square root of eigenvalues from D (standard deviations along each axis)
    d = sqrt(diag(D));  % d will be a 3x1 vector of standard deviations

    % Transformation: R * d to stretch the circle and rotate
    tr_pts = n_sigma * R * diag(d) * points;
    X = reshape(tr_pts(1,:), size(xe));
    Y = reshape(tr_pts(2,:), size(ye));
    Z = reshape(tr_pts(3,:), size(ze));

    % x and y coordinates of the ellipse
    ellipsoid_x = mu(1) + X;  % x-coordinates of the ellipse
    ellipsoid_y = mu(2) + Y;  % y-coordinates of the ellipse
    ellipsoid_z = mu(3) + Z;  % y-coordinates of the ellipse
end
