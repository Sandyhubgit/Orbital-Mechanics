% Orbital Simulation with Radial Distance, Radial Velocity, Tangential Velocity, and Orbit Shape Classification
clc;
clear;
close all;

% Input for current radial distance, radial velocity, and tangential velocity
r_current = input('Enter the current radial distance from the center of Earth in km: ') * 1e3; % Convert km to meters
v_radial = input('Enter the radial velocity (velocity towards/away from Earth) in km/s: ')*1e3;
v_tangential = input('Enter the tangential velocity (velocity perpendicular to the radial direction) in km/s: ')*1e3;

% Earth's parameters
mu_earth = 3.986e14; % Standard gravitational parameter of Earth (m^3/s^2)
radius_earth = 6371e3; % Radius of the Earth in meters

% Orbital energy equation to calculate semi-major axis (a)
v_total = sqrt(v_radial^2 + v_tangential^2); % Total velocity (magnitude)
specific_orbital_energy = (v_total^2 / 2) - (mu_earth / r_current);
semi_major_axis = -mu_earth / (2 * specific_orbital_energy);

% Eccentricity calculation
eccentricity = sqrt(1 + (2 * specific_orbital_energy * (r_current^2 * v_tangential^2)) / mu_earth^2);

% Calculate perigee and apogee distances for elliptical orbits
if eccentricity < 1
    perigee = semi_major_axis * (1 - eccentricity);  % Closest point to Earth
    apogee = semi_major_axis * (1 + eccentricity);   % Farthest point from Earth
else
    perigee = NaN; % Perigee doesn't exist for hyperbolic orbits
    apogee = NaN;  % Apogee doesn't exist for hyperbolic orbits
end

% Velocities at perigee and apogee (m/s)
if eccentricity < 1
    v_perigee = sqrt(mu_earth * (2 / perigee - 1 / semi_major_axis));
    v_apogee = sqrt(mu_earth * (2 / apogee - 1 / semi_major_axis));
else
    v_perigee = NaN;
    v_apogee = NaN;
end

% Orbital period for elliptical orbits
if eccentricity < 1
    orbital_period = 2 * pi * sqrt(semi_major_axis^3 / mu_earth);
else
    orbital_period = NaN; % No orbital period for hyperbolic orbits
end

% Display calculated parameters
fprintf('Semi-major axis: %.2f km\n', semi_major_axis / 1e3);
fprintf('Eccentricity: %.4f\n', eccentricity);

% Orbit shape classification based on eccentricity
if eccentricity == 0
    fprintf('Orbit shape: Circle\n');
elseif eccentricity > 0 && eccentricity < 1
    fprintf('Orbit shape: Ellipse\n');
elseif eccentricity == 1
    fprintf('Orbit shape: Parabola\n');
elseif eccentricity > 1
    fprintf('Orbit shape: Hyperbola\n');
end

if eccentricity < 1
    fprintf('Perigee: %.2f km\n', perigee / 1e3);
    fprintf('Apogee: %.2f km\n', apogee / 1e3);
    fprintf('Maximum velocity at perigee: %.2f m/s\n', v_perigee);
    fprintf('Minimum velocity at apogee: %.2f m/s\n', v_apogee);
    fprintf('Orbital period: %.2f seconds\n', orbital_period);
else
    fprintf('Perigee, apogee, and orbital period do not apply to hyperbolic orbits.\n');
end

% Simulation and plotting for elliptical, parabolic, and hyperbolic orbits
theta_steps = 1000; % Number of steps for the angular component

if eccentricity < 1  % Elliptical orbit
    theta = linspace(0, 2 * pi, theta_steps); % True anomaly (angle around the orbit)
    
    % Position in the orbital plane (polar coordinates to Cartesian)
    r = (semi_major_axis * (1 - eccentricity^2)) ./ (1 + eccentricity * cos(theta)); % Orbit equation
    
    x = r .* cos(theta); % x-coordinates of orbit
    y = r .* sin(theta); % y-coordinates of orbit
    
    % Plot the orbit
    figure;
    plot(x / 1e3, y / 1e3, 'b-', 'LineWidth', 1.5); % Plot in km
    hold on;
    plot(0, 0, 'ro', 'MarkerFaceColor', 'r'); % Earth's center
    title('Satellite Orbit (Elliptical)');
    xlabel('x (km)');
    ylabel('y (km)');
    axis equal;
    grid on;
    
elseif eccentricity == 1  % Parabolic orbit
    theta = linspace(-pi/2, pi/2, theta_steps); % Restricted angle for open trajectory
    
    % Position in the orbital plane (polar coordinates to Cartesian)
    r = (semi_major_axis * (1 - eccentricity^2)) ./ (1 + eccentricity * cos(theta)); % Orbit equation
    
    x = r .* cos(theta); % x-coordinates of orbit
    y = r .* sin(theta); % y-coordinates of orbit
    
    % Plot the orbit
    figure;
    plot(x / 1e3, y / 1e3, 'g-', 'LineWidth', 1.5); % Plot in km
    hold on;
    plot(0, 0, 'ro', 'MarkerFaceColor', 'r'); % Earth's center
    title('Satellite Orbit (Parabolic)');
    xlabel('x (km)');
    ylabel('y (km)');
    axis equal;
    grid on;
    
elseif eccentricity > 1  % Hyperbolic orbit
    theta = linspace(-pi/2, pi/2, theta_steps); % Restricted angle for open trajectory
    
    % Position in the orbital plane (polar coordinates to Cartesian)
    r = (semi_major_axis * (1 - eccentricity^2)) ./ (1 + eccentricity * cos(theta)); % Orbit equation
    
    x = r .* cos(theta); % x-coordinates of orbit
    y = r .* sin(theta); % y-coordinates of orbit
    
    % Plot the orbit
    figure;
    plot(x / 1e3, y / 1e3, 'r-', 'LineWidth', 1.5); % Plot in km
    hold on;
    plot(0, 0, 'ro', 'MarkerFaceColor', 'r'); % Earth's center
    title('Satellite Orbit (Hyperbolic)');
    xlabel('x (km)');
    ylabel('y (km)');
    axis equal;
    grid on;
end
