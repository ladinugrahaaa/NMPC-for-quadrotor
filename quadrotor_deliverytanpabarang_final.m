clear;
clc;

%% MODEL PARAMETERS
g = 9.81; % Gravity (m/s^2)
m = 1.343; % Mass (kg)
Ix = 0.01864; % Moment of inertia around x-axis
Iy = 0.037; % Moment of inertia around y-axis
Iz = 0.00554; % Moment of inertia around z-axis
kt = 0.1939;
Jr = 0.007588; 
l = 0.225;
Omega = 0; 

%% TIME
dt = 0.01; t = 0:dt:120; steps = length(t);

%% STATE INITIALIZATION
X = zeros(12, steps); % [x, xdot, y, ydot, z, zdot, phi, phidot, theta, thetadot, psi, psidot]
X(:, 1) = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]; % Initial state

%% INPUT INITIALIZATION
U = zeros(4, steps);
U_min = [-100; -100; -100; -100];
U_max = [100; 100; 100; 100];

%% MPC PARAMETERS
N = 50;
Q = diag([10000, 0, 10000, 0, 10000, 0, 1000, 0, 1000, 0, 1000, 0]);
R = diag([0.1, 0.1, 0.1, 0.1]);

%% DESIRED TRAJECTORY
length = 20; % Panjang lintasan (m)
width = 20; % Lebar lintasan (m)
z_height = 10; % Ketinggian lintasan (m)
side_time = 20; % Waktu untuk menempuh satu sisi lintasan (s)
total_time = 4 * side_time; % Total waktu lintasan
Xsp = zeros(6, steps); % Trajectory reference
Xspdot = zeros(3, steps); % Reference velocity

% Titik pengantaran barang
delivery_point = [10; 10; 1.5];
delivery_duration = 15; % Durasi berhenti di titik pengantaran (detik)
delivery_time = side_time; % Waktu mencapai titik pengantaran
resume_time = delivery_time + delivery_duration; % Waktu untuk kembali melanjutkan lintasan

for k = 1:steps
    current_time = t(k);
    if k <= 7000
        if current_time <= delivery_time
            % Trajektori menuju titik pengantaran
            Xsp(1:3, k) = [length * current_time / delivery_time; 0; z_height];
            Xspdot(:, k) = [length / delivery_time; 0; 0];
            Xsp(6, k) = 0; % phi, theta, psi
        elseif current_time <= resume_time
            % Drone berhenti di titik pengantaran
            Xsp(1:3, k) = delivery_point;
            Xspdot(:, k) = [0; 0; 0];
            Xsp(4:6, k) = [0; 0; 0];
        else
            % Lanjutkan lintasan setelah pengantaran
            adjusted_time = current_time - resume_time;
            if adjusted_time <= side_time
                Xsp(1:3, k) = [length - length * adjusted_time / side_time; width; z_height];
                Xspdot(:, k) = [-length / side_time; 0; 0];
                Xsp(6, k) = 0; % psi
            elseif adjusted_time <= 2 * side_time
                Xsp(1:3, k) = [0; width - width * (adjusted_time - side_time) / side_time; z_height];
                Xspdot(:, k) = [0; -width / side_time; 0];
                Xsp(6, k) = 0; % psi
            elseif adjusted_time <= 3 * side_time
                Xsp(1:3, k) = [length * (adjusted_time - 2 * side_time) / side_time; 0; z_height];
                Xspdot(:, k) = [length / side_time; 0; 0];
                Xsp(6, k) = 0; % psi
            else
                Xsp(1:3, k) = [0; 0; z_height];
                Xspdot(:, k) = [0; 0; 0];
                Xsp(6, k) = 0; % psi
            end
        end
    end

    if k >= 7000
        Xsp(1:3, k) = X(1:3, 1); 
        Xspdot(:, k) = [0; 0; 0];
        Xsp(4:6, k) = [0; 0; 0];
        if norm(X(1:3, k) - Xsp(1:3, k)) < 0.1
            break;
        end
    end
end


%% LOOP OVER TIME STEPS
figure(1);
stop_flag = false;
U_opt_array = zeros(4 * N, steps - 1);
for k = 1:steps-1
    if stop_flag
        break;
    end

    u_guess = repmat(U(:, k), 1, N);
    options = optimset('Display', 'off', 'MaxIter', 100);

    % Optimize control inputs using fmincon
    U_opt = fmincon(@(u_seq) computeCostNonlinear(X(:, k), reshape(u_seq, 4, N), ...
                Xsp(:, min(k:k+N-1, end)), Q, R, N, dt, g, m, Ix, Iy, Iz, kt, Jr, Omega), ...
                u_guess(:), [], [], [], [], repmat(U_min, N, 1), repmat(U_max, N, 1), [], options);

    U_opt_array(:, k) = U_opt;

    U(:, k) = U_opt(1:4);

    Xdot = computeStateDerivative(X(:, k), U(:, k), g, m, Ix, Iy, Iz, kt, Jr, Omega);
    X(:, k+1) = X(:, k) + dt * Xdot;

    if k * dt >= total_time
        stop_flag = true;
    end

    % Plot animasi
    if mod(k, 10) == 0
        fprintf('Step %d: X = [x: %.2f, y: %.2f, z: %.2f, phi: %.2f°, theta: %.2f°, psi: %.2f°], m = %.3f kg\n', ...
        k, X(1, k), X(3, k), X(5, k), rad2deg(X(7, k)), rad2deg(X(9, k)), rad2deg(X(11, k)), m);
        
        clf;
        hold on;

        plot3(Xsp(1, 1:k), Xsp(2, 1:k), Xsp(3, 1:k), 'k--', 'LineWidth', 1.5); % Reference trajectory
        plot3(X(1, 1:k), X(3, 1:k), X(5, 1:k), 'b', 'LineWidth', 1.5); % Actual trajectory
        plot3(X(1, k), X(3, k), X(5, k), 'ro', 'MarkerFaceColor', 'r'); % Current position

        scatter3(delivery_point(1), delivery_point(2), delivery_point(3), 100, 'g', 'filled', 'o'); % Delivery point
        scatter3([0, length, length, 0], [0, 0, width, width], [z_height, z_height, z_height, z_height], ...
            50, 'm', 'filled', 's'); % Corner points
        scatter3(X(1, 1), X(3, 1), X(5, 1), 120, 'x', 'MarkerEdgeColor', 'black', 'LineWidth', 2); % Start point

        axis equal;
        title('3D Trajectory of the Drone with Delivery');
        xlabel('x [m]');
        ylabel('y [m]');
        zlabel('z [m]');
        legend('Reference Trajectory', 'Actual Trajectory', 'Drone Position', 'Delivery Point', 'Corner Points', 'Takeoff and Landing Point');
        grid on;
        view(3);
        drawnow;
    end
end


%% PLOT STATES AND CONTROL INPUTS
figure(2);
for i = 1:12
    subplot(3, 4, i);
    plot(t(1:k), X(i, 1:k), 'LineWidth', 1.5);
    title(['State ', num2str(i)]);
    xlabel('Time [s]');
    ylabel(['X', num2str(i)]);
    grid on;
end

figure(3);
for i = 1:4
    subplot(2, 2, i);
    plot(t(1:k), U(i, 1:k), 'LineWidth', 1.5);
    title(['Control Input U', num2str(i)]);
    xlabel('Time [s]');
    ylabel(['U', num2str(i)]);
    grid on;
end

%% CALCULATE AND PLOT RMSE
rmse = zeros(6, 1);

for i = 1:6

    rmse(i) = sqrt(mean((X(2*i-1, 1:k) - Xsp(i, 1:k)).^2));
end

% Plot RMSE
figure(4);
bar(rmse, 'FaceColor', [0.2, 0.6, 0.8], 'EdgeColor', 'k');
title('RMSE of States (Position and Orientation)');
xlabel('State Index (x, y, z, φ, θ, ψ)');
ylabel('RMSE');
xticklabels({'x', 'y', 'z', 'φ', 'θ', 'ψ'});
grid on;

%% COST FUNCTION FOR NONLINEAR SYSTEM
function cost = computeCostNonlinear(X, U_seq, Xsp, Q, R, N, dt, g, m, Ix, Iy, Iz, kt, Jr, Omega)
    cost = 0;
    X_pred = X;

    for i = 1:N
        U = U_seq(:, i);
        Xdot_pred = computeStateDerivative(X_pred, U, g, m, Ix, Iy, Iz, kt, Jr, Omega);
        X_pred = X_pred + dt * Xdot_pred;

        e_pos = X_pred([1, 3, 5]) - Xsp(1:3, min(i, size(Xsp, 2)));
        e_angles = X_pred([7, 9, 11]) - Xsp(4:6, min(i, size(Xsp, 2)));

        Q_pos = Q([1, 3, 5], [1, 3, 5]);
        Q_angles = Q([7, 9, 11], [7, 9, 11]);
        cost = cost + (e_pos' * Q_pos * e_pos + e_angles' * Q_angles * e_angles + U' * R * U) * dt;
    end
end


%% NONLINEAR STATE DERIVATIVE FUNCTION
function Xdot = computeStateDerivative(X, U, g, m, Ix, Iy, Iz, kt, Jr, Omega)
    x2 = X(2); y2 = X(4); z2 = X(6); phi = X(7); phidot = X(8);
    theta = X(9); thetadot = X(10); psi = X(11); psidot = X(12);

    c7 = cos(phi); s7 = sin(phi);
    c9 = cos(theta); s9 = sin(theta);
    c11 = cos(psi); s11 = sin(psi);

    Xdot = zeros(12, 1);
    Xdot(1) = x2;
    Xdot(2) = -U(1)/m * (c7 * s9 * c11 + s7 * s11) + (kt/m) * x2;
    Xdot(3) = y2;
    Xdot(4) = -U(1)/m * (c7 * s9 * s11 - s7 * c11) + (kt/m) * y2;
    Xdot(5) = z2;
    Xdot(6) = g - U(1)/m * (c7 * c9) + (kt/m) * z2;
    Xdot(7) = phidot;
    Xdot(8) = (0.225/Ix) * U(2) + (Iy - Iz)/Ix * thetadot * psidot - kt/Ix * phidot - Jr/Ix * thetadot;
    Xdot(9) = thetadot;
    Xdot(10) = (0.225/Iy) * U(3) + (Iz - Ix)/Iy * phidot * psidot - kt/Iy * thetadot - Jr/Iy * phidot;
    Xdot(11) = psidot;
    Xdot(12) = (0.225/Iz) * U(4) + (Ix - Iy)/Iz * phidot * thetadot - kt/Iz * psidot;
end
