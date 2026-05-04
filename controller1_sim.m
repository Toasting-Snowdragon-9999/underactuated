% CONTROLLER1_SIM  Closed-loop simulation of the cart with double pendulum
% under the local LQR stabilizer in control1.m.
%
% Scenario: start with both poles slightly off vertical and at rest, then
% let control1 drive the system back to the upright equilibrium.
%
% This file exercises ONLY control1 (no swing-up). For the energy-shaping
% swing-up controller, see controller2_sim.m.

clear; close all; clc;

%% ------------------------------------------------------------------------
%  Parameters  (lab-grade cart-pendulum: realistic light viscous friction)
%  ------------------------------------------------------------------------
params.M  = 1.0;     % cart mass               [kg]
params.m1 = 0.20;    % pole-1 bob mass         [kg]
params.m2 = 0.15;    % pole-2 bob mass         [kg]
params.l1 = 0.50;    % pole-1 rod length       [m]
params.l2 = 0.40;    % pole-2 rod length       [m]
params.g  = 9.81;    % gravity                 [m/s^2]

params.b_cart  = 0.05;     % linear-bearing cart friction  [N s / m]
params.b_pend1 = 0.001;    % pole-1 joint friction         [N m s / rad]
params.b_pend2 = 0.001;    % pole-2 joint friction         [N m s / rad]

opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

%% ------------------------------------------------------------------------
%  Run -- LQR balance from a small perturbation around upright
%  ------------------------------------------------------------------------
fprintf('--- LQR balance (control1) ---\n');
x0    = [0; 0.15; -0.15; 0; 0; 0];   % small angle perturbation, at rest
tspan = [0 6];

u_fn = @(t, x) control1(t, x, params);

[t, X] = ode45(@(t, x) underactuated_model(t, x, u_fn, params), ...
               tspan, x0, opts);

U = arrayfun(@(k) control1(t(k), X(k,:).', params), 1:numel(t));

fprintf('  final |theta1| = %.3e rad,  |theta2| = %.3e rad\n', ...
        abs(X(end,2)), abs(X(end,3)));
fprintf('  peak  |u|      = %.2f N\n\n', max(abs(U)));

%% ------------------------------------------------------------------------
%  Plots and animation
%  ------------------------------------------------------------------------
figure('Name', 'control1 -- LQR balance');
subplot(4,1,1); plot(t, X(:,1), 'LineWidth', 1.1); grid on;
    ylabel('s [m]'); title('control1 -- LQR balance');
subplot(4,1,2); plot(t, X(:,2:3), 'LineWidth', 1.1); grid on;
    ylabel('\theta [rad]'); legend('\theta_1', '\theta_2');
subplot(4,1,3); plot(t, X(:,4:6), 'LineWidth', 1.1); grid on;
    ylabel('velocity'); legend('sdot', '\theta_1 dot', '\theta_2 dot');
subplot(4,1,4); plot(t, U, 'LineWidth', 1.1); grid on;
    xlabel('t [s]'); ylabel('u [N]');

animate_cartpole(t, X, params);
