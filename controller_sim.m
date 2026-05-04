% CONTROLLER_SIM  Closed-loop simulations of the cart with double pendulum.
%
% Two scenarios:
%   (S1) LQR balance     -- start near the upright equilibrium with a small
%                           perturbation and stabilize using control1.
%   (S2) Energy swing-up -- start near the hanging equilibrium and use
%                           control2, which pumps energy until the system
%                           is near upright and then hands off to LQR.
%
% Each scenario produces a state/input plot and a real-time animation.

clear; close all; clc;

%% ------------------------------------------------------------------------
%  Parameters
%  ------------------------------------------------------------------------
params.M  = 1.0;     % cart mass               [kg]
params.m1 = 0.20;    % pole-1 bob mass         [kg]
params.m2 = 0.15;    % pole-2 bob mass         [kg]
params.l1 = 0.50;    % pole-1 rod length       [m]
params.l2 = 0.40;    % pole-2 rod length       [m]
params.g  = 9.81;    % gravity                 [m/s^2]

% Mild friction here so the swing-up controller has a chance; the model
% demo (simulation.m) intentionally uses much larger ("vicious") values.
params.b_cart  = 0.10;
params.b_pend1 = 0.005;
params.b_pend2 = 0.005;

opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

%% ------------------------------------------------------------------------
%  (S1) LQR balance from a small perturbation around upright
%  ------------------------------------------------------------------------
fprintf('--- (S1) LQR balance ---\n');
x0_lqr    = [0; 0.15; -0.15; 0; 0; 0];   % small angles, at rest
tspan_lqr = [0 6];

u_lqr = @(t, x) control1(t, x, params);
[t_lqr, X_lqr] = ode45(@(t, x) underactuated_model(t, x, u_lqr, params), ...
                       tspan_lqr, x0_lqr, opts);

U_lqr = arrayfun(@(k) control1(t_lqr(k), X_lqr(k,:).', params), ...
                 1:numel(t_lqr));

fprintf('  final |theta1| = %.3e rad,  |theta2| = %.3e rad\n', ...
        abs(X_lqr(end,2)), abs(X_lqr(end,3)));
fprintf('  peak |u|       = %.2f N\n\n', max(abs(U_lqr)));

plot_run('(S1) LQR balance', t_lqr, X_lqr, U_lqr);

%% ------------------------------------------------------------------------
%  (S2) Energy-based swing-up + hand-off to LQR
%  ------------------------------------------------------------------------
fprintf('--- (S2) Swing-up + LQR ---\n');
x0_swing    = [0; pi - 0.1; pi + 0.1; 0; 0; 0];   % near hanging
tspan_swing = [0 30];

u_swing = @(t, x) control2(t, x, params);
[t_sw, X_sw] = ode45(@(t, x) underactuated_model(t, x, u_swing, params), ...
                     tspan_swing, x0_swing, opts);

U_sw = arrayfun(@(k) control2(t_sw(k), X_sw(k,:).', params), ...
                1:numel(t_sw));

fprintf('  final |theta1| = %.3e rad,  |theta2| = %.3e rad\n', ...
        abs(wrap_pi(X_sw(end,2))), abs(wrap_pi(X_sw(end,3))));
fprintf('  peak |u|       = %.2f N\n\n', max(abs(U_sw)));

plot_run('(S2) Swing-up + LQR', t_sw, X_sw, U_sw);

%% ------------------------------------------------------------------------
%  Animations -- real-time playback (one window per scenario)
%  ------------------------------------------------------------------------
animate_cartpole(t_lqr, X_lqr, params);
animate_cartpole(t_sw,  X_sw,  params);


% =========================================================================
%   Local helper functions
% =========================================================================
function plot_run(name, t, X, U)
    figure('Name', name);
    subplot(4,1,1); plot(t, X(:,1), 'LineWidth', 1.1); grid on;
        ylabel('s [m]'); title(name);
    subplot(4,1,2); plot(t, X(:,2:3), 'LineWidth', 1.1); grid on;
        ylabel('\theta [rad]'); legend('\theta_1', '\theta_2');
    subplot(4,1,3); plot(t, X(:,4:6), 'LineWidth', 1.1); grid on;
        ylabel('velocity'); legend('sdot', '\theta_1 dot', '\theta_2 dot');
    subplot(4,1,4); plot(t, U, 'LineWidth', 1.1); grid on;
        xlabel('t [s]'); ylabel('u [N]');
end

% -------------------------------------------------------------------------
function a = wrap_pi(a)
    a = mod(a + pi, 2*pi) - pi;
end
