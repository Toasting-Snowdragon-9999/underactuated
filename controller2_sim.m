% CONTROLLER2_SIM  Closed-loop simulation of the cart with double pendulum
% under the energy-shaping swing-up controller in control2.m.
%
% Scenario: start near the hanging equilibrium and let control2 pump
% mechanical energy until the chain reaches the upright neighborhood, at
% which point control2 internally hands off to the LQR (control1).
%
% This file exercises control2 (which itself calls control1 once near
% upright). For a pure LQR-balance run that never leaves the upright
% basin, see controller1_sim.m.

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
%  Run -- swing-up from near hanging
%  ------------------------------------------------------------------------
fprintf('--- Swing-up + LQR (control2) ---\n');
x0    = [0; pi - 0.1; pi + 0.1; 0; 0; 0];   % near hanging, slight kick
tspan = [0 30];

u_fn = @(t, x) control2(t, x, params);

[t, X] = ode45(@(t, x) underactuated_model(t, x, u_fn, params), ...
               tspan, x0, opts);

U = arrayfun(@(k) control2(t(k), X(k,:).', params), 1:numel(t));

% Track total mechanical energy and the upright reference E_des.
E     = arrayfun(@(k) total_energy(X(k,:).', params), 1:numel(t));
E_des = (params.m1 + params.m2)*params.g*params.l1 ...
       + params.m2*params.g*params.l2;

fprintf('  final |theta1| = %.3e rad,  |theta2| = %.3e rad\n', ...
        abs(wrap_pi(X(end,2))), abs(wrap_pi(X(end,3))));
fprintf('  E_des = %.3f J,  final E = %.3f J\n', E_des, E(end));
fprintf('  peak  |u|      = %.2f N\n\n', max(abs(U)));

%% ------------------------------------------------------------------------
%  Plots and animation
%  ------------------------------------------------------------------------
figure('Name', 'control2 -- swing-up + LQR');
subplot(5,1,1); plot(t, X(:,1), 'LineWidth', 1.1); grid on;
    ylabel('s [m]'); title('control2 -- swing-up + LQR');
subplot(5,1,2); plot(t, X(:,2:3), 'LineWidth', 1.1); grid on;
    ylabel('\theta [rad]'); legend('\theta_1', '\theta_2');
subplot(5,1,3); plot(t, X(:,4:6), 'LineWidth', 1.1); grid on;
    ylabel('velocity'); legend('sdot', '\theta_1 dot', '\theta_2 dot');
subplot(5,1,4); plot(t, U, 'LineWidth', 1.1); grid on;
    ylabel('u [N]');
subplot(5,1,5); plot(t, E, 'LineWidth', 1.1); grid on; hold on;
    yline(E_des, 'r--', 'E_{des}');
    xlabel('t [s]'); ylabel('E [J]');

animate_cartpole(t, X, params);


% =========================================================================
%   Local helper functions
% =========================================================================
function E = total_energy(x, p)
    sd  = x(4);
    th1 = x(2);  th2 = x(3);
    th1d = x(5); th2d = x(6);

    x1d =  sd + p.l1*cos(th1)*th1d;
    y1d =     - p.l1*sin(th1)*th1d;
    x2d =  sd + p.l1*cos(th1)*th1d + p.l2*cos(th2)*th2d;
    y2d =     - p.l1*sin(th1)*th1d - p.l2*sin(th2)*th2d;

    T = 0.5*p.M*sd^2 ...
      + 0.5*p.m1*(x1d^2 + y1d^2) ...
      + 0.5*p.m2*(x2d^2 + y2d^2);
    V = (p.m1 + p.m2)*p.g*p.l1*cos(th1) + p.m2*p.g*p.l2*cos(th2);

    E = T + V;
end

% -------------------------------------------------------------------------
function a = wrap_pi(a)
    a = mod(a + pi, 2*pi) - pi;
end
