% SIMULATION  Simulate and verify the cart-pole Euler-Lagrange model.
%
% Verification strategy:
%   (V1) Energy conservation. With u = 0 and no friction, total mechanical
%        energy E = T + V must be constant along trajectories. We integrate
%        the model and check the relative drift.
%   (V2) Equilibria. Initial conditions at the upright (theta = 0) and
%        hanging (theta = pi) equilibria, with zero velocity and u = 0,
%        must produce zero acceleration. We check this from the model
%        directly and from a short integration.
%   (V3) Linearization. Numerical Jacobians around the two equilibria are
%        compared against the closed-form linearizations from the
%        manipulator equations. Eigenvalues are reported and the
%        controllability rank is checked.
%   (V4) Animation. A visual sanity check of the swing.

clear; close all; clc;

%% ------------------------------------------------------------------------
%  Parameters
%  ------------------------------------------------------------------------
params.M = 1.0;      % cart mass            [kg]
params.m = 0.2;      % pendulum bob mass    [kg]
params.l = 0.5;      % pendulum rod length  [m]
params.g = 9.81;     % gravity              [m/s^2]
params.b_cart = 4.0;   % cart viscous friction   [N s / m]    -- vicious
params.b_pend = 0.25;  % pendulum viscous fric.  [N m s / rad] -- vicious

%% ------------------------------------------------------------------------
%  (V1) Energy conservation  --  free response from a perturbed hang
%  ------------------------------------------------------------------------
x0_swing = [0; pi - 0.4; 0; 0];   % small push from the hanging equilibrium
tspan    = [0 10];
opts     = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

u_zero = @(t, x) 0;

% V1 must run friction-free regardless of the global params.
params_nf = params;
params_nf.b_cart = 0;
params_nf.b_pend = 0;

[t1, X1] = ode45(@(t, x) underactuated_model(t, x, u_zero, params_nf), ...
                 tspan, x0_swing, opts);

E1   = arrayfun(@(k) total_energy(X1(k,:).', params), 1:size(X1,1));
relE = (E1 - E1(1)) / max(abs(E1(1)), eps);

fprintf('--- (V1) Energy conservation (u = 0, friction = 0) ---\n');
fprintf('  E(0)              = %+ .6f  J\n', E1(1));
fprintf('  max |E(t) - E(0)| = % .3e J\n', max(abs(E1 - E1(1))));
fprintf('  max relative drift= % .3e\n\n', max(abs(relE)));

figure('Name', '(V1) States, free response');
subplot(2,2,1); plot(t1, X1(:,1)); grid on;
    xlabel('t [s]'); ylabel('s [m]');         title('cart position');
subplot(2,2,2); plot(t1, X1(:,2)); grid on;
    xlabel('t [s]'); ylabel('theta [rad]');   title('pendulum angle');
subplot(2,2,3); plot(t1, X1(:,3)); grid on;
    xlabel('t [s]'); ylabel('sdot [m/s]');    title('cart velocity');
subplot(2,2,4); plot(t1, X1(:,4)); grid on;
    xlabel('t [s]'); ylabel('thetadot [rad/s]'); title('pendulum rate');

figure('Name', '(V1) Energy drift');
plot(t1, E1 - E1(1), 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('E(t) - E(0)  [J]');
title(sprintf('Total energy drift  (max rel. err = %.2e)', max(abs(relE))));

%% ------------------------------------------------------------------------
%  (V2) Equilibria
%  ------------------------------------------------------------------------
fprintf('--- (V2) Equilibria ---\n');
x_up   = [0; 0;  0; 0];
x_down = [0; pi; 0; 0];

f_up   = underactuated_model(0, x_up,   0, params);
f_down = underactuated_model(0, x_down, 0, params);

fprintf('  f(x_upright,   u=0) = [% .2e % .2e % .2e % .2e]\n', f_up);
fprintf('  f(x_hanging,   u=0) = [% .2e % .2e % .2e % .2e]\n\n', f_down);

%% ------------------------------------------------------------------------
%  (V3) Linearization around upright + hanging equilibria
%  ------------------------------------------------------------------------
fprintf('--- (V3) Linearization ---\n');

f_xu = @(x, u) underactuated_model(0, x, u, params);

[A_up,   B_up  ] = numerical_linearization(f_xu, x_up,   0);
[A_down, B_down] = numerical_linearization(f_xu, x_down, 0);

% Closed-form linearization from M qddot + G = B u about an equilibrium
% q* with qdot* = 0:  qddot ~= M(q*)^{-1} ( B u - dG/dq |_{q*} dq ).
[A_up_cf,   B_up_cf  ] = analytic_linearization(x_up,   params);
[A_down_cf, B_down_cf] = analytic_linearization(x_down, params);

fprintf('  Upright eq.  ||A_num - A_analytic||_F = %.2e\n', ...
        norm(A_up - A_up_cf, 'fro'));
fprintf('  Hanging eq.  ||A_num - A_analytic||_F = %.2e\n', ...
        norm(A_down - A_down_cf, 'fro'));

fprintf('  Upright eigenvalues:\n');     disp(eig(A_up).');
fprintf('  Hanging eigenvalues:\n');     disp(eig(A_down).');
fprintf('  Controllability rank (upright) = %d  (expect 4)\n', ...
        rank(ctrb(A_up,  B_up)));
fprintf('  Controllability rank (hanging) = %d  (expect 4)\n\n', ...
        rank(ctrb(A_down,B_down)));

%% ------------------------------------------------------------------------
%  (V4) Animation -- damped free response, played back in real time
%  ------------------------------------------------------------------------
[t_anim, X_anim] = ode45(@(t, x) underactuated_model(t, x, u_zero, params), ...
                         tspan, x0_swing, opts);
animate_cartpole(t_anim, X_anim, params);


% =========================================================================
%   Local helper functions
% =========================================================================
function E = total_energy(x, p)
    s_dot     = x(3);
    theta     = x(2);
    theta_dot = x(4);

    xp_dot =  s_dot + p.l*cos(theta)*theta_dot;
    yp_dot = -p.l*sin(theta)*theta_dot;

    T = 0.5*p.M*s_dot^2 + 0.5*p.m*(xp_dot^2 + yp_dot^2);
    V = p.m*p.g*p.l*cos(theta);

    E = T + V;
end

% -------------------------------------------------------------------------
function [A, B] = numerical_linearization(f, x0, u0)
    n = numel(x0);
    m = numel(u0);
    h = 1e-6;

    A = zeros(n, n);
    for i = 1:n
        e = zeros(n, 1); e(i) = h;
        A(:, i) = (f(x0 + e, u0) - f(x0 - e, u0)) / (2*h);
    end

    B = zeros(n, m);
    for j = 1:m
        e = zeros(m, 1); e(j) = h;
        B(:, j) = (f(x0, u0 + e) - f(x0, u0 - e)) / (2*h);
    end
end

% -------------------------------------------------------------------------
function [A, B] = analytic_linearization(x_eq, p)
% Linearization of  qddot = M(q)^{-1} ( B u - C(q,qdot) qdot - G(q) )
% about an equilibrium with qdot = 0. Then C qdot = 0 and dC/d(.) qdot = 0,
% so only dG/dq matters:  dqddot/dq = -M(q*)^{-1} dG/dq |_{q*}.
    theta = x_eq(2);
    M = p.M; m = p.m; l = p.l; g = p.g;

    Mq = [ M + m,         m*l*cos(theta);
           m*l*cos(theta) m*l^2          ];

    % G(q) = [0; -m g l sin(theta)],  dG/dq = [0 0; 0 -m g l cos(theta)]
    dGdq = [ 0, 0;
             0, -m*g*l*cos(theta) ];

    Bq = [1; 0];

    Minv = Mq \ eye(2);
    A = [zeros(2), eye(2);
         -Minv*dGdq, zeros(2)];
    B = [zeros(2,1);
         Minv*Bq];
end

% -------------------------------------------------------------------------
function animate_cartpole(t, X, p)
    fig = figure('Name', 'Cart-pole animation');
    cart_w = 0.30;  cart_h = 0.20;
    s_min  = min(X(:,1)) - p.l - 0.5;
    s_max  = max(X(:,1)) + p.l + 0.5;

    hold on; axis equal; grid on;
    xlim([s_min, s_max]);
    ylim([-p.l-0.4, p.l+0.4]);
    xlabel('x [m]'); ylabel('y [m]');

    plot([s_min s_max], [-cart_h/2 -cart_h/2], 'k-', 'LineWidth', 1);  % rail
    cart = rectangle('Position', [-cart_w/2 -cart_h/2 cart_w cart_h], ...
                     'Curvature', 0.1, 'FaceColor', [0.85 0.85 0.85]);
    rod  = line([0 0], [0 p.l], 'LineWidth', 2, 'Color', [0 0.3 0.8]);
    bob  = plot(0, p.l, 'o', 'MarkerSize', 10, ...
                'MarkerFaceColor', [0.85 0.1 0.1], ...
                'MarkerEdgeColor', 'k');

    fps    = 60;
    t_play = t(1):1/fps:t(end);
    X_play = interp1(t, X, t_play);

    ttl    = title('');
    t_wall = tic;
    for k = 1:numel(t_play)
        if ~ishghandle(fig); return; end
        s  = X_play(k, 1);
        th = X_play(k, 2);

        cart.Position = [s - cart_w/2, -cart_h/2, cart_w, cart_h];
        bx = s + p.l*sin(th);
        by =     p.l*cos(th);
        set(rod, 'XData', [s, bx], 'YData', [0, by]);
        set(bob, 'XData', bx,      'YData', by);

        set(ttl, 'String', sprintf('t = %5.2f s', t_play(k)));

        dt_sleep = (t_play(k) - t_play(1)) - toc(t_wall);
        if dt_sleep > 0
            pause(dt_sleep);
        else
            drawnow limitrate;
        end
    end
end
