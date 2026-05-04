% SIMULATION  Simulate and verify the cart-double-pendulum Euler-Lagrange model.
%
% Verification strategy:
%   (V1) Energy conservation. With u = 0 and no friction, total mechanical
%        energy E = T + V must be constant along trajectories. We integrate
%        the model and check the relative drift.
%   (V2) Equilibria. Initial conditions at the upright-upright (th1=th2=0)
%        and hanging-hanging (th1=th2=pi) configurations, with zero
%        velocity and u=0, must produce zero acceleration.
%   (V3) Linearization. Numerical Jacobians around the two equilibria are
%        compared against the closed-form linearizations from the
%        manipulator equations. Eigenvalues are reported and the
%        controllability rank is checked.
%   (V4) Animation. A real-time visual sanity check of the swing.

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

params.b_cart  = 4.0;    % cart viscous friction    [N s / m]    -- vicious
params.b_pend1 = 0.25;   % pole-1 viscous friction  [N m s / rad]
params.b_pend2 = 0.20;   % pole-2 viscous friction  [N m s / rad]

%% ------------------------------------------------------------------------
%  (V1) Energy conservation  --  free response from a perturbed hang
%  ------------------------------------------------------------------------
%  state = [s; th1; th2; sdot; th1dot; th2dot]
x0_swing = [0; pi - 0.30; pi + 0.20; 0; 0; 0];
tspan    = [0 10];
opts     = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

u_zero = @(t, x) 0;

% V1 must run friction-free regardless of the global params.
params_nf = params;
params_nf.b_cart  = 0;
params_nf.b_pend1 = 0;
params_nf.b_pend2 = 0;

[t1, X1] = ode45(@(t, x) underactuated_model(t, x, u_zero, params_nf), ...
                 tspan, x0_swing, opts);

E1   = arrayfun(@(k) total_energy(X1(k,:).', params_nf), 1:size(X1,1));
relE = (E1 - E1(1)) / max(abs(E1(1)), eps);

fprintf('--- (V1) Energy conservation (u = 0, friction = 0) ---\n');
fprintf('  E(0)              = %+ .6f  J\n', E1(1));
fprintf('  max |E(t) - E(0)| = % .3e J\n', max(abs(E1 - E1(1))));
fprintf('  max relative drift= % .3e\n\n', max(abs(relE)));

figure('Name', '(V1) States, free response');
subplot(3,2,1); plot(t1, X1(:,1)); grid on;
    xlabel('t [s]'); ylabel('s [m]');           title('cart position');
subplot(3,2,3); plot(t1, X1(:,2)); grid on;
    xlabel('t [s]'); ylabel('\theta_1 [rad]');  title('pole-1 angle');
subplot(3,2,5); plot(t1, X1(:,3)); grid on;
    xlabel('t [s]'); ylabel('\theta_2 [rad]');  title('pole-2 angle');
subplot(3,2,2); plot(t1, X1(:,4)); grid on;
    xlabel('t [s]'); ylabel('sdot [m/s]');           title('cart velocity');
subplot(3,2,4); plot(t1, X1(:,5)); grid on;
    xlabel('t [s]'); ylabel('\theta_1 dot [rad/s]'); title('pole-1 rate');
subplot(3,2,6); plot(t1, X1(:,6)); grid on;
    xlabel('t [s]'); ylabel('\theta_2 dot [rad/s]'); title('pole-2 rate');

figure('Name', '(V1) Energy drift');
plot(t1, E1 - E1(1), 'LineWidth', 1.2); grid on;
xlabel('t [s]'); ylabel('E(t) - E(0)  [J]');
title(sprintf('Total energy drift  (max rel. err = %.2e)', max(abs(relE))));

%% ------------------------------------------------------------------------
%  (V2) Equilibria
%  ------------------------------------------------------------------------
fprintf('--- (V2) Equilibria ---\n');
x_up   = [0; 0;  0;  0; 0; 0];   % both upright
x_down = [0; pi; pi; 0; 0; 0];   % both hanging

f_up   = underactuated_model(0, x_up,   0, params);
f_down = underactuated_model(0, x_down, 0, params);

fprintf('  f(both-upright,  u=0) = [% .2e % .2e % .2e % .2e % .2e % .2e]\n', f_up);
fprintf('  f(both-hanging,  u=0) = [% .2e % .2e % .2e % .2e % .2e % .2e]\n\n', f_down);

%% ------------------------------------------------------------------------
%  (V3) Linearization around the two equilibria
%  ------------------------------------------------------------------------
fprintf('--- (V3) Linearization ---\n');

f_xu = @(x, u) underactuated_model(0, x, u, params);

[A_up,   B_up  ] = numerical_linearization(f_xu, x_up,   0);
[A_down, B_down] = numerical_linearization(f_xu, x_down, 0);

% Closed-form linearization from M qddot + G = B u about an equilibrium
% q* with qdot* = 0:  qddot ~= M(q*)^{-1} ( B u - dG/dq |_{q*} dq ).
[A_up_cf,   B_up_cf  ] = analytic_linearization(x_up,   params);
[A_down_cf, B_down_cf] = analytic_linearization(x_down, params);

fprintf('  Both-upright eq.  ||A_num - A_analytic||_F = %.2e\n', ...
        norm(A_up - A_up_cf, 'fro'));
fprintf('  Both-hanging eq.  ||A_num - A_analytic||_F = %.2e\n', ...
        norm(A_down - A_down_cf, 'fro'));

fprintf('  Both-upright eigenvalues:\n');     disp(eig(A_up).');
fprintf('  Both-hanging eigenvalues:\n');     disp(eig(A_down).');
fprintf('  Controllability rank (both-upright) = %d  (expect 6)\n', ...
        rank(ctrb(A_up,   B_up)));
fprintf('  Controllability rank (both-hanging) = %d  (expect 6)\n\n', ...
        rank(ctrb(A_down, B_down)));

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
    s_dot = x(4);
    th1   = x(2);  th2   = x(3);
    th1d  = x(5);  th2d  = x(6);

    % bob-1 velocity
    x1d =  s_dot + p.l1*cos(th1)*th1d;
    y1d =        - p.l1*sin(th1)*th1d;

    % bob-2 velocity
    x2d =  s_dot + p.l1*cos(th1)*th1d + p.l2*cos(th2)*th2d;
    y2d =        - p.l1*sin(th1)*th1d - p.l2*sin(th2)*th2d;

    T = 0.5*p.M*s_dot^2 ...
      + 0.5*p.m1*(x1d^2 + y1d^2) ...
      + 0.5*p.m2*(x2d^2 + y2d^2);
    V = (p.m1 + p.m2)*p.g*p.l1*cos(th1) + p.m2*p.g*p.l2*cos(th2);

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
% about an equilibrium with qdot = 0. Then C qdot = 0; further, C is
% linear in qdot, so (C qdot) is quadratic in qdot and its Jacobians
% w.r.t. q and qdot vanish at qdot = 0. So only dG/dq matters:
%   d(qddot)/dq = -M(q*)^{-1} dG/dq |_{q*}.
    th1 = x_eq(2);  th2 = x_eq(3);

    M  = p.M;   m1 = p.m1;  m2 = p.m2;
    l1 = p.l1;  l2 = p.l2;  g  = p.g;

    c1  = cos(th1);    c2  = cos(th2);
    c12 = cos(th1 - th2);

    Mq = [ M + m1 + m2,        (m1 + m2)*l1*c1,     m2*l2*c2;
           (m1 + m2)*l1*c1,    (m1 + m2)*l1^2,      m2*l1*l2*c12;
           m2*l2*c2,           m2*l1*l2*c12,        m2*l2^2 ];

    % G = [0; -(m1+m2) g l1 sin(th1); -m2 g l2 sin(th2)]
    dGdq = [ 0,  0,                        0;
             0, -(m1 + m2)*g*l1*c1,        0;
             0,  0,                       -m2*g*l2*c2 ];

    Bq = [1; 0; 0];

    Minv = Mq \ eye(3);
    A = [zeros(3),    eye(3);
         -Minv*dGdq,  zeros(3)];
    B = [zeros(3,1);
         Minv*Bq];
end
