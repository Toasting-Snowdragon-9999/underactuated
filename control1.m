function u = control1(t, x, params)
% CONTROL1  Local LQR stabilizer around the upright equilibrium of the
% cart with double pendulum.
%
%   u = control1(t, x, params)
%
%   Linearizes the nonlinear model
%       xdot = f(x, u)        (from underactuated_model.m)
%   numerically about the upright equilibrium  x* = 0,  u* = 0,  then
%   solves the continuous-time algebraic Riccati equation via lqr() and
%   returns the state-feedback control  u = -K x.
%
%   This is only valid LOCALLY near the upright (theta_i ~= 0). For
%   global swing-up, use control2.
%
%   State : x = [s; theta1; theta2; sdot; theta1dot; theta2dot]
%   Cost  : Q = diag([10, 100, 100, 1, 10, 10]),  R = 1.
%
%   The gain K depends only on params, so it is cached across calls.

    persistent K_cached params_cached

    if isempty(K_cached) || ~isequal(params_cached, params)
        x_eq = zeros(6, 1);
        u_eq = 0;

        f = @(xx, uu) underactuated_model(0, xx, uu, params);
        [A, B] = fd_linearize(f, x_eq, u_eq);

        Q = diag([10, 100, 100, 1, 10, 10]);
        R = 1;

        K_cached      = lqr(A, B, Q, R);
        params_cached = params;
    end

    u = -K_cached * x;
    u = u(1);   % ensure scalar
end

% -------------------------------------------------------------------------
function [A, B] = fd_linearize(f, x0, u0)
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
