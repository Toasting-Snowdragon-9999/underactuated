function u = control3(t, x, params)
% CONTROL3  Sliding-mode controller (SMC) that stabilizes the cart with
% double pendulum at the upright equilibrium (tip of the second pendulum
% pointing straight up).
%
%   u = control3(t, x, params)
%
%   State:  x = [s; theta1; theta2; sdot; theta1dot; theta2dot]
%   Goal :  drive (theta1, theta2, s) -> 0 and keep them there, i.e. both
%           pendulums upright and the cart at the origin.
%
%   Design. The plant is underactuated (single input, three DOF), so the
%   sliding surface cannot be defined on one coordinate alone -- a
%   collocated surface on s only would leave the pendulum angles
%   unconstrained on the sliding manifold. We therefore use a sliding
%   surface that mixes the actuated cart coordinate with both
%   pendulum-angle errors:
%
%       sigma = K * x_wrap
%
%   where K is the LQR gain for the linearization about the upright
%   equilibrium and x_wrap is x with theta1, theta2 wrapped into
%   (-pi, pi]. On the manifold sigma = 0 the reduced dynamics inherit the
%   LQR closed-loop stability around upright, so the entire state goes to
%   zero (tip upright, cart centered).
%
%   The control law has the standard equivalent-plus-switching form
%
%       u = u_eq + u_sw,
%       u_eq = -K * x_wrap,                        (linear feedback)
%       u_sw = -k_sw * sat(sigma / phi),           (boundary-layer switch)
%
%   The saturation replaces sign() to keep chattering bounded while
%   preserving the disturbance-rejection character of SMC. The control is
%   clipped to a physical actuator limit.
%
%   K is cached across calls (depends only on params).

    persistent K_cached params_cached

    if isempty(K_cached) || ~isequal(params_cached, params)
        x_eq = zeros(6, 1);
        u_eq0 = 0;

        f = @(xx, uu) underactuated_model(0, xx, uu, params);
        [A, B] = fd_linearize(f, x_eq, u_eq0);

        Q = diag([10, 200, 200, 1, 10, 10]);
        R = 1;

        K_cached      = lqr(A, B, Q, R);
        params_cached = params;
    end

    %% --- wrap pendulum angles into (-pi, pi] ---------------------------
    x_w    = x;
    x_w(2) = wrap_pi(x(2));
    x_w(3) = wrap_pi(x(3));

    %% --- sliding surface (mixes cart + both pendulums) -----------------
    sigma = K_cached * x_w;

    %% --- equivalent control (LQR baseline) -----------------------------
    u_eq = -K_cached * x_w;

    %% --- switching control with boundary layer -------------------------
    k_sw = 8;
    phi  = 0.4;
    u_sw = -k_sw * sat(sigma / phi);

    %% --- final control with actuator saturation ------------------------
    u = u_eq + u_sw;

    u_max = 100;
    if u >  u_max, u =  u_max; end
    if u < -u_max, u = -u_max; end
end

% -------------------------------------------------------------------------
function y = sat(x)
    y = max(-1, min(1, x));
end

% -------------------------------------------------------------------------
function a = wrap_pi(a)
    a = mod(a + pi, 2*pi) - pi;
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
