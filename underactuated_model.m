function dxdt = underactuated_model(t, x, u, params)
% UNDERACTUATED_MODEL  Cart-pole dynamics derived via Euler-Lagrange.
%
%   dxdt = underactuated_model(t, x, u, params)
%
%   System: a pendulum (point mass m, massless rod length l) hinged on a
%   cart (mass M) that slides horizontally. The cart is actuated by a
%   horizontal force u; the pendulum joint is passive -> 2 DOF, 1 input.
%
%   Generalized coordinates:
%       q = [s; theta]
%   where  s     = horizontal cart position [m]
%          theta = pendulum angle from upward vertical [rad]
%                  (theta = 0   -> upright, unstable equilibrium)
%                  (theta = pi  -> hanging, stable equilibrium)
%
%   State vector (used by ode45):
%       x = [s; theta; sdot; thetadot]    (4x1)
%
%   Input:
%       u : horizontal force on the cart [N]
%           (scalar, or function_handle u = @(t,x) ...)
%
%   --- Euler-Lagrange derivation ---
%   Cart position    : (s, 0)
%   Pendulum bob pos : (s + l sin(theta),  l cos(theta))
%
%   Kinetic energy:
%     T = 1/2 M sdot^2
%       + 1/2 m [ (sdot + l cos(theta) thetadot)^2 + (l sin(theta) thetadot)^2 ]
%
%   Potential energy:
%     V = m g l cos(theta)
%
%   Lagrangian L = T - V. Applying
%     d/dt( dL/dqdot ) - dL/dq = Q,    Q = [u; 0]
%   yields the standard manipulator form
%
%       M(q) qddot + C(q,qdot) qdot + G(q) = B u
%
%   with
%       M(q)      = [ M+m,           m l cos(theta);
%                     m l cos(theta) m l^2          ]
%       C(q,qdot) = [ 0, -m l sin(theta) thetadot;
%                     0,  0                        ]
%       G(q)      = [ 0; -m g l sin(theta) ]
%       B         = [ 1; 0 ]    (rank(B) < dim(q)  =>  underactuated)
%
%   Optional viscous friction is added as F = diag([b_cart, b_pend]) qdot
%   if those fields are present in params (defaults to zero otherwise).

    q    = x(1:2);
    qdot = x(3:4);

    if isa(u, 'function_handle')
        u_val = u(t, x);
    else
        u_val = u;
    end

    [Mq, Cq, Gq, Bq] = manipulator_matrices(q, qdot, params);

    b_cart = get_or_default(params, 'b_cart', 0);
    b_pend = get_or_default(params, 'b_pend', 0);
    F_fric = [b_cart; b_pend] .* qdot;

    qddot = Mq \ (Bq * u_val - Cq * qdot - Gq - F_fric);

    dxdt = [qdot; qddot];
end

% -------------------------------------------------------------------------
function [Mq, Cq, Gq, Bq] = manipulator_matrices(q, qdot, p)
    theta    = q(2);
    thetadot = qdot(2);

    M = p.M;  m = p.m;  l = p.l;  g = p.g;
    s_th = sin(theta);  c_th = cos(theta);

    Mq = [ M + m,      m*l*c_th;
           m*l*c_th,   m*l^2     ];

    Cq = [ 0, -m*l*s_th*thetadot;
           0,  0                  ];

    Gq = [ 0;
          -m*g*l*s_th ];

    Bq = [ 1; 0 ];
end

% -------------------------------------------------------------------------
function v = get_or_default(s, name, default)
    if isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default;
    end
end
