function dxdt = underactuated_model(t, x, u, params)
% UNDERACTUATED_MODEL  Cart with two-link (serial) pendulum dynamics,
% derived via Euler-Lagrange.
%
%   dxdt = underactuated_model(t, x, u, params)
%
%   System: a cart (mass M) slides horizontally and is actuated by a
%   horizontal force u. A first pendulum (point mass m1, massless rod
%   length l1) is hinged on the cart. A second pendulum (point mass m2,
%   massless rod length l2) is hinged at the tip of the first.
%   3 DOF, 1 input  ->  underactuated.
%
%   Generalized coordinates:
%       q = [s; theta1; theta2]
%   where  s      = horizontal cart position [m]
%          theta1 = angle of pole 1 from upward vertical [rad]
%          theta2 = angle of pole 2 from upward vertical [rad]
%                   (theta_i = 0  -> upright; theta_i = pi -> hanging)
%
%   State vector (used by ode45):
%       x = [s; theta1; theta2; sdot; theta1dot; theta2dot]    (6x1)
%
%   Input:
%       u : horizontal force on the cart [N]
%           (scalar, or function handle u = @(t,x) ...)
%
%   --- Euler-Lagrange derivation ---
%   Cart            : (s, 0)
%   Pole-1 bob (m1) : (s + l1 sin th1,            l1 cos th1)
%   Pole-2 bob (m2) : (s + l1 sin th1 + l2 sin th2,
%                          l1 cos th1 + l2 cos th2)
%
%   Kinetic energy:
%     T = 1/2 M sdot^2
%       + 1/2 m1 (xdot1^2 + ydot1^2)
%       + 1/2 m2 (xdot2^2 + ydot2^2)
%
%   Potential energy (gravity acts in -y, so V = +sum m_i g y_i):
%     V = (m1 + m2) g l1 cos(th1)  +  m2 g l2 cos(th2)
%
%   Lagrangian L = T - V. Applying the Euler-Lagrange equations gives
%   the manipulator form
%
%       M(q) qddot + C(q,qdot) qdot + G(q) = B u
%
%   with
%       M(q) = [ M+m1+m2,            (m1+m2) l1 c1,        m2 l2 c2;
%                (m1+m2) l1 c1,      (m1+m2) l1^2,         m2 l1 l2 c12;
%                m2 l2 c2,           m2 l1 l2 c12,         m2 l2^2 ]
%
%       C(q,qdot) = [ 0,  -(m1+m2) l1 s1 th1d,    -m2 l2 s2 th2d;
%                     0,   0,                      m2 l1 l2 s12 th2d;
%                     0,  -m2 l1 l2 s12 th1d,      0 ]
%
%       G(q) = [ 0;
%               -(m1+m2) g l1 s1;
%               -m2 g l2 s2 ]
%
%       B = [1; 0; 0]    (rank(B) < dim(q)  =>  underactuated)
%
%   where s1 = sin(th1), c1 = cos(th1), s2 = sin(th2), c2 = cos(th2),
%   s12 = sin(th1 - th2), c12 = cos(th1 - th2), th1d = theta1dot, etc.
%
%   Optional viscous friction is added as
%       F = diag([b_cart, b_pend1, b_pend2]) qdot
%   if those fields are present in params (defaults to zero otherwise).

    q    = x(1:3);
    qdot = x(4:6);

    if isa(u, 'function_handle')
        u_val = u(t, x);
    else
        u_val = u;
    end

    [Mq, Cq, Gq, Bq] = manipulator_matrices(q, qdot, params);

    b_cart  = get_or_default(params, 'b_cart',  0);
    b_pend1 = get_or_default(params, 'b_pend1', 0);
    b_pend2 = get_or_default(params, 'b_pend2', 0);
    F_fric  = [b_cart; b_pend1; b_pend2] .* qdot;

    qddot = Mq \ (Bq * u_val - Cq * qdot - Gq - F_fric);

    dxdt = [qdot; qddot];
end

% -------------------------------------------------------------------------
function [Mq, Cq, Gq, Bq] = manipulator_matrices(q, qdot, p)
    th1  = q(2);     th2  = q(3);
    th1d = qdot(2);  th2d = qdot(3);

    M  = p.M;   m1 = p.m1;  m2 = p.m2;
    l1 = p.l1;  l2 = p.l2;  g  = p.g;

    s1  = sin(th1);          c1  = cos(th1);
    s2  = sin(th2);          c2  = cos(th2);
    s12 = sin(th1 - th2);    c12 = cos(th1 - th2);

    Mq = [ M + m1 + m2,        (m1 + m2)*l1*c1,     m2*l2*c2;
           (m1 + m2)*l1*c1,    (m1 + m2)*l1^2,      m2*l1*l2*c12;
           m2*l2*c2,           m2*l1*l2*c12,        m2*l2^2 ];

    Cq = [ 0,  -(m1 + m2)*l1*s1*th1d,    -m2*l2*s2*th2d;
           0,   0,                        m2*l1*l2*s12*th2d;
           0,  -m2*l1*l2*s12*th1d,        0 ];

    Gq = [ 0;
          -(m1 + m2)*g*l1*s1;
          -m2*g*l2*s2 ];

    Bq = [ 1; 0; 0 ];
end

% -------------------------------------------------------------------------
function v = get_or_default(s, name, default)
    if isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default;
    end
end
