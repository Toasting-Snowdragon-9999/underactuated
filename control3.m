function u = control3(t, x, params)
% CONTROL3  Sliding-mode controller (SMC) for the cart with double
% pendulum, using a collocated sliding surface on the actuated cart
% coordinate s. The two pendulums are stabilized only indirectly,
% through the cart motion.
%
%   u = control3(t, x, params)
%
%   State:  x = [s; theta1; theta2; sdot; theta1dot; theta2dot]
%
%   Sliding surface:  sigma = (sdot - sdot_des) + lambda * (s - s_des)
%   Control law:      u = u_eq + u_sw
%                     u_eq = -(Kp s + Kd sdot)        (PD-style baseline)
%                     u_sw = -k * sat(sigma / phi)    (smoothed switching)
%
%   The boundary-layer saturation replaces sign() so chattering is bounded
%   and the integrator stays well-behaved.

    s   = x(1);
    sd  = x(4);

    %% --- desired trajectory --------------------------------------------
    s_des  = 0;
    sd_des = 0;

    %% --- sliding surface -----------------------------------------------
    lambda = 2;
    e      = s  - s_des;
    edot   = sd - sd_des;
    sigma  = edot + lambda * e;

    %% --- equivalent (PD-baseline) control ------------------------------
    Kp = 10;
    Kd = 5;
    u_eq = -(Kp * s + Kd * sd);

    %% --- switching control (with boundary layer to avoid chatter) ------
    k   = 20;
    phi = 0.1;     % phi > 0 by construction -- no division-by-zero risk
    u_sw = -k * sat(sigma / phi);

    %% --- final control -------------------------------------------------
    u = u_eq + u_sw;
end

% -------------------------------------------------------------------------
function y = sat(x)
    y = max(-1, min(1, x));
end
