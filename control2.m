function u = control2(t, x, params)
% CONTROL2  Energy-shaping swing-up controller for the cart with double
% pendulum, with cart regulation and a hand-off to the local LQR
% stabilizer near upright.
%
%   u = control2(t, x, params)
%
%   Strategy:
%     - Far from upright: pump total mechanical energy E toward its
%       upright value E_des, using
%           u_E   = k_E * (E - E_des) * sign( th1d cos(th1) + th2d cos(th2) )
%       Add cart regulation
%           u_reg = -k_s * s - k_v * sdot
%       so the cart cannot drift off while pumping. (Without this the
%       chaotic double-pendulum swing-up easily walks the cart out of
%       any useful range and stalls.)
%     - Near upright (small wrapped angles AND modest velocities): hand
%       off to the LQR controller from control1, which exponentially
%       stabilizes the upright equilibrium. Angles are unwrapped to
%       [-pi, pi] before being passed to the LQR -- otherwise a chain
%       that has rotated through 2*pi during pumping would feed a huge
%       (and wrong-direction) error into the linear feedback.
%
%   State : x = [s; theta1; theta2; sdot; theta1dot; theta2dot]

    s    = x(1);
    th1  = x(2);   th2  = x(3);
    sd   = x(4);   th1d = x(5);   th2d = x(6);

    M  = params.M;   m1 = params.m1;  m2 = params.m2;
    l1 = params.l1;  l2 = params.l2;  g  = params.g;

    %% --- Switch to LQR when near upright -------------------------------
    th1_w = wrap_pi(th1);
    th2_w = wrap_pi(th2);
    near_upright = abs(th1_w) < 0.35 && abs(th2_w) < 0.35 ...
                && abs(th1d)  < 2.0  && abs(th2d)  < 2.0;
    if near_upright
        % IMPORTANT: pass wrapped angles so LQR sees a small error even
        % if the chain has rotated through full revolutions during pumping.
        x_lqr = x;
        x_lqr(2) = th1_w;
        x_lqr(3) = th2_w;
        u = control1(t, x_lqr, params);
        return;
    end

    %% --- Total mechanical energy ---------------------------------------
    x1d =  sd + l1*cos(th1)*th1d;
    y1d =     - l1*sin(th1)*th1d;
    x2d =  sd + l1*cos(th1)*th1d + l2*cos(th2)*th2d;
    y2d =     - l1*sin(th1)*th1d - l2*sin(th2)*th2d;

    T = 0.5*M*sd^2 ...
      + 0.5*m1*(x1d^2 + y1d^2) ...
      + 0.5*m2*(x2d^2 + y2d^2);
    V = (m1 + m2)*g*l1*cos(th1) + m2*g*l2*cos(th2);
    E = T + V;

    % Upright (kinetic = 0, cos = 1)
    E_des = (m1 + m2)*g*l1 + m2*g*l2;

    %% --- Energy-pumping law + cart regulation --------------------------
    k_E = 12;     % energy-pumping gain      (was 5)
    k_s = 2;      % cart-position regulator
    k_v = 2;      % cart-velocity damping

    % Configuration / damping shaping gains (bias toward upright,
    % suppress residual pendulum oscillations).
    k_theta = 5;
    alpha   = 0.5;
    k_d     = 2;

    direction = sign(th1d*cos(th1) + th2d*cos(th2));
    if direction == 0
        direction = 1;     % avoid being stuck at sign = 0 at t = 0
    end

    u_energy = k_E * (E - E_des) * direction;
    u_reg    = -k_s * s - k_v * sd;
    u_shape  = -k_theta * (sin(th1) + alpha * sin(th2));
    u_damp   = -k_d * (th1d + th2d);

    u = u_energy + u_reg + u_shape + u_damp;

    % Soft saturation keeps swing-up well-behaved for ode45.
    u_max = 200;
    if u >  u_max, u =  u_max; end
    if u < -u_max, u = -u_max; end
end

% -------------------------------------------------------------------------
function a = wrap_pi(a)
    a = mod(a + pi, 2*pi) - pi;
end
