function u = control2(t, x, params)
% CONTROL2  Energy-shaping swing-up controller for the cart with double
% pendulum, with a hand-off to the local LQR stabilizer near upright.
%
%   u = control2(t, x, params)
%
%   Strategy:
%     - Far from upright: apply an energy-pumping cart force that drives
%       the total mechanical energy E toward its value at the upright
%       equilibrium, E_des. This swings the pendulum chain up.
%     - Near upright (small angles AND small velocities): hand off to the
%       LQR controller from control1, which exponentially stabilizes the
%       upright equilibrium.
%
%   Energy-pumping law (gain k, default 5):
%       u = k * (E - E_des) * sign( th1d cos(th1) + th2d cos(th2) )
%
%   Sign of  d/dt (m1 y1 + m2 y2)  is approximated by the sign argument
%   above; pumping in this direction raises the bobs when energy is low
%   and lowers them when energy is high.
%
%   State : x = [s; theta1; theta2; sdot; theta1dot; theta2dot]

    th1  = x(2);   th2  = x(3);
    sd   = x(4);   th1d = x(5);   th2d = x(6);

    M  = params.M;   m1 = params.m1;  m2 = params.m2;
    l1 = params.l1;  l2 = params.l2;  g  = params.g;

    %% --- Switch to LQR when near upright -------------------------------
    th1_w = wrap_pi(th1);
    th2_w = wrap_pi(th2);
    near_upright = abs(th1_w) < 0.2 && abs(th2_w) < 0.2 ...
                && abs(th1d)  < 1.0 && abs(th2d)  < 1.0;
    if near_upright
        u = control1(t, x, params);
        return;
    end

    %% --- Total mechanical energy ---------------------------------------
    % bob velocities (same expressions used in simulation.m total_energy)
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

    %% --- Energy-pumping law --------------------------------------------
    k = 5;
    direction = sign(th1d*cos(th1) + th2d*cos(th2));
    if direction == 0
        direction = 1;     % avoid being stuck at sign = 0
    end

    u = k * (E - E_des) * direction;

    % Soft saturation keeps the swing-up well-behaved for ode45.
    u_max = 100;
    if u >  u_max, u =  u_max; end
    if u < -u_max, u = -u_max; end
end

% -------------------------------------------------------------------------
function a = wrap_pi(a)
    a = mod(a + pi, 2*pi) - pi;
end
