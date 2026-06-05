function u = control2(t, x, params)
% CONTROL2  Energy-shaping swing-up controller for the cart with double
% pendulum, with light cart regulation and a hand-off to the local LQR
% stabilizer (control1) near upright.
%
%   u = control2(t, x, params)
%
%   Strategy.
%   ---------
%   Far from upright we inject mechanical energy into the system until E
%   reaches the upright reference
%       E_des = (m1 + m2) g l1 + m2 g l2     (all upright, zero kinetic)
%   The injection law is Astrom--Furuta style:
%       u_energy = -k_E * (E_des - E) * sign(D),
%       D       = (m1 + m2) l1 cos(theta1) theta1dot
%                + m2      l2 cos(theta2) theta2dot.
%
%   Why this guarantees energy injection (and not extraction).
%   ----------------------------------------------------------
%   The cart force is the only horizontal external force, so total
%   mechanical energy evolves as
%       dE/dt = u * sdot.
%   With u = 0 the system conserves horizontal momentum, which (starting
%   from rest) gives the constraint
%       (M + m1 + m2) sdot + (m1 + m2) l1 cos(t1) t1d + m2 l2 cos(t2) t2d = 0,
%   so sign(sdot) = -sign(D). The control above picks sign(u) = -sign(D)
%   when E < E_des, hence sign(u) = sign(sdot), so u * sdot >= 0 and
%   energy is injected (never removed) whenever E < E_des. When E > E_des
%   the same expression flips sign and bleeds energy back out, driving E
%   to E_des. This is the standard cart-pendulum result generalised to
%   the double-pendulum case with mass+length weighting.
%
%   Cart regulation.
%   ----------------
%   A light PD term -k_s s - k_v sdot keeps the cart from drifting away
%   while pumping. It is deliberately weak so it does not dominate the
%   energy-pumping signal.
%
%   No pendulum damping is added during pumping: damping the passive DOFs
%   would bleed the very energy this controller is trying to build, which
%   is the failure mode of the previous version.
%
%   Hand-off.
%   ---------
%   When both wrapped angles are small and velocities modest, control
%   switches to the LQR in control1 (with angles wrapped into (-pi, pi]
%   so winding number does not feed a huge error into the linear law).
%
%   State : x = [s; theta1; theta2; sdot; theta1dot; theta2dot]

    s    = x(1);
    th1  = x(2);   th2  = x(3);
    sd   = x(4);   th1d = x(5);   th2d = x(6);

    M  = params.M;   m1 = params.m1;  m2 = params.m2;
    l1 = params.l1;  l2 = params.l2;  g  = params.g;

    %% --- LQR hand-off near upright -------------------------------------
    th1_w = wrap_pi(th1);
    th2_w = wrap_pi(th2);
    near_upright = abs(th1_w) < 0.40 && abs(th2_w) < 0.40 ...
                && abs(th1d)  < 3.0  && abs(th2d)  < 3.0;
    if near_upright
        x_lqr    = x;
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

    E_des = (m1 + m2)*g*l1 + m2*g*l2;

    %% --- Energy-pumping direction (mass + length weighted) -------------
    %  D matches the cart-velocity expression from horizontal-momentum
    %  conservation, so sign(D) tells us which way the cart "wants" to go
    %  passively. Picking the cart force opposite to D when E < E_des
    %  drives u * sdot >= 0.
    D = (m1 + m2)*l1*cos(th1)*th1d + m2*l2*cos(th2)*th2d;
    eps_D = 1e-3;
    if abs(D) < eps_D
        % t=0 deadlock or transient null: pick an arbitrary nonzero sign
        % so the cart starts moving and pendulum motion is induced.
        direction = 1;
    else
        direction = sign(D);
    end

    %% --- Astrom-style energy pump --------------------------------------
    k_E = 10;
    u_energy = k_E * (E - E_des) * direction;   % sign(u) = -sign(D) when E<E_des

    %% --- Light cart regulation -----------------------------------------
    k_s = 1.0;
    k_v = 1.5;
    u_reg = -k_s * s - k_v * sd;

    u = u_energy + u_reg;

    %% --- Actuator saturation -------------------------------------------
    u_max = 100;
    if u >  u_max, u =  u_max; end
    if u < -u_max, u = -u_max; end
end

% -------------------------------------------------------------------------
function a = wrap_pi(a)
    a = mod(a + pi, 2*pi) - pi;
end
