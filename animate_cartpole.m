function animate_cartpole(t, X, p)
% ANIMATE_CARTPOLE  Real-time animation of the cart with double pendulum.
%
%   animate_cartpole(t, X, p)
%
%   Inputs:
%     t : N-by-1 time vector from ode45
%     X : N-by-6 state trajectory  [s  th1  th2  sdot  th1dot  th2dot]
%     p : params struct with fields  l1, l2  (rod lengths)
%
%   Frames are paced to wall-clock time using tic/toc + pause so the
%   playback matches real-world duration. When a frame falls behind,
%   `drawnow limitrate` is used to catch up.

    fig = figure('Name', 'Cart with double pendulum');
    cart_w = 0.30;  cart_h = 0.20;

    L_total = p.l1 + p.l2;
    s_min   = min(X(:,1)) - L_total - 0.3;
    s_max   = max(X(:,1)) + L_total + 0.3;

    hold on; axis equal; grid on;
    xlim([s_min, s_max]);
    ylim([-L_total - 0.3, L_total + 0.3]);
    xlabel('x [m]'); ylabel('y [m]');

    plot([s_min s_max], [-cart_h/2 -cart_h/2], 'k-', 'LineWidth', 1);  % rail
    cart = rectangle('Position', [-cart_w/2 -cart_h/2 cart_w cart_h], ...
                     'Curvature', 0.1, 'FaceColor', [0.85 0.85 0.85]);

    rod1 = line([0 0], [0 p.l1], 'LineWidth', 2, 'Color', [0.00 0.30 0.80]);
    bob1 = plot(0, p.l1, 'o', 'MarkerSize', 9, ...
                'MarkerFaceColor', [0.95 0.65 0.10], ...
                'MarkerEdgeColor', 'k');

    rod2 = line([0 0], [p.l1 p.l1 + p.l2], 'LineWidth', 2, ...
                'Color', [0.55 0.10 0.55]);
    bob2 = plot(0, p.l1 + p.l2, 'o', 'MarkerSize', 10, ...
                'MarkerFaceColor', [0.85 0.10 0.10], ...
                'MarkerEdgeColor', 'k');

    fps    = 60;
    t_play = t(1):1/fps:t(end);
    X_play = interp1(t, X, t_play);

    ttl    = title('');
    t_wall = tic;
    for k = 1:numel(t_play)
        if ~ishghandle(fig); return; end
        s   = X_play(k, 1);
        th1 = X_play(k, 2);
        th2 = X_play(k, 3);

        cart.Position = [s - cart_w/2, -cart_h/2, cart_w, cart_h];

        b1x = s   + p.l1*sin(th1);
        b1y =       p.l1*cos(th1);
        b2x = b1x + p.l2*sin(th2);
        b2y = b1y + p.l2*cos(th2);

        set(rod1, 'XData', [s,   b1x], 'YData', [0,   b1y]);
        set(bob1, 'XData', b1x,        'YData', b1y);
        set(rod2, 'XData', [b1x, b2x], 'YData', [b1y, b2y]);
        set(bob2, 'XData', b2x,        'YData', b2y);

        set(ttl, 'String', sprintf('t = %5.2f s', t_play(k)));

        dt_sleep = (t_play(k) - t_play(1)) - toc(t_wall);
        if dt_sleep > 0
            pause(dt_sleep);
        else
            drawnow limitrate;
        end
    end
end
