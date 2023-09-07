% Adam-Bashford 4, by David Ahnlund and Emil Gestsson SF2520

% u0 should be a column vector i.e. [ - ; - ; ...]

function [t, u] = AB4(dudt, T, u0, h)
    N = ceil(T/h);
    t = zeros(N, 1);
    u = zeros(N, length(u0));

    [~, u_init] = rk3(dudt, h*4, u0, h);
    u(1:4, :) = u_init(1:4, :);

    for k = 4:N
        u_new = ab4step(dudt, t, u, h, k);

        u(k+1,:) = u_new;
        t(k+1) =  k*h;
    end

    
    function u_new = ab4step(dudt, t, u, h, k)
        f0 = dudt(t(k - 0), u(k - 0, :)');
        f1 = dudt(t(k - 1), u(k - 1, :)');
        f2 = dudt(t(k - 2), u(k - 2, :)');
        f3 = dudt(t(k - 3), u(k - 3, :)');

        u_new = u(k, :) + h/24*(55*f0 - 59*f1 + 37*f2 - 9*f3)';
    end

end