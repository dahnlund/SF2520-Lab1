% Runge-Kutta 3, by David Ahnlund and Emil Gestsson SF2520

% u0 should be a column vector i.e. [ - ; - ; ...]

function [tk, uk] = rk3_noplot(dudt, T, u0, h)
    N = ceil(T/h);
    uk = u0;
    
    for k = 1:N
        uk = rk3step(dudt, (k - 1)*h, uk, h);
    end

    uk = uk';
    tk = N * h;

    function u_new = rk3step(dudt, tk, uk, h)
        k1 = dudt(tk, uk);
        k2 = dudt(tk + h, uk + h*k1);
        k3 = dudt(tk + h/2, uk + h*k1/4 + h*k2/4);
        u_new = uk + h/6*(k1 + k2 + 4*k3);
    end
end