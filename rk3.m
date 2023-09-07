% Runge-Kutta 3, by David Ahnlund and Emil Gestsson SF2520

% u0 should be a column vector i.e. [ - ; - ; ...]

function [t, u] = rk3(dudt, T, u0, h)
    
    t = zeros(T/h, 1);
    u = zeros(T/h, length(u0));
    u(1,:) = u0';
    N = T/h;

    for k = 1:N
        tk = t(k);
        uk = u(k,:);
        u_new = rk3step(dudt, tk, uk', h);

        u(k+1,:) = u_new;
        t(k+1) =  k*h;
    end


    function u_new = rk3step(dudt, tk, uk, h)
        k1 = dudt(tk, uk);
        k2 = dudt(tk + h, uk + h*k1);
        k3 = dudt(tk + h/2, uk + h*k1/4 + h*k2/4);
        u_new = uk + h/6*(k1 + k2 + 4*k3);
    end

end