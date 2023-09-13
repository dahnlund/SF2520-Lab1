% explicit euler, by David Ahnlund and Emil Gestsson SF2520

% u0 should be a column vector i.e. [ - ; - ; ...]

function [t, u] = expeuler(dudt, T, u0, h)
    N = ceil(T/h);
    t = zeros(N, 1);
    u = zeros(N, length(u0));
    u(1,:) = u0';

    for k = 1:N
        tk = t(k);
        uk = u(k,:);
        u_new = uk + h * dudt(tk, uk')';

        u(k+1,:) = u_new;
        t(k+1) =  k*h;
    end
end