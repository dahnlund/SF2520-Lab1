% Part 1
a = 1/4 * [1, sqrt(11), 2]';
alpha = 0.07;

%-----------------------------------------------
%Example ode to test rk3 funciton
dudt = @(t,u) [u(2);normpdf(t,0,0.01)-0.2*u(1)];
h = 0.01;

[t, u_rk] = rk3(dudt, 10, [1;1], h);

[~, u_ode] = ode45(dudt, 0:h:10, [1;1]);

disp(norm(u_rk(:,1)-u_ode(:,1)))

plot(t, u_rk(:,1));
hold on
plot(t, u_ode(:,1));
%-----------------------------------------------

%% C)

% standard matrix form using cross product matrix
a_cross = [0, -a(3), a(2) ;
           a(3), 0 , -a(1);                 
           -a(2), a(1), 0];

A = a_cross + alpha * a_cross ^ 2;
lambdas = eig(A);

% find stable region root along eigenvalue line with h in (0, T] 
s_region = @(z) abs(1 + z + z^2/2 + z^3/6) - 1;
find_h0 = @(lambda) fzero(@(h) s_region(h*lambda), [eps(1), T]);

h0_candidates = arrayfun(find_h0, lambdas(abs(lambdas) > eps(1)));
h0 = min(h0_candidates);


