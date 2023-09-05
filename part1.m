% Part 1



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