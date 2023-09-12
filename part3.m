% Part 3

%Define rate constants

r1 = 5e-2; r2 = 1.2e-4; r3 = 4e7;

%% A

dxdt = @(t,x) [-r1*x(1) + r2*x(2)*x(3); r1*x(1) - r2*x(2)*x(3) - r3*x(2)^2; r3*x(2)^2];

T = 100;

[t, u_rk] = rk3(dxdt, 1000, [1;0;0], 0.0001);

plot(t, u_rk)

