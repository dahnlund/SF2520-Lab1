% Part 2

mu = 1/82.45;
r0 = [-mu; 0];
r1 = [1 - mu; 0];
B = [0, 1;-1, 0];

drrdt= @(r, dr) -(1-mu) * (r - r0) / norm(r - r0)^3 ...
                - mu * (r - r1) / norm(r - r1)^3  ...
                + 2 * B * dr;

% first order ode in the form
% x, y, x', y'
dudt = @(t, u) [
    u(3:4);
    drrdt(u(1:2), u(3:4))
];

u0 = [
% position
  1.15; 0;
% velocity
  0; -0.975
];

T = 10;
[t, u] = AB4(dudt, T, u0, 0.001);

rx = u(:, 1); ry = u(:, 2);
drx= u(:, 3); dry= u(:, 4);




