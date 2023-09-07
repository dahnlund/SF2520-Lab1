% Part 2

mu = 1/82.45;
r0 = [-mu; 0];
r1 = [1 - mu; 0];
B = [0, 1;-1, 0];

drrdt= @(r, dr) -(1-mu) * (r - r0) / norm(r - r0)^3 ...
                - mu * (r - r1) / norm(r - r1)^3  ...
                + 2 * B * dr + r;

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

T = 20;
[t, u] = AB4(dudt, T, u0, 0.0001);

disp("Number of time steps AB4:")
disp(length(t))

rx = u(:, 1); ry = u(:, 2);
drx= u(:, 3); dry= u(:, 4);

timesteps_AB4 = diff(t); % Equals h

%% C
%Comparing solution with the adaptive - built in - ode23

options = odeset(Reltol = 1e-4);

[t, u] = ode23(dudt, [0 T], u0, options);
disp("Number of time steps ode23:")
disp(length(t))

rx_23 = u(:, 1); ry_23 = u(:, 2);
drx_23 = u(:, 3); dry_23 = u(:, 4);

error_23 = norm([-0.2186; -0.2136] - [rx_23(end);ry_23(end)]); %Verify accuracy

timesteps_23 = diff(t);

min_timestep23 = min(timesteps_23);
max_timestep23 = max(timesteps_23);

plot(t(2:end), timesteps_23)
title("Timestep size in ode23")
