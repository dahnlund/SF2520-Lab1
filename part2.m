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

T = 10;
[t, u] = AB4(dudt, T, u0, 0.0001);

fprintf("Number of time steps AB4: %.0f\n", length(t))

rx = u(:, 1); ry = u(:, 2);
drx= u(:, 3); dry= u(:, 4);

[t2, u2] = AB4(dudt, T, u0, 0.001);
rx2 = u2(:, 1); ry2 = u2(:, 2);
drx2= u2(:, 3); dry2= u2(:, 4);

figure
subplot(1,2, 1)
plot(rx, ry);
hold on 
plot(rx2, ry2);
hold off
legend("h = 10^{-4}", "h = 10^{-3}")
title("2a) position")
xlabel("x"); ylabel("y");

subplot(1,2, 2)
plot(drx, dry);
hold on 
plot(drx2, dry2);
hold off
legend("h = 10^{-4}", "h = 10^{-3}")
title("2a) velocity")
xlabel("x'_t"); ylabel("y'_t");

timesteps_AB4 = diff(t); % Equals h

%% B

% time, correct values, upper limits of N for each method, and tolerence.
T  = [5, 20, 40];
V  = [0.4681 0.6355;-0.2186 -0.2136;-1.4926 -0.3339];
tol= 0.1;

% functions and thier upper limits on N.
funcs = {@expeuler, @rk3, @AB4};

% store errors and resulting N
res = zeros(3);
errs= zeros(3);
for r = 1:3
    f = funcs{r};
    for c = 1:3
        t = T(c);
        

        % bisection method implementation of lower_bound
        % finds smallest N with error > tol
        lo = 10;
        hi = 1e6;
        while hi - lo > 1
            N = floor((lo + hi) / 2);
            [~, u]= f(dudt, t, u0, t/N);
            err = norm(u(end, 1:2) - V(c, :));
    
            if err < tol
                hi = N;
            else
                lo = N;
            end
        end
        
        res(r, c) = N;
        errs(r,c) = err;

        fprintf("%d &", N);
    end
    fprintf("\\\\\\hline\n");
end

%res =
%      940662      999999      999999
%        3489       40046      646973
%        2975       24739      175505

%errs =
%    0.1000    3.1647    1.5956
%    0.1000    0.1012    0.1000
%    0.0999    0.1000    0.1000

%% C
%Comparing solution with the adaptive - built in - ode23

options = odeset(Reltol = 1e-4);

T = 20;
[t, u] = ode23(dudt, [0 T], u0, options);
disp("Number of time steps ode23:")
disp(length(t))

rx_23 = u(:, 1); ry_23 = u(:, 2);
drx_23 = u(:, 3); dry_23 = u(:, 4);

error_23 = norm([-0.2186; -0.2136] - [rx_23(end);ry_23(end)]); %Verify accuracy
fprintf("Error (L2) from ODE23: %.4f\n", error_23)

timesteps_23 = diff(t);

min_timestep23 = min(timesteps_23);
max_timestep23 = max(timesteps_23);
fprintf("ode23 max time step: %.4f. Min time step: %.6f\n", max_timestep23, min_timestep23)
plot(t(2:end), timesteps_23)
title("Timestep size in ode23")
xlabel("Time, t")
ylabel("Step size, h")

plot(t, rx_23); hold on;plot(t, ry_23)
legend("x-component of r", "y-component of r")
title("r-components as a funciton of time")
ylabel("r")
xlabel("Time, t")