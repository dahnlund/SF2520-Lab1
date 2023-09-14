% Part 3

%Define rate constants

r1 = 5e-2; r2 = 1.2e4; r3 = 4e7;

%% A

dxdt = @(t,x) [-r1*x(1) + r2*x(2)*x(3); r1*x(1) - r2*x(2)*x(3) - r3*x(2)^2; r3*x(2)^2];

T = 10;

h = 7e-4:1e-6:8e-4;    % For T = 10, this is an appropriate interval to test: h = 7e-4:1e-6:8e-4. For T = 100: h = 3e-4:5e-5:6e-4;
max_vals = zeros(length(h), 1);
for i = 1:length(h)
    [t, u_rk] = rk3(dxdt, T, [1;0;0], h(i));

    delta = 0.01;
    [~, u_rk_tilde] = rk3(dxdt, T, [1+delta;0;0], h(i));  %Introduce a small petrubation
    
    e = abs(u_rk - u_rk_tilde)./u_rk;
    
    if norm(max(e)) > 10*delta   % Condition for when the solution is considered unstable
        fprintf("Unstable for h = %d\n", h(i))
        h_max = h(i-1);
        fprintf("Max allowed h is %d\n", h_max)
        break
    end   
    h_max = h(i);
end

[t, u_rk] = rk3(dxdt, T, [1;0;0], h_max);
plot(e)
title("Relative error when \delta is introduced")
ylabel("RelError")
xlabel("n")
figure
plot(t, u_rk)
title("Plot of u")
figure
semilogy(t, u_rk)
title("Semilog-plot for u")


%% B

J = @(x) [-r1, r2*x(3) r2*x(2); r1, -r2*x(3)-2*r3*x(2), -r2*x(2); 0, 2*r3*x(2), 0];

eigs = zeros(size(u_rk));

for i = 1:length(u_rk)
    
    eigs(i,:) = eig(J(u_rk(i,:)));

end

plot(t,eigs(:,1))
figure
plot(t,eigs(:,3))

max_eig = eigs(end,1);  % Seen from the plots

s_condition = @(h) 2 + (h*max_eig) + (h*max_eig).^2/2 + (h*max_eig).^3/6;
h_max_theoretical = fzero(s_condition, 1e-5);
fprintf("Theoretical h_max: %d\n", h_max_theoretical)

%% C

T = 1000;

h = 2e-4;

tic
[t, u_rk] = rk3(dxdt, T, [1;0;0], h);
time_rk = toc;
fprintf("Time to run RK3 for T=1000 where h = %d: %f seconds\n", h,time_rk)

delta = 0.01;
[~, u_rk_tilde] = rk3(dxdt, T, [1+delta;0;0], h);  %Introduce a small petrubation
e = abs(u_rk - u_rk_tilde)./u_rk;
if norm(max(e)) > 10*delta
    disp("Unstable solution")
else
    disp("Stable solution")
end

plot(t, u_rk)
figure
semilogy(t, u_rk)

eigs = zeros(size(u_rk));

for i = 1:length(u_rk)
    
    eigs(i,:) = eig(J(u_rk(i,:)));

end
plot(t,eigs(:,1))
figure
plot(t,eigs(:,3))

max_eig = eigs(end,1);  % Seen from the plots

s_condition = @(h) 2 + (h*max_eig) + (h*max_eig).^2/2 + (h*max_eig).^3/6;
h_max_theoretical = fzero(s_condition, 1e-5);
fprintf("Theoretical h_max: %d\n", h_max_theoretical)

%% D
T = 10;
h_rk3 = 2e-4;
h_ie = 0.2;

[t_imp, u_imp] = impeuler(T, [1;0;0], h_ie);
[t_rk3, u_rk] = rk3(dxdt, T, [1;0;0], h_rk3);

figure
title("3d)")
subplot(2,1,1);
hold on
plot(t_imp, u_imp);
plot(t_rk3, u_rk);
title("comparison between impeuler and rk3 with T = 10.");
hold off

subplot(2,1,2);
T = 1000;
[t_imp, u_imp] = impeuler(T, [1;0;0], h_ie);
plot(t_imp, u_imp);
title("impeuler T = 1000.");

%% E
T = 1000;
correct_u = [0.293414227164 0.000001716342048 0.706584056494];
h_ie = [1, 0.1, 0.01, 0.001];
h_rk3 = 2e-4;

tic
[~, u_rk3] = rk3_noplot(dxdt, T, [1;0;0], h_rk3);
t = toc;

fprintf("RK & %f & %.13f & %f \\\\\\hline\n", h_rk3, norm(u_rk3 - correct_u), t);

for h = h_ie
    tic
    [~, u] = impeuler_noplot(T, [1;0;0], h);
    t = toc;
    
    fprintf("IE & %f & %.13f & %f \\\\\\hline\n", h, norm(u - correct_u), t);
end
