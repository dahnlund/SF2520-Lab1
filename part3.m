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
        sprintf("Unstable for h = %d", h(i))
        h_max = h(i-1);
        sprintf("Max allowed h is %d", h_max)
        break
    end   
    h_max = h(i);
end

[t, u_rk] = rk3(dxdt, T, [1;0;0], h_max);
plot(e)

figure
plot(t, u_rk)
figure
semilogy(t, u_rk)


%% B

J = @(x) [-r1, r2*x(3) r2*x(2); r1, -r2*x(3)-2*r3*x(2), -r2*x(2); 0, 2*r3*x(2), 0];

eigs = zeros(size(u_rk));

for i = 1:length(u_rk)
    
    eigs(i,:) = eig(J(u_rk(i,:)));

end

plot(t,eigs(:,1))
figure
plot(t,eigs(:,3))

max_eig = eigs(end,1);

% find stable region root along eigenvalue line with h in (0, T] 
s_condition = @(h) 2 + (h*max_eig) + (h*max_eig).^2/2 + (h*max_eig).^3/6;

h_max_theoretical = fzero(s_condition, 1e-5);

