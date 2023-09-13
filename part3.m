% Part 3

%Define rate constants

r1 = 5e-2; r2 = 1.2e4; r3 = 4e7;

%% A

dxdt = @(t,x) [-r1*x(1) + r2*x(2)*x(3); r1*x(1) - r2*x(2)*x(3) - r3*x(2)^2; r3*x(2)^2];

T = 100;

h = 3e-4:5e-5:6e-4;    % For T = 10, this is an appropriate interval to test: h = 7e-4:1e-6:8e-4
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

J = @(a,b,c) [-r1, r2*c r2*b; r1, -r2*c-2*r3*b, -r2*b; 0, 2*r3*b, 0];

eigs = zeros(size(u_rk));

for i = 1:length(u_rk)
    
    eigs(i,:) = eig(J(u_rk(i,1),u_rk(i,2),u_rk(i,3)));

end

plot(eigs(:,1))
figure
plot(eigs(:,3))
