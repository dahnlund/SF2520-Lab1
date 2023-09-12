% Part 3

%Define rate constants

r1 = 5e-2; r2 = 1.2e4; r3 = 4e7;

%% A

dxdt = @(t,x) [-r1*x(1) + r2*x(2)*x(3); r1*x(1) - r2*x(2)*x(3) - r3*x(2)^2; r3*x(2)^2];


h = 8e-4:1e-6:1e-3;
max_vals = zeros(length(h), 1);
for i = 1:length(h)
    [t, u_rk] = rk3(dxdt, 10, [1;0;0], h(i));
    
    max_vals(i) = max(norm(u_rk));
    
    if isnan(max_vals(i))
        disp("h breakpoint")
        disp(h(i))
        break
    end
end
plot(t, u_rk)
figure
semilogy(t,u_rk)

