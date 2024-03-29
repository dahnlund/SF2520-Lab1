% Part 1
clear
a = 1/4 * [1, sqrt(11), 2]';
alpha = 0.07;

%% A)

m0 = [0;0;1];

dmdt = @(t,m) cross(a,m) + alpha*cross(a, cross(a,m));
[t, m_rk] = rk3(dmdt, 50, m0, 0.001);   % Run the scheme using the rk3.m file where Runge Kutta 3 is defined

%Plot the three components as a function of t

plot(t, m_rk)
title('The three components as a funtion of t')
xlabel('time in seconds')
legend('First','Second', 'Third')

figure 
plot3(m_rk(:,1), m_rk(:,2), m_rk(:,3))
axis equal
xlabel('m(:,1)')
ylabel('m(:,2)')
zlabel('m(:,3)')
title("Plot3 of m")
hold on
plot3([0 a(1)], [0 a(2)], [0 a(3)])


%% B)

N = [50 100 200 400 800 1600 3200 6400 12800];
T = 50;
saved_m_end = zeros(length(N),3);

for i = 1:length(N)
    h = T/N(i);
    
    [~, m_rk] = rk3(dmdt, T, m0, h);
    
    saved_m_end(i,:) = m_rk(end,:);
end

d = diff(saved_m_end);
errors = sqrt(sum(d.^2, 2));

loglog(T./N(1:end-1),errors)

hold on

loglog(T./N,(T./N).^(3))
title("Log-log plot")
legend("|m_N(T) - m_{2N}(T)|", "h^{3}")
xlabel("h")

diff(log2(flip(errors)))

%% C)

% standard matrix form using cross product matrix
a_cross = [0, -a(3), a(2) ;
           a(3), 0 , -a(1);                 
           -a(2), a(1), 0];

A = a_cross + alpha * a_cross ^ 2;
lambdas = eig(A);

T = 50;

% find stable region root along eigenvalue line with h in (0, T] 
s_region = @(z) abs(1 + z + z^2/2 + z^3/6) - 1;
find_h0 = @(lambda) fzero(@(h) s_region(h*lambda), [0.00001, T]);

h0_candidates = arrayfun(find_h0, lambdas(abs(lambdas) > eps(1)));
h0 = min(h0_candidates);


%% D
[t1, m_rk1] = rk3(dmdt, T, m0, h0 - 0.1);
[t3, m_rk3] = rk3(dmdt, T, m0, h0 + 0.1);
[t, m_rk] = rk3(dmdt, 50, m0, 0.001);

figure
title("Solutions with h surrouding the stability limit");
subplot(1, 2, 1);
hold on 
plot3(m_rk(:,1), m_rk(:,2), m_rk(:,3))
plot3(m_rk1(:,1), m_rk1(:,2), m_rk1(:,3))
hold off
title("h = h_0 - 0.1");
legend("h = 0.001", "h = h_0 - 0.1");
xlabel("x"); ylabel("y"); zlabel("z");

subplot(1, 2, 2);
hold on 
plot3(m_rk(:,1), m_rk(:,2), m_rk(:,3))
plot3(m_rk3(:,1), m_rk3(:,2), m_rk3(:,3))
hold off
title("h = h_0 + 0.1")
legend("h = 0.001", "h = h_0 + 0.1");
xlabel("x"); ylabel("y"); zlabel("z");


