function [t y]=impeuler(Tend, y0, h)
%IMPEULER  Solve Robertson problem with the Implicit Euler method.
%   [TOUT,YOUT] = IMPEULER(T,Y0,h) computes the concentrations xA, xB, xC from 
%   time t=0 to t=T using timestep h, and initial data Y0=[xA(0),xB(0),xC(0)]. 
%   The solution is returned in the array YOUT, where each row contains the 
%   concentrations xA, xB, xC at at the time given in the accompanying column 
%   vector TOUT.

r1 = 0.05;  % Rate constants
r2 = 1.2e4;
r3 = 4e7;

f = @(u) [-r1*u(1)+r2*u(2)*u(3); r1*u(1)-r2*u(2)*u(3)-r3*u(2)^2; r3*u(2)^2];  % RHS function
Jf = @(u) [-r1 r2*u(3) r2*u(2); r1 -r2*u(3)-2*r3*u(2) -r2*u(2); 0 2*r3*u(2) 0]; % Jacobian of f
I = eye(3,3);

tol = 1e-10;  % Tolerance in Newton

t = 0:h:Tend;
u = y0(:);    % Make sure u is a column vector
y = u';


for t1=t(2:end)
    d = 1;  % Dummy value
    v = u;  % Starting guess for u_{n+1} taken as u_n
    
    while (norm(d)>tol)    % Solve F(v) = v - h*f(v) - u_n = 0, with Newton's method
        F = v - h*f(v) - u;
        J = I - h*Jf(v);   % Jacobian of F
        d = -J\F;
        v = v + d;
    end
    
    u = v;        % u_{n+1} = v, solution to F(v)=0.
    y = [y; u'];  % Save values for plotting
end
t = t';
end
