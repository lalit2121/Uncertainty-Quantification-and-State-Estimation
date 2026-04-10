clear; close all; clc;
x0 = [0.4 0.3]'+0.1*randn(2,1);
mu = 5;
options = odeset('abstol',1e-8,'reltol',1e-8);
[t,x] = ode45(@(t,x) VanderPol(t,x,mu),[0 100],x0,options);
plot(x(:,1),x(:,2)); hold on;
axis equal


function dx = VanderPol(t,x,mu)
    dx(1,1) = x(2);
    dx(2,1) = mu*(1-x(1)^2)*x(2) - x(1);
end