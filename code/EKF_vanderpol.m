clear; close all; clc;
m0 = [0.4 0.3]';
P0 = 0.1^2*eye(2);
S0 = sqrtm(P0);
mu = 2;
R = 0.1^2;
SR = sqrt(R);
H = [1 0];
time = linspace(0,20,201);
options = odeset('abstol',1e-8,'reltol',1e-8);



xtrue = zeros(2,length(time));
xtrue(:,1) = m0 + S0*randn(2,1);
xhat = zeros(2,length(time));
xhat(:,1) = m0;
Phat = P0;

Ps = zeros(2,length(time));
Ps(:,1) = diag(P0);

for i = 2:length(time)
    [~,x] = ode45(@(t,x) VanderPol(t,x,mu),[time(i-1) time(i)],xtrue(:,i-1),options);
    xtrue(:,i) = x(end,:)';
    ytrue = H*xtrue(:,i) + SR*randn(1);

    state = [xhat(:,i-1); reshape(Phat,4,1)];
    [~,state] = ode45(@(t,x) VanderPol_EKF(t,x,mu),[time(i-1) time(i)],state,options);
    state = state(end,:)';
    xhat(:,i) = state(1:2);
    Phat = reshape(state(3:6),2,2);

    yhat = H*xhat(:,i);
    Pxy = Phat*H';
    Pyy = H*Phat*H'+R;
    K = Pxy/Pyy;
    xhat(:,i) = xhat(:,i) + K*(ytrue-yhat);
    Phat = Phat - K*Pyy*K';
    
    Ps(:,i) = diag(Phat);
end


figure;
subplot(2,1,1);
plot(time,xhat(1,:)-xtrue(1,:),'k'); hold on;
plot(time,3*sqrt(Ps(1,:)),'r'); plot(time,-3*sqrt(Ps(1,:)),'r');
subplot(2,1,2);
plot(time,xhat(2,:)-xtrue(2,:),'k'); hold on;
plot(time,3*sqrt(Ps(2,:)),'r'); plot(time,-3*sqrt(Ps(2,:)),'r');

figure;
plot(xtrue(1,:),xtrue(2,:),'b'); hold on;
plot(xhat(1,:),xhat(2,:),'r'); 
axis equal



function dx = VanderPol(t,x,mu)
    dx(1,1) = x(2);
    dx(2,1) = mu*(1-x(1)^2)*x(2) - x(1);
end

function xdot = VanderPol_EKF(t,x,mu)

    % state
    state = x(1:2);
    xdot(1,1) = state(2);
    xdot(2,1) = mu*(1-state(1)^2)*state(2) - state(1);

    % covariance
    F = [     0,                          1;
        -2*mu*state(1)*state(2)-1,   mu-mu*state(1)^2];
    covariance = x(3:6);
    P = reshape(covariance,2,2);
    Pdot = F*P + P*F';
    xdot(3:6,1) = reshape(Pdot,4,1);
end
























































