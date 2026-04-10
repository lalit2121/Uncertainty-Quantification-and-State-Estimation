clear; close all; clc;
m0 = [0.4 0.3]';
P0 = 0.2^2*eye(2);
S0 = sqrtm(P0);
mu = 3.1;
R = 0.1^2;
SR = sqrt(R);
H = [1 0];
time = linspace(0,15,151);
options = odeset('abstol',1e-8,'reltol',1e-8);


xtrue = zeros(2,length(time));
xtrue(:,1) = m0 + S0*randn(2,1);
xhat = zeros(2,length(time));
xhat(:,1) = m0;
xest = zeros(2,length(time));
xest(:,1) = m0;
Phat = P0;

Pest = zeros(2,length(time));
Pest(:,1) = diag(P0);

par = 0.2:0.2:4;
n = length(par);
w = 1/n*ones(1,n);

means = cell(1,n);
covs = cell(1,n);

for j = 1:n
    means{j} = xhat;
    covs{j} = Phat;
end

for i = 2:length(time)

    [~,x] = ode45(@(t,x) VanderPol(t,x,mu),[time(i-1) time(i)],xtrue(:,i-1),options);
    xtrue(:,i) = x(end,:)';
    ytrue = H*xtrue(:,i) + SR*randn(1);
    
    for j = 1:n
        
        xhat = means{j};
        Phat = covs{j};
        
        state = [xhat(:,i-1); reshape(Phat,4,1)];
        [~,state] = ode45(@(t,x) VanderPol_EKF(t,x,par(j)),[time(i-1) time(i)],state,options);
        state = state(end,:)';
        xhat(:,i) = state(1:2);
        Phat = reshape(state(3:6),2,2);
        
        yhat = H*xhat(:,i);
        Pxy = Phat*H';
        Pyy = H*Phat*H'+R;
        K = Pxy/Pyy;
        xhat(:,i) = xhat(:,i) + K*(ytrue-yhat);
        Phat = Phat - K*Pyy*K';

        w(j) = w(j) * (1/(2*pi*sqrt(Pyy))) * exp(-1/2*((ytrue - yhat)*(ytrue - yhat)/Pyy));

        means{j} = xhat;
        covs{j} = Phat;
       
    end
    w = w/sum(w);

    xest(:,i) = zeros(2,1);
    P = zeros(2,2);
    for j = 1:n
        xest(:,i) = xest(:,i) + w(j) * means{j}(:,i);
    end
    for j = 1:n
        P = P + w(j)*(means{j}(:,i)*means{j}(:,i)' + covs{j});
    end
    P = P - xest(:,i)*xest(:,i)';
    Pest(:,i) = diag(P);

    figure(1);
    set(figure(1),'Position', [100, 100, 800, 600]);
    bar(w); ylim([0 1]); title(sprintf('Step %d', i-1));
    set(gca, 'XTick', 1:n);
    set(gca, 'XTickLabel', string(par));

    pause(1/50);

end


figure;
subplot(2,1,1);
plot(time,xest(1,:)-xtrue(1,:),'k'); hold on;
plot(time,3*sqrt(Pest(1,:)),'r'); plot(time,-3*sqrt(Pest(1,:)),'r');
subplot(2,1,2);
plot(time,xest(2,:)-xtrue(2,:),'k'); hold on;
plot(time,3*sqrt(Pest(2,:)),'r'); plot(time,-3*sqrt(Pest(2,:)),'r');

figure;
plot(xtrue(1,:),xtrue(2,:),'b'); hold on;
plot(xest(1,:),xest(2,:),'r'); 
axis equal


%% Auxiliarly Functions
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
























































