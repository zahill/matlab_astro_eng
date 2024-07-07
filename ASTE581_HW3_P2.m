clear all
close all

%  ARBITRARY TRAJECTORIES
mu = 0.1;
CPR3BP(mu, [0.5, 0, 0, 0.5], 2*pi);
CPR3BP(mu, [0.25, 1, 0.75, 0.5], 2*pi);
CPR3BP(mu, [0, 0.75, 1.625, 0], 2*pi);
CPR3BP(mu, [-0.652, -0.652, 0.5, 0.25], 2*pi);
CPR3BP(mu, [1, -1, -0.5, -0.5], 2*pi);

% LAGRANGE POINTS 4 & 5
mu = 0.01;
CPR3BP(mu, [1/2 - mu + 0.01, sqrt(3)/2 + 0.01, 0, 0], 100, "L4 Stability, mu=0.01");
CPR3BP(mu, [1/2 - mu + 0.01, -sqrt(3)/2 + 0.01, 0, 0], 100, "L5 Stability, mu=0.01");

mu = 0.1;
CPR3BP(mu, [1/2 - mu + 0.01, sqrt(3)/2 + 0.01, 0, 0], 100, "L4 Stability, mu=0.1");
CPR3BP(mu, [1/2 - mu + 0.01, -sqrt(3)/2 + 0.01, 0, 0], 100, "L5 Stability, mu=0.1");

% LAGRANGE POINTS 1 & 2
mu = 0.01;
CPR3BP(mu, [1 - (mu/3)^(1/3) + 0.01, 0 + 0.01, -0.25, 0.25], 10, "L1 Stability, mu=0.01");
CPR3BP(mu, [1 + (mu/3)^(1/3) + 0.01, 0 + 0.01, -0.25, 0.25], 10, "L2 Stability, mu=0.01");

function [X,J] = CPR3BP(mu,X0,t, txt)
tspan = linspace(0,t,1000);
[t,X] = ode45(@(t,X) ode_func(t,mu,X),tspan,X0);

J = zeros(1,length(X));
for i = 1:length(J)
    x = X(i,1);
    y = X(i,2);
    xdot = X(i,3);
    ydot = X(i,4);

    r1 = sqrt((x + mu)^2 + y^2);
    r2 = sqrt((x - 1 + mu)^2 + y^2);
    V = (x^2 + y^2)/2 + (1 - mu)/r1 + mu/r2;
    J(i) = (xdot^2 + ydot^2)/2 - V;
end

if exist('txt',"var") ~= 1
    txt = sprintf('mu = %.2f, X0 = [%.2f,%.2f,%.2f,%.2f]',mu,X0);
end

fig = figure();
plot(t,J)
axis equal
box on
title('Jacobian Integral')
subtitle(txt)
xlabel('Time')
ylabel('Jacobian Integral')
set(gcf,'Color','w');
set(gca,'FontSize',12,'FontName','Courier New','Color','w');

fig = figure();
hold on
scatter(-mu, 0, "filled")
scatter(1 - mu, 0, "filled")
scatter(1 - (mu/3)^(1/3), 0, "filled")
scatter(1 + (mu/3)^(1/3), 0, "filled")
scatter(-1 - (sqrt(2) - 1)/3*mu, 0, "filled")
scatter(1/2 - mu, sqrt(3)/2, "filled")
scatter(1/2 - mu, -sqrt(3)/2, "filled")
plot(X(:,1),X(:,2))
axis equal
box on
title('X-Y Plane')
subtitle(txt)
xlabel('X-Position [-]')
ylabel('Y-Position [-]')
legend('P1','P2','L1','L2','L3','L4','L5','Trajectory','Location','eastoutside')
set(gcf,'Color','w');
set(gca,'FontSize',12,'FontName','Courier New','Color','w');
hold off
end

function dX = ode_func(~,mu,X)
x = X(1);
y = X(2);
xdot = X(3);
ydot = X(4);

r1 = sqrt((x + mu)^2 + y^2);
r2 = sqrt((x - 1 + mu)^2 + y^2);

dX = zeros(length(X),1);
dX(1) = xdot;
dX(2) = ydot;
dX(3) = 2*ydot + x - (1 - mu)*(x + mu)/r1^3 - mu*(x + mu - 1)/r2^3;
dX(4) = -2*xdot + y*(1 - (1 - mu)/r1^3 - mu/r2^3);
end