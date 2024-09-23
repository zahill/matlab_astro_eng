% PROBLEM 3
mu = 398600;

x0 = 7000;
y0 = 0;
z0 = 0;

xdot0 = 0;
ydot0 = 9;
zdot0 = 0;

phi0 = eye(6);

P0 = diag([10^2, 10^2, 10^2, (1e-3)^2, (1e-3)^2, (1e-3)^2]);

% PART A
n = 1000;
dx = 10*randn(3,n);
figure()
hold on
theta = 0:pi/50:2*pi;
plot(10*cos(theta) + x0, 10*sin(theta) + y0, 'LineWidth', 2)
plot(20*cos(theta) + x0, 20*sin(theta) + y0, 'LineWidth', 2)
plot(30*cos(theta) + x0, 30*sin(theta) + y0, 'LineWidth', 2)
for i=1:n
    scatter(x0 + dx(1,i), y0 + dx(2,i))
end
title('Random Sample Distribution')
subtitle('Position')
xlabel('x [km]')
ylabel('y [km]')
hold off

dxdot = 1e-3*randn(3,n);
figure()
hold on
theta = 0:pi/50:2*pi;
plot(1e-3*cos(theta) + xdot0, 1e-3*sin(theta) + ydot0, 'LineWidth', 2)
plot(2e-3*cos(theta) + xdot0, 2e-3*sin(theta) + ydot0, 'LineWidth', 2)
plot(3e-3*cos(theta) + xdot0, 3e-3*sin(theta) + ydot0, 'LineWidth', 2)
for i=1:n
    scatter(xdot0 + dxdot(1,i), ydot0 + dxdot(2,i))
end
title('Random Sample Distribution')
subtitle('Velocity')
xlabel('vx [km/s]')
ylabel('vy [km/s]')
hold off

% PART B
E = (xdot0^2 + ydot0^2 + zdot0^2)/2 - mu/sqrt(x0^2 + y0^2 + z0^2);
a = -mu/2/E;
T = 2*pi/sqrt(mu)*a^(3/2);
tspan = [0,T/4,T/2,3*T/4,T];
[~,y] = ode45(@(t,X) odewithstm(t,X,mu),tspan,[x0,y0,z0,xdot0,ydot0, ...
    xdot0,reshape(phi0,1,[])]);
fprintf('| %10s | %10s | %10s | %10s | %10s | %10s | %10s |\n', ...
    'Time','Sigma x','Sigma y','Sigma z','Sigma vx','Sigma vy','Sigma vz')
for i=1:5
    t = tspan(i);
    phi = reshape(y(i,7:end),6,6);
    P = phi*P0*transpose(phi);
    sigma_x = sqrt(P(1,1));
    sigma_y = sqrt(P(2,2));
    sigma_z = sqrt(P(3,3));
    sigma_vx = sqrt(P(4,4));
    sigma_vy = sqrt(P(5,5));
    sigma_vz = sqrt(P(6,6));
    fprintf('| %10.4f | %10.4f | %10.4f | %10.4f | %10.4f | %10.4f | %10.4f |\n' ...
        ,t,sigma_x,sigma_y,sigma_z,sigma_vx,sigma_vy,sigma_vz)
end

% PART C & E
for j = 1:5
    t = tspan(j);
    figure()
    hold on
    for i=1:n
        X0 = [x0 + dx(1,i); y0 + dx(2,i); z0 + dx(3,i); xdot0 + dxdot(1,i);
            ydot0 + dxdot(2,i); zdot0 + dxdot(3,i)];
        phi = reshape(y(j,7:end),6,6);
        X = phi*X0;
        scatter(X(1), X(2))
    end
    title('Linear Map Position')
    subtitle(sprintf('t=%.4f',t))
    xlabel('x [km]')
    ylabel('y [km]')
    hold off
end

for j = 1:5
    t = tspan(j);
    figure()
    hold on
    for i=1:n
        X0 = [x0 + dx(1,i); y0 + dx(2,i); z0 + dx(3,i); xdot0 + dxdot(1,i);
            ydot0 + dxdot(2,i); zdot0 + dxdot(3,i)];
        phi = reshape(y(j,7:end),6,6);
        X = phi*X0;
        scatter(X(1), X(2))
    end
    title('Linear Map Velocity')
    subtitle(sprintf('t=%.4f',t))
    xlabel('vx [km/s]')
    ylabel('vy [km/s]')
    hold off
end

% PART D & E
for j = 1:5
    t = tspan(j);
    figure()
    hold on
    for i=1:n
        X0 = [x0 + dx(1,i); y0 + dx(2,i); z0 + dx(3,i); xdot0 + dxdot(1,i);
            ydot0 + dxdot(2,i); zdot0 + dxdot(3,i)];
        if t > 0
            [~,y] = ode45(@(t,X) odefunc(t,X,mu),[0 t],X0);
            X = y(end,:);
        else
            X = X0;
        end
        scatter(X(1), X(2))
    end
    title('Non-Linear Map Position')
    subtitle(sprintf('t=%.4f',t))
    xlabel('x [km]')
    ylabel('y [km]')
    hold off
end

for j = 1:5
    t = tspan(j);
    figure()
    hold on
    for i=1:n
        X0 = [x0 + dx(1,i); y0 + dx(2,i); z0 + dx(3,i); xdot0 + dxdot(1,i);
            ydot0 + dxdot(2,i); zdot0 + dxdot(3,i)];
        if t > 0
            [~,y] = ode45(@(t,X) odefunc(t,X,mu),[0 t],X0);
            X = y(end,:);
        else
            X = X0;
        end
        scatter(X(4), X(5))
    end
    title('Non-Linear Map Velocity')
    subtitle(sprintf('t=%.4f',t))
    xlabel('vx [km/s]')
    ylabel('vy [km/s]')
    hold off
end

function dX = odewithstm(~, X, mu)
x = X(1);
y = X(2);
z = X(3);
xdot = X(4);
ydot = X(5);
zdot = X(6);
phi = reshape(X(7:end),6,6);

r = sqrt(x^2 + y^2 + z^2);

A = zeros(6);
A(1,4) = 1;
A(2,5) = 1;
A(3,6) = 1;
A(4,1) = 3*mu*x^2/r^5 - mu/r^3;
A(4,2) = 3*mu*x*y/r^5;
A(4,3) = 3*mu*x*z/r^5;
A(5,1) = 3*mu*x*y/r^5;
A(5,2) = 3*mu*y^2/r^5 - mu/r^3;
A(5,3) = 3*mu*y*z/r^5;
A(6,1) = 3*mu*x*z/r^5;
A(6,2) = 3*mu*y*z/r^5;
A(6,3) = 3*mu*z^2/r^5 - mu/r^3;

dX = zeros(6+6*6,1);
dX(1) = xdot;
dX(2) = ydot;
dX(3) = zdot;
dX(4) = -mu/r^3*x;
dX(5) = -mu/r^3*y;
dX(6) = -mu/r^3*z;
dX(7:end) = reshape(A*phi,[],1);
end

function dX = odefunc(~, X, mu)
x = X(1);
y = X(2);
z = X(3);
xdot = X(4);
ydot = X(5);
zdot = X(6);

r = sqrt(x^2 + y^2 + z^2);

dX = zeros(6,1);
dX(1) = xdot;
dX(2) = ydot;
dX(3) = zdot;
dX(4) = -mu/r^3*x;
dX(5) = -mu/r^3*y;
dX(6) = -mu/r^3*z;
end