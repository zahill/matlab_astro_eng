clc
clear all
close all

mu = 398600;
R = 6400;
J2 = 1e-3;

r0 = [7000; 0; 0];
v0 = [0; 4.9; 5.8];

% PROBLEM 2a
[a, e, i, Omega, omega, M] = state2coes([r0; v0],mu);
fprintf('a: %.4f km\n',a)
fprintf('e: %.4f\n',e)
fprintf('i: %.4f(%.4f) rad(deg)\n',i,rad2deg(i))
fprintf('RAAN: %.4f(%.4f) rad(deg)\n',Omega,rad2deg(Omega))
fprintf('omega: %.4f(%.4f) rad(deg)\n',omega,rad2deg(omega))
fprintf('M: %.4f(%.4f) rad(deg)\n',M,rad2deg(M))

% PROBLEM 2b
T = 2*pi/sqrt(mu)*a^(3/2);
tspan = linspace(0,10*T,1000);
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t,X1] = ode45(@(t,X) twobody(t, mu, X), tspan, [r0; v0], opts);

[x y z] = sphere;
x = x*R;
y = y*R;
z = z*R;

fig = figure();
plot3(X1(:,1),X1(:,2),X1(:,3),'LineWidth',2,'Color','#0072BD')
hold on
surf(x,y,z,'FaceColor','#77AC30','EdgeColor','black')
title('3D Orbit')
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal
hold off

fig = figure();
plot(X1(:,1),X1(:,2),'LineWidth',2,'Color','#0072BD')
hold on
scatter(0,0,'filled')
title('X-Y Plane')
xlabel('X')
ylabel('Y')
axis equal
hold off

% PROBLEM 2c
[t,X2] = ode45(@(t,X) J2effects(t, mu, R, J2, X), tspan, [r0; v0], opts);

fig = figure();
plot(t,X2(:,1)-X1(:,1),'LineWidth',2)
hold on
plot(t,X2(:,2)-X1(:,2),'LineWidth',2)
plot(t,X2(:,3)-X1(:,3),'LineWidth',2)
title('Position Deviation')
xlabel('Time [s]')
ylabel('\delta r [km]')
legend('\delta r_x','\delta r_y','\delta r_z')
hold off

fig = figure();
plot(t,X2(:,4)-X1(:,4),'LineWidth',2)
hold on
plot(t,X2(:,5)-X1(:,5),'LineWidth',2)
plot(t,X2(:,6)-X1(:,6),'LineWidth',2)
title('Velocity Deviation')
xlabel('Time [s]')
ylabel('\delta v [km/s]')
legend('\delta v_x','\delta v_y','\delta v_z')
hold off

coes = zeros(length(X2),6);
for k = 1:length(X2)
    [a, e, i, Omega, omega, M] = state2coes(X2(k,:),mu);
    coes(k,1) = a;
    coes(k,2) = e;
    coes(k,3) = i;
    coes(k,4) = Omega;
    coes(k,5) = omega;
    coes(k,6) = M;
end

% PROBLEM 3
[a, e, i, Omega, omega, M] = state2coes([r0; v0],mu);
[t,X3] = ode45(@(t,X) lagrange_eqs(t, mu, R, J2, X), tspan, [a, e, i, Omega, omega, M], opts);
fig = figure();
subplot(5,1,1)
plot(t,X3(:,1))
hold on
plot(t,coes(:,1))
hold off
title('Semi-Major Axis')
legend('P3','P2')

subplot(5,1,2)
plot(t,X3(:,2))
hold on
plot(t,coes(:,2))
hold off
title('Eccentricity')
legend('P3','P2')

subplot(5,1,3)
plot(t,X3(:,3))
hold on
plot(t,coes(:,3))
hold off
title('Inclination')
legend('P3','P2')

subplot(5,1,4)
plot(t,X3(:,4))
hold on
plot(t,coes(:,4))
hold off
title('RAAN')
legend('P3','P2')

subplot(5,1,5)
plot(t,X3(:,5))
hold on
plot(t,coes(:,5))
hold off
title('Argument of Periapsis')
legend('P3','P2')

function [a, e, i, Omega, omega, M] = state2coes(state,mu)
    r = transpose(state(1:3));
    v = transpose(state(4:6));

    I = [1, 0, 0];
    J = [0, 1, 0];
    K = [0, 0, 1];

    a = -mu/2*(norm(v)^2/2 - mu/norm(r))^(-1);
    h = cross(r,v);
    e = cross(v,h)/mu - r/norm(r);
    i = acos(h(3)/norm(h));
    n = cross(K,h);
    if n(2) >= 0
        Omega = acos(n(1)/norm(n));
    else
        Omega = -acos(n(1)/norm(n));
    end
    if e(3) >= 0
        omega = acos(dot(e,n)/(norm(e)*norm(n)));
    else
        omega = -acos(dot(e,n)/(norm(e)*norm(n)));
    end
    if dot(r,v) >= 0
        nu = acos(dot(e,r)/(norm(e)*norm(r)));
    else
        nu = -acos(dot(e,r)/(norm(e)*norm(r)));
    end
    E = acos((norm(e) + cos(nu))/(1 + norm(e)*cos(nu)));
    M = E - norm(e)*sin(E);
    e = norm(e);
end

function dX = twobody(t, mu, X)
    r = X(1:3);
    v = X(4:6);

    dX = zeros(6,1);
    dX(1:3) = v;
    dX(4:6) = -mu/norm(r)^3*r;
end

function dX = J2effects(t, mu, R, J2, X)
    r = X(1:3);
    v = X(4:6);

    dX = zeros(6,1);
    dX(1:3) = v;
    dX(4:6) = -mu/norm(r)^3*r - 3*mu/(2*norm(r)^5)*R^2*J2*([1, 0, 0; 0, 1, 0; 0, 0, 3] - 5*(r(3)/norm(r))^2*eye(3))*r;
end

function dX = lagrange_eqs(t, mu, R, J2, X)
    a = X(1);
    e = X(2);
    i = X(3);
    Omega = X(4);
    omega = X(5);
    M = X(6);
    n = sqrt(mu/a^3);

    dRdi = -3*mu*R^2*J2/(2*a^3*(1 - e^2)^(3/2))*sin(i)*cos(i);
    dRda = -3*mu*R^2*J2/(2*a^3*(1 - e^2)^(3/2))*(1 - 3/2*sin(i)^2);
    dRde = 3*mu*R^2*J2/(2*a^3*(1 - e^2)^(5/2))*e*(1 - 3/2*sin(i)^2);
    
    dX = zeros(6,1);
    dX(1) = 0;
    dX(2) = 0;
    dX(3) = 0;
    dX(4) = 1/(n*a^2*sqrt(1 - e^2)*sin(i))*dRdi;
    dX(5) = sqrt(1 - e^2)/(n*a^2*e)*dRde - cot(i)/(n*a^2*sqrt(1 - e^2))*dRdi;
    dX(6) = n - (1 - e^2)/(n*a^2*e)*dRde - 2/(n*a)*dRda;
end
