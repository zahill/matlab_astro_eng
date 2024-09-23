clc;
clf;
close all;
clear all;

% Satellite initial values
r0 = [757700.0; 5222607.0; 4851500.0]; %m
v0 = [2213.21; 4678.34; -5371.30]; %m/s
phi0 = eye(18);

% Constants
R = 6378136.3; %m
mu = 3.986004415e14; %m^3/s^2
J2 = 1.082626925638815e-3; %-
CD = 2; %-
X1 = [-5127510.0; -3794160.0; 0.0]; %m
X2 = [3860910.0; 3238490.0; 3898094.0]; %m
X3 = [549505.0; -1380872.0; 6182197.0]; %m

% Initialize the vectors and matricies
Xstar0 = [r0; v0; mu; J2; CD; X1; X2; X3];
sigma_rho = 0.01; %m
sigma_rhodot = 0.001; %m
W = [1/sigma_rho^2 0; 0 1/sigma_rhodot^2];
x0bar = zeros(18,1);
P0 = diag([ones(1,6)*1e6, 1e20, ones(1,2)*1e6, ones(1,3)*1e-10, ones(1,6)*1e6]);

% Import observation data
Y = readmatrix("obs_data.xlsx");

% Set the integrator absolute and relative tolerances
opts = odeset('AbsTol',1e-12,'RelTol',1e-12);

% Integrate and plot the nominal orbit
[t,X] = ode89(@(t,X) ode_func(t, X),Y(:,1),[Xstar0; reshape(phi0,[],1)],opts);
plot_orbit("Nominal Orbit",X);

% Run the batch processor using all observation data
data1 = batch_processor("Batch Processor, All Data",Xstar0,Y,x0bar,P0,W,[1 2],1e-3,opts);

% Run the batch processor using only range data
data2 = batch_processor("Batch Processor, Range Data",Xstar0,Y,x0bar,P0,W,[1],1e-3,opts);

% Run the batch processor using only the range-rate data
data3 = batch_processor("Batch Processor, Range-Rate Data",Xstar0,Y,x0bar,P0,W,[2],1e-3,opts);

% Run the Kalman Filter
sigma = kalman_filter("Kalman Filter",Xstar0,Y,x0bar,P0,inv(W),opts);

% Write iteration data to a file
writematrix(data1,"iteration_data.xlsx","Sheet","Batch Processor All Data");
writematrix(data2,"iteration_data.xlsx","Sheet","Batch Processor Range Data");
writematrix(data3,"iteration_data.xlsx","Sheet","Batch Processor Range-Rate Data");
writematrix(sigma,"iteration_data.xlsx","Sheet","Kalman Filter");
function dX = ode_func(~, X)
% ODE_FUNC Calculte the derivative of the state matrix for the current 
% timestep.
% 
% INPUTS:
%     - X  : The value of the state matrix for the previous time step.
%
% OUTPUTS:
%     - dX : The value of the derivative of the state matrix.

% Constants
R = 6378136.3; %m
S = 3; %m^2
m = 970; %kg
rho0 = 3.614e-13; %kg/m^3
r0 = 700000.0 + R; %m
H = 88667; %m
thetadot = 7.2921158543e-5; %rad/s

% Define the previous state variables
x = X(1);
y = X(2);
z = X(3);
xdot = X(4);
ydot = X(5);
zdot = X(6);
mu = X(7);
J2 = X(8);
CD = X(9);
phi = reshape(X(19:end),18,18);

% Calculate additional needed variables
r = sqrt(x^2 + y^2 + z^2);
rho = rho0*exp(-(r - r0)/H);
v_a = norm([xdot;ydot;zdot] - cross([0;0;thetadot],[x;y;z]));
A = Amatrix(X);

% Define the derivative of the state vector
dX = zeros(18+18*18,1);
dX(1) = xdot;
dX(2) = ydot;
dX(3) = zdot;
dX(4) = -mu*x*(1/r^3 + 3/2*J2*R^2/r^5*(1 - 5*(z/r)^2)) - 1/2*CD*S/m*rho*v_a*(xdot + thetadot*y);
dX(5) = -mu*y*(1/r^3 + 3/2*J2*R^2/r^5*(1 - 5*(z/r)^2)) - 1/2*CD*S/m*rho*v_a*(ydot - thetadot*x);
dX(6) = -mu*z*(1/r^3 + 3/2*J2*R^2/r^5*(3 - 5*(z/r)^2)) - 1/2*CD*S/m*rho*v_a*zdot;
dX(19:end) = reshape(A*phi,[],1);
end

function A = Amatrix(X)
% AMATRIX Calculate the A matrix for the current timestep.
% 
% INPUTS:
%     - X : The current state matrix.
%
% OUTPUTS:
%     - A : The A matrix for the current timestep.

% Constants
R = 6378136.3; %m
S = 3; %m^2
m = 970; %kg
rho0 = 3.614e-13; %kg/m^3
r0 = 700000.0 + R; %m
H = 88667; %m
thetadot = 7.2921158543e-5; %rad/s

% Define the current state variables
x = X(1);
y = X(2);
z = X(3);
xdot = X(4);
ydot = X(5);
zdot = X(6);
mu = X(7);
J2 = X(8);
CD = X(9);

% Calculate additional needed variables
r = sqrt(x^2 + y^2 + z^2);
rho = rho0*exp(-(r - r0)/H);
v_a = norm([xdot;ydot;zdot] - cross([0;0;thetadot],[x;y;z]));

% Calculate the A matrix
A = zeros(18,18);
A(1,4) = 1;
A(2,5) = 1;
A(3,6) = 1;
A(4,1) = -mu/r^3*(1 - 3/2*J2*(R/r)^2*(5*(z/r)^2 - 1)) + 3*mu*x^2/r^5*(1 - 5/2*J2*(R/r)^2*(7*(z/r)^2 - 1)) + 1/2*CD*S/m*rho*(v_a*x*(xdot + thetadot*y)/(r*H) - (-thetadot*ydot + thetadot^2*x)*(xdot + thetadot*y)/v_a);
A(4,2) = 3*mu*x*y/r^5*(1 - 5/2*J2*(R/r)^2*(7*(z/r)^2 - 1)) + 1/2*CD*S/m*rho*(v_a*y*(xdot + thetadot*y)/(r*H) - (thetadot*xdot + thetadot^2*y)*(xdot + thetadot*y)/v_a - v_a*thetadot);
A(4,3) = 3*mu*x*z/r^5*(1 - 5/2*J2*(R/r)^2*(7*(z/r)^2 - 3)) + 1/2*CD*S/m*rho*v_a*z*(xdot + thetadot*y)/(r*H);
A(4,4) = -1/2*CD*S/m*rho*((xdot + thetadot*y)^2/v_a + v_a);
A(4,5) = -1/2*CD*S/m*rho*(ydot - thetadot*x)*(xdot + thetadot*y)/v_a;
A(4,6) = -1/2*CD*S/m*rho*zdot*(xdot + thetadot*y)/v_a;
A(4,7) = -x/r^3*(1 - 3/2*J2*(R/r)^2*(5*(z/r)^2 - 1));
A(4,8) = 3/2*mu*x/r^3*(R/r)^2*(5*(z/r)^2 - 1);
A(4,9) = -1/2*S/m*rho*v_a*(xdot + thetadot*y);
A(5,1) = 3*mu*x*y/r^5*(1 - 5/2*J2*(R/r)^2*(7*(z/r)^2 - 1)) + 1/2*CD*S/m*rho*(v_a*(ydot - thetadot*x)/(r*H) - (thetadot^2*x - thetadot*ydot)*(ydot - thetadot*x)/v_a + v_a*thetadot);
A(5,2) = -mu/r^3*(1 - 3/2*J2*(R/r)^2*(5*(z/r)^2 - 1)) + 3*mu*y^2/r^5*(1 - 5/2*J2*(R/r)^2*(7*(z/r)^2 - 1)) + 1/2*CD*S/m*rho*(v_a*y*(ydot - thetadot*x)/(r*H) - (thetadot*xdot + thetadot^2*y)*(ydot - thetadot*x)/v_a);
A(5,3) = 3*mu*y*z/r^5*(1 - 5/2*J2*(R/r)^2*(7*(z/r)^2 - 3)) + 1/2*CD*S/m*rho*v_a*z*(ydot - thetadot*x)/(r*H);
A(5,4) = -1/2*CD*S/m*rho*(ydot - thetadot*x)*(xdot + thetadot*y)/v_a;
A(5,5) = -1/2*CD*S/m*rho*((ydot - thetadot*x)^2/v_a + v_a);
A(5,6) = -1/2*CD*S/m*rho*zdot*(ydot - thetadot*x)/v_a;
A(5,7) = -y/r^3*(1 - 3/2*J2*(R/r)^2*(5*(z/r)^2 - 1));
A(5,8) = 3/2*mu*y/r^3*(R/r)^2*(5*(z/r)^2 - 1);
A(5,9) = -1/2*S/m*rho*v_a*(ydot - thetadot*x);
A(6,1) = 3*mu*x*z/r^5*(1 - 5/2*J2*(R/r)^2*(7*(z/r)^2 - 3)) + 1/2*CD*S/m*rho*(v_a*zdot*x/(r*H) - zdot*(thetadot^2*x - thetadot*ydot)/v_a);
A(6,2) = 3*mu*y*z/r^5*(1 - 5/2*J2*(R/r)^2*(7*(z/r)^2 - 3)) + 1/2*CD*S/m*rho*(v_a*zdot*y/(r*H) - zdot*(thetadot*xdot + thetadot^2*y)/v_a);
A(6,3) = -mu/r^3*(1 - 3/2*J2*(R/r)^2*(5*(z/r)^2 - 3)) + 3*mu*z^2/r^5*(1 - 5/2*J2*(R/r)^2*(7*(z/r)^2 - 5)) + 1/2*CD*S/m*rho*v_a*z*zdot/(r*H);
A(6,4) = -1/2*CD*S/m*rho*zdot*(xdot + thetadot*y)/v_a;
A(6,5) = -1/2*CD*S/m*rho*zdot*(ydot - thetadot*x)/v_a;
A(6,6) = -1/2*CD*S/m*rho*(zdot^2/v_a + v_a);
A(6,7) = -z/r^3*(1 - 3/2*J2*(R/r)^2*(5*(z/r)^2 - 3));
A(6,8) = 3/2*mu*z/r^3*(R/r)^2*(5*(z/r)^2 - 3);
A(6,9) = -1/2*S/m*rho*v_a*zdot;
end

function H_tilde = H_tilde(r,v,rs,rho,rhodot,theta,thetadot,si)
% H_TILDE Calculate the H_tilde matrix for the current timestep.
% 
% INPUTS:
%     - r        : The radius vector of the satellite.
%     - v        : The velocity vector of the satellite.
%     - rs       : The radius vector to the current station.
%     - rho      : The calculated range of the satellite from the current
%     station.
%     - rhodot   : The calculated range-rate of the satellite.
%     - theta    : The angle between the ECI and ECEF coordinate systems.
%     - thetadot : The rotational velocity of the earth.
%     - si       : The indices for the current station in the state vector.
%
% OUTPUTS:
%     - H_tilde  : The H_tilde matrix for the current timestep.

% Define the current state variables
x = r(1);
y = r(2);
z = r(3);
xdot = v(1);
ydot = v(2);
zdot = v(3);
xs = rs(1);
ys = rs(2);
zs = rs(3);

% Calculate the H_tilde matrix
H_tilde = zeros(2,18);
H_tilde(1,1) = (x - xs*cos(theta) + ys*sin(theta))/rho;
H_tilde(1,2) = (y - ys*cos(theta) - xs*sin(theta))/rho;
H_tilde(1,3) = (z - zs)/rho;
H_tilde(1,si(1)) = (xs - x*cos(theta) - y*sin(theta))/rho;
H_tilde(1,si(2)) = (ys - y*cos(theta) + x*sin(theta))/rho;
H_tilde(1,si(3)) = (zs - z)/rho;
H_tilde(2,1) = (xdot + thetadot*xs*sin(theta) + thetadot*ys*cos(theta))/rho - rhodot*(x - xs*cos(theta) + ys*sin(theta))/rho^2;
H_tilde(2,2) = (ydot + thetadot*ys*sin(theta) - thetadot*xs*cos(theta))/rho - rhodot*(y - ys*cos(theta) - xs*sin(theta))/rho^2;
H_tilde(2,3) = zdot/rho - rhodot*(z - zs)/rho^2;
H_tilde(2,4) = (x - xs*cos(theta) + ys*sin(theta))/rho;
H_tilde(2,5) = (y - ys*cos(theta) - xs*sin(theta))/rho;
H_tilde(2,6) = (z - zs)/rho;
H_tilde(2,si(1)) = (-xdot*cos(theta) + thetadot*x*sin(theta) - ydot*sin(theta) - thetadot*y*cos(theta))/rho - rhodot*(xs - x*cos(theta) - y*sin(theta))/rho^2;
H_tilde(2,si(2)) = (-ydot*cos(theta) + thetadot*y*sin(theta) + xdot*sin(theta) + thetadot*x*cos(theta))/rho - rhodot*(ys - y*cos(theta) + x*sin(theta))/rho^2;
H_tilde(2,si(3)) = -zdot/rho - rhodot*(zs - z)/rho^2;
end

function H = Hmatrix(r,v,rs,rho,rhodot,theta,thetadot,si,phi)
% HMATRIX Calculate the H matrix for the current timestep.
% 
% INPUTS:
%     - r        : The radius vector of the satellite.
%     - v        : The velocity vector of the satellite.
%     - rs       : The radius vector to the current station.
%     - rho      : The calculated range of the satellite from the current
%     station.
%     - rhodot   : The calculated range-rate of the satellite.
%     - theta    : The angle between the ECI and ECEF coordinate systems.
%     - thetadot : The rotational velocity of the earth.
%     - si       : The indices for the current station in the state vector.
%     - phi      : The state-transition matrix.
%
% OUTPUTS:
%     - H        : The H matrix for the current timestep.

H = H_tilde(r,v,rs,rho,rhodot,theta,thetadot,si)*phi;
end

function [x, P] = cholesky_decomp(M,N)
% CHOLESKY_DECOMP Solve the problem Mx=N using the Cholesky decomposition
% method.
% 
% INPUTS:
%     - M : An nxn matrix
%     - N : An nx1 vector
%
% OUTPUTS:
%     - x : The solution to the given problem
%     - P : The covariance matrix for the given problem

% Determine the height of the M matrix
n = height(M);

% Compute the upper-triangular matrix R
R = zeros(size(M));
for i = 1:n
    sum = 0;
    for k = 1:(i-1)
        sum = sum + R(k,i)^2;
    end
    R(i,i) = sqrt(M(i,i) - sum);
    for j = (i+1):n
        sum = 0;
        for k = 1:(i-1)
            sum = sum + R(k,i)*R(k,j);
        end
        R(i,j) = (M(i,j) - sum)/R(i,i);
    end
end

% Compute the Z vector
Z = zeros(n,1);
for i = 1:n
    sum = 0;
    for j = 1:(i-1)
        sum = sum + R(j,i)*Z(j);
    end
    Z(i) = (N(i) - sum)/R(i,i);
end

% Solve for the x vector
x = zeros(n,1);
for i = n:-1:1
    sum = 0;
    for j = i+1:n
        sum = sum + R(i,j)*x(j);
    end
    x(i) = (Z(i) - sum)/R(i,i);
end

% Compute the S matrix
S = zeros(size(M));
for i = 1:n
    S(i,i) = 1/R(i,i);
end
for i = 1:n
    for j = i+1:n
        sum = 0;
        for k = i:(j-1)
            sum = sum + R(k,j)*S(i,k);
        end
        S(i,j) = -S(j,j)*sum;
    end
end

% Solve for the covariance matrix
P = S*S';
end

function sigma = std_dev(P)
% STD_DEV Calculate the standard deviations from the given covariance
% matrix.
% 
% INPUTS:
%     - P     : The covariance matrix
%
% OUTPUTS:
%     - sigma : The calculated standard deviations

% Compute the standard deviation for each variable by taking the squareroot
% of the elements on the diagonal.
sigma = zeros(length(P),1);
for i = 1:length(P)
    sigma(i) = sqrt(P(i,i));
end
end

function plot_orbit(name,X)
% PLOT_ORBIT Plots the orbit for the given state time-series
% 
% INPUTS:
%     - name : The title for the plot
%     - X    : The state time-series to plot

% Constants
R = 6378136.3; %m

% Calculate the average altitude and the orbit inclination
r = sqrt(X(:,1).^2 + X(:,2).^2 + X(:,3).^2);
h_avg = mean(r - R)/1000;
h = cross([X(1,1);X(1,2);X(1,3)],[X(1,4);X(1,5);X(1,6)]);
i = acosd(h(3)/norm(h));

% Create and format the plot
fig = figure();
hold on;
plot3(X(:,1),X(:,2),X(:,3),'b','LineWidth', 2);
[x,y,z] = sphere(100);
x = x*R;
y = y*R;
z = z*R;
surf(x,y,z,'EdgeColor','k','FaceColor','#808080');
axis equal;
grid on;
set(gcf,'Color','w');
set(gca,'Color','w');
set(gca,'FontSize',12,'FontName','Courier New');
title(name);
subtitle(sprintf("Average Altitude: %.2f km & Inclination: %.2f deg",h_avg,i));
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
view(-45,30);
hold off;
end

function data = batch_processor(name,Xstar0,Y,x0bar,P0,W,obs,tol,opts)
% BATCH_PROCESSOR Run the batch processor
% 
% INPUTS:
%     - name   : The name to use for the output data and plots
%     - XStar0 : The initial state vector
%     - Y      : The observation matrix
%     - x0bar  : The a priori x data
%     - P0     : The a priori covariance matrix 
%     - W      : The weight matrix
%     - obs    : A list of which observation data to use
%     - tol    : The required tolerance for the batch processor
%     - opts   : The options for the integrator
%
% OUTPUTS:
%     - data   : A table containing the state vector, sigmas, and RMS
%     values for each iteration

% Constants
thetadot = 7.2921158543e-5; %rad/s

% Save the location of station 1
X1 = Xstar0(10:12);

% Initialize needed variables for the processor
new = 0;
old = 0;
err = 1;
j = 1;

% Run the processor until it has converged to the required tolerance or it
% hits a maximum number of iterations
while err > tol && j<=20
    
    % Initialize this iteration
    old = new;
    M = pinv(P0);
    N = pinv(P0)*x0bar;
    epsilon = zeros(length(Y),2);
    [~,X] = ode89(@(t,X) ode_func(t, X),Y(:,1),[Xstar0; reshape(eye(18),[],1)],opts);
    for i = 1:size(Y,1)

        % Read the next observation
        t = Y(i,1);
        station = Y(i,2);
        rho_obs = Y(i,3);
        rhodot_obs = Y(i,4);

        % Determine the current station
        if station == 101
            i_xs = 10;
            i_ys = 11;
            i_zs = 12;
        elseif station == 337
            i_xs = 13;
            i_ys = 14;
            i_zs = 15;
        else
            i_xs = 16;
            i_ys = 17;
            i_zs = 18;
        end

        % Read the next state values from the integrator
        x = X(i,1);
        y = X(i,2);
        z = X(i,3);
        xdot = X(i,4);
        ydot = X(i,5);
        zdot = X(i,6);
        mu = X(i,7);
        J2 = X(i,8);
        CD = X(i,9);
        xs = X(i,i_xs);
        ys = X(i,i_ys);
        zs = X(i,i_zs);
        phi = reshape(X(i,19:end),18,18);
    
        % Compute needed values
        theta = thetadot*t; %rad
        rho = sqrt(x^2 + y^2 + z^2 + xs^2 + ys^2 + zs^2 - 2*(x*xs + y*ys)*cos(theta) + 2*(x*ys - y*xs)*sin(theta) - 2*z*zs); %m
        rhodot = (x*xdot + y*ydot + z*zdot - (xdot*xs + ydot*ys)*cos(theta) + thetadot*(x*xs + y*ys)*sin(theta) + (xdot*ys - ydot*xs)*sin(theta) + thetadot*(x*ys - y*xs)*cos(theta) - zdot*zs)/rho; %m/s
    
        % Accumulate current observations
        H = Hmatrix([x;y;z],[xdot;ydot;zdot],[xs;ys;zs],rho,rhodot,theta,thetadot,[i_xs;i_ys;i_zs],phi);
        y = [rho_obs; rhodot_obs] - [rho; rhodot];
        M = M + H'*W*H;
        N = N + H'*W*y(obs);
        epsilon(i,:) = y';
    end
    
    % Plot the residuals
    figure()
    s = stackedplot(Y(:,1),epsilon);
    s.FontSize = 12;
    s.Title = sprintf("%s Residuals, Iteration %i", name, j);
    s.DisplayLabels = ["Range [m]", "Range Rate [m/s]"];
    s.XLabel = "Time [s]";
    s.LineStyle = "-";
    s.LineWidth = 2;
    set(gcf,'Color','w');
    set(gca,'FontSize',12,'FontName','Courier New');
    new = rms(epsilon);

    % Solve the normal equation
    [xhat0, P] = cholesky_decomp(M,N);

    % Compute the standard deviation
    sigma = std_dev(P);

    % Add a row to the data table
    data(j,:) = [Xstar0' sigma' new];

    % Update the values for the next iteration
    Xstar0 = Xstar0 + xhat0;
    x0bar = x0bar - xhat0;
    err = max(abs(new - old));
    j = j + 1;

    % Keep station 1 from changing
    Xstar0(10) = X1(1);
    Xstar0(11) = X1(2);
    Xstar0(12) = X1(3);
end

% Determine the current position of the satellite and plot the updated
% orbit
[~,X] = ode89(@(t,X) ode_func(t,X),Y(:,1),[Xstar0; reshape(eye(18),[],1)],opts);
plot_orbit(name,X);
fprintf("\n%s, Total Iterations: %2i\n",name, j-1)
fprintf("Calculated Final Values\n")
fprintf("r0=[%.1f, %.1f, %.1f]\n",X(end,1),X(end,2),X(end,3))
fprintf("v0=[%.2f, %.2f, %.2f]\n",X(end,4),X(end,5),X(end,6))
fprintf("mu=%.9e\n",X(end,7))
fprintf("J2=%.15e\n",X(end,8))
fprintf("CD=%.4f\n",X(end,9))
fprintf("X1=[%.1f, %.1f, %.1f]\n",X(end,10),X(end,11),X(end,12))
fprintf("X2=[%.1f, %.1f, %.1f]\n",X(end,13),X(end,14),X(end,15))
fprintf("X3=[%.1f, %.1f, %.1f]\n",X(end,16),X(end,17),X(end,18))
end

function sigma = kalman_filter(name,Xstar0,Y,x0bar,P0,R,opts)
% BATCH_PROCESSOR Run the batch processor
% 
% INPUTS:
%     - name   : The name to use for the output data and plots
%     - XStar0 : The initial state vector
%     - Y      : The observation matrix
%     - x0bar  : The a priori x data
%     - P0     : The a priori covariance matrix 
%     - R      : The error matrix
%     - opts   : The options for the integrator
%
% OUTPUTS:
%     - sigma  : The standard deviation vector

% Constants
thetadot = 7.2921158543e-5; %rad/s

% Initialize the values needed for the filter
xhat = x0bar;
P = P0;
epsilon = zeros(length(Y),2);
[~,X] = ode89(@(t,X) ode_func(t, X),Y(:,1),[Xstar0; reshape(eye(18),[],1)],opts);
for i = 1:length(Y)

    % Read the next observation
    t = Y(i,1);
    station = Y(i,2);
    rho_obs = Y(i,3);
    rhodot_obs = Y(i,4);

    % Determine the current station
    if station == 101
        i_xs = 10;
        i_ys = 11;
        i_zs = 12;
    elseif station == 337
        i_xs = 13;
        i_ys = 14;
        i_zs = 15;
    else
        i_xs = 16;
        i_ys = 17;
        i_zs = 18;
    end

    % Read the next state values from the integrator
    x = X(i,1);
    y = X(i,2);
    z = X(i,3);
    xdot = X(i,4);
    ydot = X(i,5);
    zdot = X(i,6);
    mu = X(i,7);
    J2 = X(i,8);
    CD = X(i,9);
    xs = X(i,i_xs);
    ys = X(i,i_ys);
    zs = X(i,i_zs);
    phi = reshape(X(i,19:end),18,18);

    % Compute needed values
    theta = thetadot*t; %rad
    rho = sqrt(x^2 + y^2 + z^2 + xs^2 + ys^2 + zs^2 - 2*(x*xs + y*ys)*cos(theta) + 2*(x*ys - y*xs)*sin(theta) - 2*z*zs); %m
    rhodot = (x*xdot + y*ydot + z*zdot - (xdot*xs + ydot*ys)*cos(theta) + thetadot*(x*xs + y*ys)*sin(theta) + (xdot*ys - ydot*xs)*sin(theta) +thetadot*(x*ys - y*xs)*cos(theta) - zdot*zs)/rho; %m/s

    % Time update
    xbar = phi*xhat;
    P = phi*P*phi';

    % Compute the observation deviation, observation-state matrix, and gain
    % matrix
    y = [rho_obs;rhodot_obs] - [rho;rhodot];
    H = H_tilde([x;y;z],[xdot;ydot;zdot],[xs;ys;zs],rho,rhodot,theta,thetadot,[i_xs;i_ys;i_zs]);
    K = P*H'*pinv(H*P*H' + R);

    % Measurement update
    xhat = xbar + K*(y - H*xbar);
    P = (eye(18) - K*H)*P*(eye(18) - K*H)' + K*R*K';
    Xhat = X(end,1:18)' + xhat;
    epsilon(i,:) = transpose(y);
end

% Calculate the standard deviations
sigma = std_dev(P);

% Plot the residuals
figure()
s = stackedplot(Y(:,1),epsilon);
s.FontSize = 12;
s.Title = sprintf("%s Residuals", name);
s.DisplayLabels = ["Range [m]", "Range Rate [m/s]"];
s.XLabel = "Time [s]";
s.LineWidth = 2;
set(gcf,'Color','w');
set(gca,'FontSize',12,'FontName','Courier New');

% Plot the updated orbit using the computed current location
[~,X] = ode89(@(t,X) ode_func(t, X),[0 5*3600],[Xhat; reshape(eye(18),[],1)],opts);
plot_orbit(name,X)

% Print the current state data
fprintf("\n%s\n",name)
fprintf("Calculated Final Values\n")
fprintf("r0=[%.1f, %.1f, %.1f] m\n",Xhat(1),Xhat(2),Xhat(3))
fprintf("v0=[%.2f, %.2f, %.2f] m/s\n",Xhat(4),Xhat(5),Xhat(6))
fprintf("mu=%.9e m^3/s^2\n",Xhat(7))
fprintf("J2=%.15e\n",Xhat(8))
fprintf("CD=%.4f\n",Xhat(9))
fprintf("X1=[%.1f, %.1f, %.1f] m\n",Xhat(10),Xhat(11),Xhat(12))
fprintf("X2=[%.1f, %.1f, %.1f] m\n",Xhat(13),Xhat(14),Xhat(15))
fprintf("X3=[%.1f, %.1f, %.1f] m\n",Xhat(16),Xhat(17),Xhat(18))
end