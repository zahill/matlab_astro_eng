clc;
clf;

%Earth Properties
mu = 3.986e14; %m^3/s^2
R = 6378e3; %m
g0 = 9.80665; %m/s^2
L = 28.45; %deg
C_D = 0.2; %-
d = 20; %m
S = pi*(d/2)^2; %m^2
h0 = 7194; %m
rho0 = 1.225; %kg/m^3
rho = @(h) rho0*exp(-h/h0); %kg/m^3
g = @(h) g0/(1 + h/R)^2; %m/s

%Orbit Properties
v_s = 2006; %m/s
v_f = 7944; %m/s

%SSTO Properties
m_tot = 61398.9e3; %kg
m_f = 2340.6e3; %kg
T = 782752.77e3; %N
Isp = 301; %s

%Variables
gamma_kick = 0.1593; %deg
delta_t = 1; %s

%Initial Values
t = 0; %s
h = 0; %m
v = 0; %m/s
x = 0; %m
gamma = 90; %deg
a = T/m_tot - g0; %m/s^2
m = m_tot; %kg
m_p = 0; %kg
W = m*g0; %N
D = 1/2*C_D*S*rho0*v^2; %N
X = [t,h,rho(h),v,x,gamma,a,m,m_p,W,D,g0];
kick = false;

%Iterate until the velocity meets the needed condition, all of the
%propellant is used up, or the LV falls back to the Earth's surface
while v < v_s && h >= 0
    %Implement the kick once the altitude is greater than 1km
    if h >= 1000 && kick == false
        gamma = gamma - gamma_kick;
        kick = true;
    end

    %Update the values
    if v ~= 0
        gammadot = 180/pi*(v^2/(R + h) - g(h))/v*cosd(gamma); %deg/s
    else
        gammadot = 0; %deg/s
    end
    t = t + delta_t; %s
    h = h + v*sind(gamma); %m
    v = v + a*delta_t; %m/s
    x = x + v*cosd(gamma); %m
    gamma = gamma + gammadot*delta_t; %deg
    a = (T - D)/m - g(h)*sind(gamma); %m/s^2
    mdot = T/(Isp*g0); %kg/s
    m = m - mdot*delta_t; %kg
    m_p = m_p + mdot*delta_t; %kg
    W = m*g(h); %N
    D = 1/2*C_D*S*rho(h)*v^2; %N
    X = [X;[t,h,rho(h),v,x,gamma,a,m,m_p,W,D,g(h)]];
end

m_tot = 13826.4e3; %kg
m_f = 783.7e3; %kg
m_p = 0; %kg
T = 94913e3; %N
Isp = 363; %s

while v < v_f && h >= 0
    %Update the values
    gammadot = 180/pi*(v^2/(R + h) - g(h))/v*cosd(gamma); %deg/s
    t = t + delta_t; %s
    h = h + v*sind(gamma); %m
    v = v + a*delta_t; %m/s
    x = x + v*cosd(gamma); %m
    gamma = gamma + gammadot*delta_t; %deg
    a = (T - D)/m - g(h)*sind(gamma); %m/s^2
    mdot = T/(Isp*g0); %kg/s
    m = m - mdot*delta_t; %kg
    m_p = m_p + mdot*delta_t; %kg
    W = m*g(h); %N
    D = 0; %N
    X = [X;[t,h,rho(h),v,x,gamma,a,m,m_p,W,D,g(h)]];
end

%Values at shutdown
tb = t %s
vb = v/1000 %km/s
hb = h/1000 %km
gammab = gamma %deg

%Specify plot colors
blue = [57 106 177]./255;
black = [204 37 41]./255;
green = [62 150 81]./255;
brown = [146 36 40]./255;
purple = [107 76 154]./255;

%Create Figure 201
figure(201);
hold on;
grid on;
set(gcf,'Color','w');
set(gca,'Color','w');
set(gca,'FontSize',22,'FontName','Courier New');
title('SSTO Performance');
xlabel('Time [s]');
legend('Location','southoutside');
yyaxis left
%ylim([0 120]);
plot(X(:,1),X(:,7),'DisplayName','Acceleration [m/s^2]','LineWidth',2,'Color',blue);
plot(X(:,1),X(:,6),'DisplayName','Flight Path Angle [deg]','LineWidth',2,'Color',black);
plot(X(:,1),X(:,2)/1000,'DisplayName','Altitude [km]','LineWidth',2,'Color',green);
yyaxis right
plot(X(:,1),X(:,4)/1000,'DisplayName','Velocity [km/s]','LineWidth',2,'Color',brown);
plot(X(:,1),X(:,3).*X(:,4).^2/1000,'DisplayName','Dynamic Pressure [kPa]','LineWidth',2,'Color',purple);
hold off;

t = X(67,1) %s
h = X(67,2)/1000 %km
v = X(67,4) %m/s
rho = X(67,3) %kg/m^3
M = v/303.32 %-
m = X(67,8)/1000 %T
q_max = max(X(:,3).*X(:,4).^2) %Pa
m_p = X(67,9)/1000 %T
a = 782752.77/m/g0 %g