close all;
clc;

%Lunar Properties
mu = 4902.8e9; %m^3/s^2
R = 1740e3; %m
T = 655.7*3600; %s

%Orbit Properties
h1 = 30e3; %m
h2 = 250e3; %m
r1 = R + h1; %m
r2 = R + h2; %m
v_LS = 2*pi*R/T %m/s
v1 = sqrt(2*mu*(1/r1 - 1/(r1+r2))) %m/s
v2 = sqrt(2*mu*(1/r2 - 1/(r1+r2))) %m/s
v_CM = sqrt(mu/r2) %m/s
delta_v2 = v2 - v_CM %m/s
a = (r1 + r2)/2; %m
t_coast = pi*sqrt(a^3/mu) %s
t_LW = pi*sqrt(r2^3/mu) - t_coast %s
theta_LW = rad2deg(v_CM*t_LW/r2) %deg

%LM Properties
m_tot = 4780; %kg
m_p = 2375; %kg
T = 15.57e3; %N
Isp = 311; %s
g0 = 9.8067; %m/s^2

%Variables
gamma_kick = 7.565; %deg
delta_t = 1; %s

%Initial Values
t = 0; %s
a = T/m_tot - mu/R^2; %m/s^2
gammadot = 0;
v = 0; %m/s
gamma = 90; %deg
v_x = 0;
x = 0; %m
v_y = 0;
h = 0; %m
g = mu/R^2; %m/s^2
m = m_tot; %kg
X = [t,a,gammadot,v,gamma,v_x,x,v_y,h,g,m];
kick = false;

%Iterate until the velocity meets the needed condition or the LV falls back
%to the surface of the moon
while v < (v1 - v_LS) && h >= 0
    %Implement the kick once when the altitude is greater than 100m
    if h >= 100 && kick == false
        gamma = gamma - gamma_kick;
        kick = true;
    end
    
    %Update the values
    t = t + delta_t; %s
    mdot = T/(g0*Isp); %kg/s
    a = T/m - g*sind(gamma); %m/s^2
    if v ~= 0
        gammadot = 180/pi*(v^2/(R + h) - g)/v*cosd(gamma); %deg/s
    else
        gammadot = 0; %deg/s
    end
    v = v + a*delta_t; %m/s
    gamma = gamma + gammadot*delta_t; %deg
    v_x = v*cosd(gamma); %m/s
    x = x + v_x; %m
    v_y = v*sind(gamma); %m/s
    h = h + v_y; %m
    g = mu/(R + h)^2; %m/s^2
    m = m - mdot*delta_t; %kg
    X = [X;[t,a,gammadot,v,gamma,v_x,x,v_y,h,g,m]];
end

%Values at shutdown
tb = t %s
vb = v/1000 %km/s
hb = h/1000 %km
gammab = gamma %deg

%Allow for the flight simulation to continue for 50s after burnout
for i = 1:50*delta_t
    t = t + delta_t; %s
    a = -g*sind(gamma); %m/s^2
    gammadot = 180/pi*(v^2/(R + h) - g)/v*cosd(gamma); %deg/s
    v = v + a*delta_t; %m/s
    gamma = gamma + gammadot*delta_t; %deg
    v_x = v*cosd(gamma); %m/s
    x = x + v_x; %m
    v_y = v*sind(gamma); %m/s
    h = h + v_y; %m
    g = mu/(R + h)^2; %m/s^2
    X = [X;[t,a,gammadot,v,gamma,v_x,x,v_y,h,g,m]];
end

%Specify plot colors
blue = [57 106 177]./255;
black = [204 37 41]./255;
green = [62 150 81]./255;
brown = [146 36 40]./255;
purple = [107 76 154]./255;

%Create Figure 101
figure(101);
hold on;
grid on;
set(gcf,'Color','w');
set(gca,'Color','w');
set(gca,'FontSize',22,'FontName','Courier New');
title('LM Properties');
xlabel('Time [s]');
legend('Location','southoutside');
yyaxis left
ylim([0 100]);
plot(X(:,1),X(:,2),'DisplayName','Acceleration [m/s^2]','LineWidth',2,'Color',blue);
plot(X(:,1),X(:,5),'DisplayName','Flight Path Angle [deg]','LineWidth',2,'Color',black);
plot(X(:,1),X(:,9)/1000,'DisplayName','Altitude [km]','LineWidth',2,'Color',green);
plot(X(:,1),X(:,7)/1000,'DisplayName','Downrange Dist. [km]','LineWidth',2,'Color',brown);
yyaxis right
plot(X(:,1),X(:,4),'DisplayName','Velocity [m/s]','LineWidth',2,'Color',purple);
hold off;

%Create Figure 102
figure(102);
hold on;
grid on;
set(gcf,'Color','w');
set(gca,'Color','w');
set(gca,'FontSize',22,'FontName','Courier New');
title('Altitude vs. Downrange Distance');
xlabel('Downrange Distance [km]');
ylabel('Altitude [km]');
plot(X(:,7)/1000,X(:,9)/1000,'LineWidth',2,'Color',blue);
hold off;

%Calculate gravitational losses
g_loss = sum(X(:,10).*sind(X(:,5))*delta_t) %m/s
g_loss = sqrt(2*mu*(1/R - 1/r1))*0.8 %m/s

%Calculate the max acceleration
a_max = max(X(:,2))/g0 %g