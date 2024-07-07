%Earth Properties
mu = 3.986e14; %m^3/s^2
R = 6378e3; %m
g0 = 9.80665; %m/s^2

%Orbit Properties
h1 = 250e3; %m
h2 = 926e3; %m
r1 = R + h1; %m
r2 = R + h2; %m

%Launch Properties
L = 28.45 %deg
h_LS = 3 %m
v_LS = 465.1*cosd(L) %m/s
beta = 90 %deg

%Flight Values
drag_loss = 200 %m/s
g_loss = sqrt(2*mu*(1/R - 1/r1))*0.8 %m/s
v_cs1 = sqrt(mu/r1) %m/s
v1 = sqrt(2*mu*(1/r1 - 1/(r1 + r2))) %m/s
delta_v1 = v1 - v_cs1 %m/s
v_cs2 = sqrt(mu/r2) %m/s
v2 = sqrt(2*mu*(1/r2 - 1/(r1 + r2))) %m/s
delta_v2 = v2 - v_cs2 %m/s

v_design = v1 - v_LS + g_loss + drag_loss %m/s