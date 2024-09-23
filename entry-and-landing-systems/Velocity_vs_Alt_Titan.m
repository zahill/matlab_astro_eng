clc
clear
close all


% Define Entry Conditions
entryVel = 4 * 1000; % m/s
entryAlt = 1000 * 1000; % m
FPA = -50; % deg
Cd_capsule = 1;
d_capsule = 4; % m
entryMass = 1100; % kg


% Define Parachute Paramiters
number_of_chutes = 2;
cluster_performance = 0.05; % reduction per chute
d_chute = 3; % m
h_open_chute = 140*1000; % m
theta_trim = 5; % deg
Cd_chute = 1.65;


% call functions to calculate vel vs alt durring each phase
[V_ball,S_ball,h_ball,t_ball] = ballistic_entry_phase(entryVel,entryAlt,h_open_chute,FPA,Cd_capsule,d_capsule,entryMass);
[V_chute,S_chute,h_chute,t_chute] = parachute_phase(h_open_chute,V_ball(1),d_capsule,Cd_capsule,entryMass,number_of_chutes,cluster_performance,d_chute,theta_trim,Cd_chute,S_ball(end));
V_chute = [V_ball(1),V_chute(2:end)];
S_chute = [S_chute(1),S_chute(2:end)];
h_chute = [h_ball(1),h_chute(2:end)];
S = calculate_downrange(flip(h_ball),h_chute,FPA,theta_trim);

function [V,S,h,t] = parachute_phase(h_open_chute,v_initial,d_capsule,Cd_capsule,entryMass,number_of_chutes,cluster_performance,d_chute,theta_trim,Cd_chute,S_start)
% Define Titan atmospheric conditions
    % mollar mass of atmosphere species
    M = 28.02/1000; %kg/mol
    % mass in kg
    m_molecule = M./(6.02*(10^23));
    % boltzmans contatnt 
    k = 1.38*(10^(-23)); %J/K
    %specific gas constant
    R = 8.31./M; %J/kgK
    % mean surface temp
    T0 = -179 + 273.15; %K
    % surface pressure
    P0 = 1.6*101325; %Pa
    % atmospheric density at surface
    rho0 = P0/(R*T0); %kg/m^3
    % surface gravitational constant
    g = 1.3; %m/s^2


% calculate the scale height
A = (g.*m_molecule)./(k.*T0); %scale height [1/m]

% calculate area of parachute and entry vehicle
S_p = pi.* ((d_chute ./ 2).^2); % m^2
S_EV = pi.* ((d_capsule ./ 2).^2); % m^2 

% find height and velocity as a function of time
delta_t = 0.1; % s
% calculate drag
if number_of_chutes > 1
    dragComp = (Cd_chute*S_p*number_of_chutes*(1-cluster_performance))+(Cd_capsule*S_EV);
else
    dragComp = (Cd_chute*S_p)+(Cd_capsule*S_EV);
end
%calculate gravity component
gravityComp = g*cosd(theta_trim);
% set initial conditions
h_old = h_open_chute;
v_old = v_initial;
s_old = S_start;
j = 1;
h(j) = h_old;
V(j) = v_old;
S(j) = s_old;
t(j) = delta_t;
while h_old >= 0
j = j + 1;
h_new = h_old - (v_old* cosd(theta_trim) * delta_t);
s_new = s_old - (v_old* sind(theta_trim) * delta_t);
density = rho0*exp(-A*h_old);
v_new = v_old + ((gravityComp - (1/(2*entryMass))*density*(v_old^2)*dragComp)*delta_t);
dragForce(j) = 0.5*density*(v_old^2)*(Cd_chute*S_p);
h(j) = h_new;
S(j) = s_new;
V(j) = v_new;
h_old = h_new;
v_old = v_new;
t(j) = t(j-1) + delta_t;
end
end



function [V,X,h,t,N_max,Hn_max] = ballistic_entry_phase(entryVel,entryAlt,h_open_chute,FPA,Cd,d_capsule,entryMass)

% shorten variable names
Ve = entryVel;
d = d_capsule;
me = entryMass;

% Area of entry vehicle
S = pi*((d./2).^2); % m^2

% Define Titan atmospheric conditions
    % mollar mass of atmosphere species
    M = 28.02/1000; %kg/mol
    % mass in kg
    m_molecule = M./(6.02*(10^23));
    % boltzmans contatnt 
    k = 1.38*(10^(-23)); %J/K
    %specific gas constant
    R = 8.31./M; %J/kgK
    % mean surface temp
    T0 = -179 + 273.15; %K
    % surface pressure
    P0 = 1.6*101325; %Pa
    % atmospheric density at surface
    rho0 = P0/(R*T0); %kg/m^3
    % surface gravitational constant
    g0 = 1.3; %m/s^2

% create a vector for hight above planet surface
h = h_open_chute:1:entryAlt; % m

% calculate the scale height
A = (g0.*m_molecule)./(k.*T0); %scale height [1/m]

%calculate ballistic coefficient
beta = me./(Cd.*S);

% create vector of velocity vs alt
C = rho0 ./ (2.*beta.*A.*sind(FPA));
V = Ve.*exp(C.*exp(-A.*h));

t = 0;
V2 = flip(V);
X(1) = 0;
for i = 1:1:length(h)
    if i<length(h)
        t_pass = (abs(h(i) - h(i+1))) / V(i);
    end
    if i>1 
        X(i) = X(i-1) - (V2(i-1)* cosd(FPA) * t_pass);
    end
    t = t + t_pass;
end

end

function S = calculate_downrange(h_ball,h_chute,FPA,theta_trim)
    h = [h_ball, h_chute];
    FPA = [FPA*ones(1,length(h_ball)),(theta_trim - 90)*ones(1,length(h_chute))];
    S(1) = 0;
    for i = 2:length(h)
        S(i) = S(i-1) + cosd(FPA(i-1))*abs(h(i) - h(i-1));
    end

    fig = figure();
    pl = plot(h/1000,S/1000);
    pl.LineWidth = 2;
    ax = gca;
    ax.XLim = [0 h_ball(1)/1000];
    ax.FontName = "Times New Roman";
    ax.FontSize = 16;
    box(ax,"on");
    title(ax,"Downrange Dist vs. Altitude");
    subtitle(ax,sprintf("Total Downrange Dist: %.2f km",S(end)/1000));
    ylabel(ax,"Downrange Distance [km]");
    xlabel(ax,"Altitude [km]");
    xl = xline(h_chute(1)/1000,'--',{'Parachute','Deployed'});
    xl.LineWidth = 2;
    xl.LabelVerticalAlignment = 'bottom';
    xl.FontSize = 16;
    xl.LabelOrientation = 'horizontal';
end


function formatPlot
set(gca,'FontSize',16,"Fontname","Times New Roman",'XminorTick','off')
box on
end 














