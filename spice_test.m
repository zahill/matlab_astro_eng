clc
clf
addpath("MICE/src/mice");
addpath("MICE/lib")
cspice_furnsh({'MICE/kernels/de440.bsp', 'MICE/kernels/latest_leapseconds.tls','MICE/kernels/pck00011.tpc', 'MICE/kernels/gm_de440.tpc'})
cspice_furnsh({'MICE/kernels/de421.bsp', 'MICE/kernels/lrorg_2023258_2023349_v01.bsp'})

n = 10000;
et = cspice_str2et('2023 Sep 15 00:00:00.000');
tspan = et + linspace(0,10*24*60*60,n);

TARGET = cspice_bods2c('LRO')
OBS = cspice_bods2c('MOON')

[pos, t] = cspice_spkpos('LRO', tspan, 'J2000', 'NONE', 'MOON');
[state, ~] = cspice_spkez(TARGET,et,'J2000','NONE',OBS);

r = zeros(n,1);
for i = 1:n
    r(i) = sqrt(pos(1,i)^2 + pos(2,i)^2 + pos(3,i)^2);
end

moon = dictionary('mass',0.07346e24,'mu',4902,'radius',1738,'J2',202.7e-6);
mean(r) - moon('radius')

cspice_gdpool('BODY301_GM',0,1)

figure(3);
hold on;
plot3(pos(1,:),pos(2,:),pos(3,:),'b','LineWidth',2);
[x,y,z] = sphere(100);
x = x*moon('radius');
y = y*moon('radius');
z = z*moon('radius');
surf(x,y,z,'EdgeColor','k','FaceColor','#808080');
axis equal;
grid on;
set(gcf,'Color','w');
set(gca,'Color','w');
set(gca,'FontSize',12,'FontName','Courier New');
title('LRO')
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
view(-45,30);
legend('SPICE')
hold off;

cspice_kclear;