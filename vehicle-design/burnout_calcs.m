h = 640; %km
R = 6400; %km
mu = 3.986e5; %km^3/s^2
v_bo = 10.4; %km/s

r_bo = h + R; %km
Q_bo = v_bo^2 * r_bo / mu; %-

phi_bo = linspace(0,90,1.0e3);
psi = 2 * acosd((1-Q_bo .* cosd(phi_bo).^2) ./ sqrt(1+Q_bo.*(Q_bo-2).*cosd(phi_bo).^2));

fig = figure(101);
plot(phi_bo,psi);
title(sprintf('Range vs. Flight Angle\nQ = %0.2f',Q_bo),'Interpreter','latex');
xlabel('\phi_{bo} [deg]');
ylabel('\Psi [deg]');

a = r_bo/(2 - Q_bo);
e = sqrt(1 + Q_bo*(Q_bo - 2).*cosd(phi_bo).^2);
r_a = a.*(1+e);
h_a = r_a - R;

fig = figure(102);
plot(phi_bo,h_a);
title('Altitude @ Apogee vs. Flight Angle','Interpreter','latex');
xlabel('\phi_{bo} [deg]');
ylabel('Altitude [km]');

altitude(90);
altitude(88);

function altitude(phi_bo)
h = 640; %km
R = 6400; %km
mu = 3.986e5; %km^3/s^2
v_bo = 10.4; %km/s

r_bo = h + R; %km
Q_bo = v_bo^2 * r_bo / mu; %-

a = r_bo/(2 - Q_bo);
e = sqrt(1 + Q_bo*(Q_bo - 2).*cosd(phi_bo).^2);
b = a * sqrt(1 - e^2);
r_a = a.*(1+e);
r_p = a.*(1-e);
h_a = r_a - R;
h_p = r_p - R;
fprintf('---------------------------------\n');
fprintf('Flight angle @ %i deg\n',phi_bo);
fprintf('Altitude @ Apogee: %.2f km\n', h_a);
fprintf('Altitude @ Perigee: %.2f km\n', h_p);
fprintf('---------------------------------\n');
end