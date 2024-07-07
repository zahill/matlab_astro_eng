mu = 4.9028001184575496e3;
R = 1738;
J2 = 202.7e-6;
a = 8500;
e = linspace(0,0.75,100);
i = 65;

x = -3/2.*mu.*R^2.*J2./(a.^3.*(1 - e.^2).^(3/2)).*sin(i).*cos(i);
plot(e,x)
