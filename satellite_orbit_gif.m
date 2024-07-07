global mu R;

mu = 398600;
R = 6400;

% Inputs
x = 10000;
y = 1000;
z = 0;
xdot = 1;
ydot = sqrt(mu/x);
zdot = 1;

tspan = linspace(0,12*3600,1.0e3);

options = odeset('RelTol',1.0e-6,'InitialStep',1.0e-6,'AbsTol',1.0e-6);

[t,X] = ode45(@ode_func,tspan,[x,y,z,xdot,ydot,zdot],options);

filename = 'orbit.gif';
fig2 = figure(101);
hold on;
[x,y,z] = sphere(100);
x = x*R;
y = y*R;
z = z*R;
surf(x,y,z,'EdgeColor','k','FaceColor','#808080');
axis equal;
grid on;
set(gcf,'Color','w');
set(gca,'Color','w');
set(gca,'FontSize',22,'FontName','Courier New');
title(sprintf('Orbit\nTime: %0.2f sec', t(1)),'Interpreter','Latex');
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
plot3(X(:,1),X(:,2),X(:,3),'Color','none');
p = plot3(X(1,1),X(1,2),X(1,3),'b', 'LineWidth', 2);
m = scatter3(nan,nan,nan,'filled','b');
view(-45,30);

for k = 1:length(t)
    p.XData = X(1:k,1);
    p.YData = X(1:k,2);
    p.ZData = X(1:k,3);

    m.XData = X(k,1);
    m.YData = X(k,2);
    m.ZData = X(k,3);

    title(sprintf('Orbit\nTime: %0.2f sec', t(k)),'Interpreter','Latex');

    %pause(1.0e-3);

    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if k == 1
        imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',0.1);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
end
hold off;

function Xdot = ode_func(~,X)
global mu

x = X(1);
y = X(2);
z = X(3);
xdot = X(4);
ydot = X(5);
zdot = X(6);

r=sqrt(x^2+y^2+z^2);

Xdot(1,1) = xdot;
Xdot(2,1) = ydot;
Xdot(3,1) = zdot;
Xdot(4,1) = -mu * x / r^3;
Xdot(5,1) = -mu * y / r^3;
Xdot(6,1) = -mu * z / r^3;
end