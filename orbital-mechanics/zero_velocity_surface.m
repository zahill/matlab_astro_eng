close all
clear all

C = -2;

dinc = 0.005;
xvec = -1.5:dinc:1.5;
yvec = -1.5:dinc:1.5;
[xsm,ysm] = meshgrid(xvec,yvec);
z = 0;

nx = length(xvec);
ny = length(yvec);

velocity_surf = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        x = xsm(i,j);
        y = ysm(i,j);

        r = sqrt(x^2 + y^2);
        V = 1/r + 3/2*x^2;
        if (V + C) > 0
            velocity_surf(i,j) = 1;
        else
            velocity_surf(i,j) = 0;
        end
    end
end

fig = figure();
hold on
contourf(xsm,ysm,velocity_surf,100)
colormap('gray')
hold off

M1 = 5.97219e24;
M2 = 7.34767e22;
mu = M2/(M1 + M2);
velocity_surf = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        x = xsm(i,j);
        y = xsm(i,j);

        r1 = sqrt((x + mu)^2 + y^2);
        r2 = sqrt((x - 1 + mu)^2 + y^2);
        V = (x^2 + y^2)/2 + (1 - mu)/r1 + mu/r2;
        if (V + C) > 0
            velocity_surf(i,j) = 1;
        else
            velocity_surf(i,j) = 0;
        end
    end
end

fig = figure();
hold on
contourf(xsm,ysm,velocity_surf,100)
colormap('gray')
hold off