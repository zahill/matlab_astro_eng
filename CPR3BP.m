function dX = CPR3BP(~,mu,X)
x = X(1);
y = X(2);
xdot = X(3);
ydot = X(4);

r1 = sqrt((x + mu)^2 + y^2);
r2 = sqrt((x - 1 + mu)^2 + y^2);

dX = zeros(length(X),1);
dX(1) = xdot;
dX(2) = ydot;
dX(3) = 2*ydot + x - (1 - mu)*(x + mu)/r1^3 - mu*(x + mu - 1)/r2^3;
dX(4) = -2*xdot + y*(1 - (1 - mu)/r1^3 - mu/r2^3);
end