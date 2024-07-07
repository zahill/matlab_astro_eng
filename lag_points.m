syms x y mu;

V = 1/2*(x^2 + y^2) + (1 - mu)/sqrt((x + mu)^2 + y^2) + mu/sqrt((x - 1 + mu)^2 + y^2);
dVdx = diff(V,x);
dVdx = subs(dVdx,y,0);
fzero(dVdx,x0)