function dxdt = f(t,y,K,m,niu,alpha)


dxdt = zeros(2, 1);   % this vector of two elements represents:  [dy/dt; dv/dt]

x = y(1);
v = y(2);


dxdt(1) = v;   % derivative of "y"
dxdt(2) =(-K*(x+alpha*x.^2)+niu*sin(v)*v)/m; %derivative of "v"
