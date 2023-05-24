function dxdt = f(t,y,K,m,alpha)
m=1.5;
K=2;


dxdt = zeros(2, 1);   % this vector of two elements represents:  [dy/dt; dv/dt]

x = y(1);
v = y(2);


dxdt(1) = v;   % derivative of "y"
dxdt(2) = -K*x*(1+(3/2)*alpha*x)*1/m; %derivative of "v"
