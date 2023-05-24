function dxdt = f(t,x)

% x(1) means "y"
% x(2) means "v"

% physical parameters:
epsilon = 1;


dxdt = zeros(2, 1);   % this vector of two elements represents:  [dy/dt; dv/dt]

dxdt(1) = x(2);   % derivative of "y"
dxdt(2) = -epsilon * (x(1)^2 - 1) * x(2) - x(1);   % derivative of "v"
