% function to solve:

function dydt = func(t, y)
% "dydt" means the derivative "dy/dt"

% physical parameters:
w0 = 1;
q = 1/2;
wD = 2/3;
FD = 1.2;

dydt = zeros(2, 1);

% the derivative of x (theta)
dydt(1) = y(2);   % the derivative of "x" is just "v"

% the derivative of v (theta's derivative)
dydt(2) = -w0 * sin( y(1) ) - q * y(2) + FD * sin(wD * t);

