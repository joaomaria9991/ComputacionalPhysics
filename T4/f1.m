function dxdt = f1(t,x)

% x(1)-> y
% x(2)-> v


epsilon = 1;
F0=1;

dxdt = zeros(2, 1);  

dxdt(1) = x(2);  
dxdt(2) = F0*cos(1.7*t)-epsilon * (x(1)^2 - 1) * x(2) - x(1);   
