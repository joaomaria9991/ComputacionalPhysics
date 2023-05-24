%João Maria Machado, NMEC-89132, PL4
function dxdt = f(t,x,F0,w0,niu,alpha,K)
% x(1)-> y
% x(2)-> v

m=1;
dxdt = zeros(2, 1);  
dxdt(1) = x(2);  
dxdt(2) = (niu*cos(x(2))*x(2)+F0*cos(w0*t)-K*(x(1)+alpha*x(1)^3))/m;