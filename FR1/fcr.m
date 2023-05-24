function F = fcr(xv,xold,vold,const)
% const(1), const(2) e const(3) est?o definidas no programa principal.
% xold ? x(k) e vold ? vx(k).
% xv(1) ? x(k+1) e xv(2) ? vx(k+1).

x=xv(1);
v=xv(2);

F(1)= x - xold - const(1) * (v + vold);

F(2)= v - vold + const(2) * (x + xold + const(3) * (x^3 + xold^3));


end