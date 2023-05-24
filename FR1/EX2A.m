close all
clear all
clc

%Constantes

m = 1;
K= 1;
alpha = 0.1;




ti=0;

tf=10;

x0 = 1;
v0 = 1;
t0 = 0;
tf = 10;





hvals = [0.001, 0.002, 0.005, 0.01, 0.02];


for i=1:length(hvals)
    
    
    h=hvals(i);
    t=ti:h:tf;
    
    N=length(t);
    
    
    x = zeros(N, 1);
    v = zeros(N, 1);
    
    
    x(1) = x0;
    v(1) = v0;
    
    
    for k=1:N-1
       
        
        v(k+1) = v(k) - h * K / m * (x(k) + 2 * alpha * x(k)^3);
        x(k+1) = x(k) + h * v(k+1);
    end
    xf(i) = x(end);
end

%subplot(2,1,1)
plot(hvals, xf, 'rs-')
title('Erro-Euler-Cromer')












