close all
clear all
clc

%Constantes

m=1;
K=1;
alpha=0.1;

ti=0;
tf=10;


x0 = 1;
v0 = 1;

hvals = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5];


for i=1:length(hvals)
    h=hvals(i);
    
    t=ti:h:tf;
    
    
    options=optimset('Display','off','Tolx', 1e-10,'TolFun',1e-10);
    
    const = [h/2, K*h/(2*m), 2*alpha];
    N=length(t);
    
    x=zeros(1,N);
    vx=zeros(1,N);
    
    x(1)=1;
    vx(1)=1;
    
    for k=1:N-1
        func= @(xv) fcr(xv,x(k),vx(k),const);
        
        
        
        xv0=[x(k), vx(k)];
        aux=fsolve(func,xv0,options);
        
        x(k+1)=aux(1);
        vx(k+1)=aux(2);
        
        
    end
    hsquare(i) = h^2;
    xf1(i) = x(end);
    
    
end


% subplot(2,1,2)
plot(hsquare, xf1, 'rs-')
title('Erro-Crank-Nicolson')

















