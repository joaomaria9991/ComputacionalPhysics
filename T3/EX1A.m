close all
clear all
clc

%Constantes

k=16;
m=1;
w=1;

%Vetores e condições iniciais
ti=0;
tf=15;  %500; 12000;
h=0.01;

t=ti:h:tf;

x=zeros(length(t),1);
v=zeros(length(t),1);

x(1)=1;
v(1)=0;

%Funções Implicitas
fv= @(x) -k*x/m;
fx= @(v) v;


%Integrador Runge-Kutta 2ª Ordem
% 
for i=1:length(t)-1
    r1v=fv(x(i));
    r1x=fx(v(i));

    r2v=fv(x(i)+r1x*(h/2));
    r2x=fx(v(i)+r1v*(h/2));


    v(i+1)=v(i)+r2v*h;
    x(i+1)=x(i)+r2x*h;


    Em_RK(i+1)=0.5*m*v(i).^2+0.5*k*x(i)^2;
end

for i=1:length(t)-1
    a=-k/m*x(i);
    r1v=a;
    r1x=v(i);
    
    r2v= -k/m *(x(i)+h/2*r1x);
    r2x=v(i)+h/2*r1v;
    
    
    v(i+1)=v(i)+r2v*h;
    x(i+1)=x(i)+r2x*h;
    
    Em_RK(i+1)=0.5*m*v(i).^2+0.5*k*x(i)^2;
    
end


subplot(2,1,1)
plot(t,Em_RK,'-b')
title('Energia Mecânica-Runge-Kutta')


%Integrador Euler
for i=1:length(t)-1
    a=-k/m*x(i);
    v(i+1)=v(i)+a*h;
    x(i+1)=x(i)+v(i)*h;
    Em(i+1)=0.5*m*abs(v(i))^2+0.5*k*x(i)^2;
    
end

subplot(2,1,2)
plot(t,Em,'-r')
title('Energia Mecânica-Euler')






