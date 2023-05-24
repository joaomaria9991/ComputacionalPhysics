close all
clear all
clc

ti=0;
tf=5000;
h=0.82;

t=ti:h:tf;

%vetores e condi??es iniciais 

x=zeros(1,length(t));
y=zeros(1,length(t));

y(1)=0.01;
x(1)=0.01;


%Integrador Euler 


for i=1:length(t)-1
    dxdt=y(i);
    dydt=(1-x(i).^2-y(i).^2)*y(i)-x(i);
    
    y(i+1)=y(i)+dydt*h;
    x(i+1)=x(i)+dxdt*h;
    
end

subplot(1,2,1)
plot(x,y,'.')
title('M?todo de Euler')


%fun??es an?nimas RK_4

 fx= @(y) y;
 fy= @(x,y) (1-x^2-y^2)*y-x;


%Integrador Runge-Kutta 4? Ordem

for i=1:length(t)-1
  r1v=fv(y(i),v(i),t(i));
  r1y=fy()
  
    
end


subplot(1,2,2)
plot(x,y,'.')
title('M?todo de R.K 4? Ordem')



