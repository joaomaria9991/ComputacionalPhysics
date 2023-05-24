close all
clear all
clc

%Constantes 

k=1;
m=1;
ti=0;
tf=50;
h=0.01;
w=1; %vem da física, o resto vem da computação.

%Vetores e condições iniciais

t=ti:h:tf;
x=zeros(1,length(t));

x(1)=1;
v(1)=0;


%Integrador Euler 

for i=1:length(t)-1
    a=-k*x(i);
    v(i+1)=a*h+v(i);
    x(i+1)=v(i)*h+x(i);
    EM1(i+1)=1/2*v(i)^2+1/2*k*x(i)^2;

    
end
figure(1)
subplot(2,2,1)
plot(t,x);
title('Método de Euler')

figure(2)
subplot(2,2,1);
plot(t,EM1);
title('Energia Mecânica segundo Método de Euler')




%Integrador Euler-Cromer


for i=1:length(t)-1
       a=-k*x(i);
       v(i+1)=v(i)+a*h;
       x(i+1)=v(i+1)*h+x(i);
       EM2(i+1)=1/2*v(i)^2+1/2*k*x(i)^2;

      
end

figure(1)
subplot(2,2,2);
plot(t,x);
title('Método de Euler-Cromer')

figure(2)
subplot(2,2,2);
plot(t,EM2);
title('Energia Mecânica segundo Método de Euler-Cromer')


%Integrador Euler Implícito

for i=1:length(t)-1
    a=-k*x(i);
    x(i+1)=(x(i)+v(i)*h)/(1+h^2);
    v(i+1)=v(i)-x(i+1)*h;
    EM3(i+1)=1/2*v(i)^2+1/2*k*x(i)^2;


end

figure(1)
subplot(2,2,3);
plot(t,x);
title('Método de Euler-Implicito')

figure(2)
subplot(2,2,3);
plot(t,EM3);
title('Energia Mecânica segundo Método de Euler-Implícito')



%Integrador Crank-Nicolson

A=[1 -h/2; w^2*h/2 1];


for i=1:length(t)-1
B=[x(i)+(h/2)*v(i); v(i)-(h/2)*x(i)];
Z=linsolve(A,B);

x(i+1)=Z(1);
v(i+1)=Z(2);
EM4(i+1)=(1/2)*v(i+1)^2+(1/2)*k*x(i+1)^2;
end


figure(1)
subplot(2,2,4);
plot(t,x);
title('Método de Crank-Nicolson')

figure(2)
subplot(2,2,4);
plot(t,EM4);
title('Energia Mecânica segundo Método de Crank-Nicolson')
%axis([0 tf min(EM4) max(EM4)])

%%%%%%%%%%%%%%%%%%%%%%%%%
num=0;
for i=2:length(t)-1
    if (x(i-1)<=x(i) && x(i)>=x(i+1))
        num=num+1;
        aux=lagr(t(i-1:i+1),x(i-1:i+1));
        tvals(num)=aux(1);
        xvals(num)=aux(2);
    end

end

%Amplitude 
Amp_num=mean(xvals);
Amp_an=x(1);

%Period

T_num=(tvals(end)-tvals(1))/(num-1)






