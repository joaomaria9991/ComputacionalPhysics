close all 
clear all
clc

%Constantes

g=9.8;
G=6.67*10^-11;
ms=1.989*10^33;
mm=3.285*10^23;
Gms=4*pi^2;

UA=149597870700;

%Vetores e condições iniciais 

h=0.0001; %anos

ti=0;
tf=1;
t=ti:h:tf;

x=zeros(length(t),1);
y=zeros(length(t),1);
vx=zeros(length(t),1);
vy=zeros(length(t),1);
ang=zeros(length(t),1);
ang2=zeros(length(t),1);


vx(1)=0;
x(1)=0.47;%*UA;
vy(1)=8.2; %UA por ano
y(1)=0;
ang(1)=0;
ang2(1)=0;
nr=0; %número de rotações a volta do sol


for i=1:length(t)-1
       r=sqrt(x(i)^2+y(i)^2);
    
       ax=-4*pi^2/(r^2)*x(i)/r;
       ay=-4*pi^2/(r^2)*y(i)/r;

       
    
       vx(i+1)=vx(i)+ax*h;
       vy(i+1)=vy(i)+ay*h;

       
       x(i+1)=vx(i+1)*h+x(i);      
       y(i+1)=vy(i+1)*h+y(i);
       
       ang(i+1)=atan2(y(i+1),x(i+1));       
       if(ang(i+1)<ang(i))
        nr=nr+1;
       end
     ang2(i+1)=ang(i+1)+nr*2*pi;
             
end


axis([-0.5 0.5 -0.5 0.5])
set(gca,'PlotBoxAspectRatio',[1 1 1])
plot(x,y)
hold on
plot(0,0,'r*');

figure(2)
plot(t,ang,'k');
hold on
plot(t,ang2,'b-');

%calcular o periodo

T1=interp1(ang2,t,2*pi);
T2=interp1(ang2,t,4*pi);
T=T2-T1





