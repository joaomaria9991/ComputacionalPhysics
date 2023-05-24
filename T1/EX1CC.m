close all
clear all
clc



%Constnates 
g=9.8;
vlim=6.8;
m=1;
alpha=m*g/(vlim^2);

%vetores
ti=0;
tf=5;
h=0.001;
t=ti:h:tf;
vi_an=-vlim*tanh((g/vlim)*t);
v(1)=16;

z=zeros(length(t),1);
z(1)=1;

%Euler

for i=1:length(t)-1
     a(i)=-g-((g*v(i)*abs(v(i))/vlim^2));    
     v(i+1)=v(i)+a(i)*h;
     z(i+1)=z(i)+v(i)*h;
     
     if z(i+1)<0
         break;
     end
     
end


ft=interp1(z(2072:2073),t(2072:2073),0);


plot(t,z)

