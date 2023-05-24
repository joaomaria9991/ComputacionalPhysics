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
tf=2;
h=0.1;
t=ti:h:tf;
vi_an=-vlim*tanh((g/vlim)*t);
v=zeros(1,length(t));


%Euler

for i=1:length(t)-1
     a(i)=-g-((g*v(i)*abs(v(i))/vlim^2));    
     v(i+1)=v(i)+a(i)*h;
     
end



plot(t,v)

