close all
clear all
clc

%Constantes 

k=1;
m=1;
ti=0;
tf=10;
h=0.001;
w=1;

%Vetores e condições iniciais

t=ti:h:tf;
x=zeros(1,length(t));

x(1)=1;
v(1)=0;
%Integrador Crank-Nicolson

for i=1:length(t)-1
       a=-k*x(i);
       v(i+1)=v(i)+a*h;
       x(i+1)=v(i+1)*h+x(i);
       EM2(i+1)=1/2*v(i)^2+1/2*k*x(i)^2;

      
end

EM4=1/2*v(end)^2+1/2*k*x(end)^2;

plot(t,x);
title('Método de Crank-Nicolson')


A=lagr(x,t); %Amplitude

%8819
T=t(8819)-t(1);

