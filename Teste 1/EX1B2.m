%João Maria Machado, NMEC-89132, PL4
close all
clear all
clc

%Constantes

m=1;
K=1;
alpha=0.2;
niu=0.8;
F0=0.8;
w0=2;




%Inicializar Vetores

ti=0;
tf=150;
h=0.01;

y0=1.5;
v0=0;


%Parametros da função ODE45


reltol = 3*10^-14;
abstol_1=1*10^-13;
abstol_2=1*10^-13;


options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);
[t,y] = ode45(@f,[ti tf],[y0 v0],options,F0,w0,niu,alpha,K);


x = y(:,1);
v = y(:,2);

figure(1)
subplot(3,1,1)
plot(t,x)
title('Posição do Oscilador em função do Tempo')


subplot(3,1,2)
plot(t,v)
title('Velocidade do Oscilador em função do Tempo')



subplot(3,1,3)
plot(x,v)
title('Position-Velocity graph')





