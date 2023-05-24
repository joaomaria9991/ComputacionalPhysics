close all 
clear all
clc


reltol = 3*10^-14;
abstol_1=1*10^-13;
abstol_2=1*10^-13;


t0=0;
tf=15;
v0=0;
x0=1;
k=16;
m=1;

options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);

%solving
[t,y] = ode45(@f,[t0 tf],[x0 v0],options);

x = y(:,1);
v = y(:,2);

plot(t,x)

