close all
clear all
clc

%Oscilador Van der Pol

ti=0;
tf=100;
 

%Condi??es iniciais 

y0=2;
v0=0;


%ODE45


reltol = 3E-14;
abstol_1 = 1E-13;
abstol_2 = 1E-13;
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);

[t, x] = ode45(@f,[ti tf],[y0 v0],options);

y = x(:,1);
v = x(:,2);


plot(y, v, 'b-')
hold on



