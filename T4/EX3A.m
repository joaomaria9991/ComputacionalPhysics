close all
clear all
clc


%Constantes 

t0=0;
tf=100;

%Vetores

x0=0;
y0=2;
z0=-1;




%ODE45
reltol = 3E-14;
abstol_1 = 1E-13;
abstol_2 = 1E-13;
abstol_3 = 1E-13;
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2 abstol_3]);



[t, r] = ode45(@f2,[t0 tf],[x0 y0 z0],options);

min_index = 60000;   % start plotting from here to see limit cycles clearly

x = r(min_index:end, 1);
y = r(min_index:end, 2);
z = r(min_index:end, 3);


% do plotting:
figure(1)
plot3(x, y, z, 'r-')
hold on



