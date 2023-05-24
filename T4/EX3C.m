close all
clear all
clc


%Constantes 

t0=0;
tf=100;

%Vetores

x0a = 1; x0b = 1.001;
y0a = 2; y0b = 2;
z0a = -1; z0b = -1;




%ODE45
reltol = 3E-14;
abstol_1 = 1E-13;
abstol_2 = 1E-13;
abstol_3 = 1E-13;
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2 abstol_3]);



% solve:
[t, ra] = ode45(@f2,[t0 tf],[x0a y0a z0a],options);
[t, rb] = ode45(@f2,[t0 tf],[x0b y0b z0b],options);


max_index = 100000;   % plot until here ("max_index" must be smaller than the length of both "ra" and "rb")

xa = ra(1:max_index, 1);
ya = ra(1:max_index, 2);
za = ra(1:max_index, 3);

xb = rb(1:max_index, 1);
yb = rb(1:max_index, 2);
zb = rb(1:max_index, 3);

t = t(1:max_index);


% calculate distance between 2 trajectories:
d = sqrt( (xa-xb).^2 + (ya-yb).^2 + (za-zb).^2 );


% plot distance versus time:
figure(1)
semilogy(t, d, 'r-')
hold on


