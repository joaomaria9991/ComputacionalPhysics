% Using ode45 to study damped and forced oscillator
% (We can use a lot of the program for Problem 3.3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

% parameters:
x0 = 0.2;
v0 = 0;
% for simplicity:
% I'm using x for "theta"
% and v for "the derivative of theta"

t0 = 0;
tf = 500;

% the physical parameters are defined in the function "func.m"

% define options (same as in Problem 3.3):
reltol = 3E-14;
abstol_1 = 1E-13;
abstol_2 = 1E-13;
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);

% solve:
[t, y] = ode45(@func,[t0 tf],[x0 v0],options);

x = y(:,1);
v = y(:,2);


% do plotting:
figure(1)
plot(t, x, 'r-')

figure(2)
plot(t, v, 'b-')

figure(3)
plot(x, v, 'g-')