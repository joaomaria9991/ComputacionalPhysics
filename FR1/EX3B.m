% Badminton - 3rd order Runge-Kutta method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

% define parameters:
g = 9.8;
v_lim = 6.8;
m = 1;
alpha = m*g / v_lim^2;
z0 = 1;
v0 = 16;

% set variables:
t0 = 0;
tf = 2;
h = 0.01;

% define vectors:
t = [t0:h:tf]';
N = length(t);
z = zeros(N, 1);
v = zeros(N, 1);

% initial conditions:
z(1) = 1;
v(1) = 16;


% Run 3rd order Runge-Kutta:
for i=1:N-1
    
    % the derivatives of v and z at time t = t(i):
    r1v = -g - alpha * sign(v(i)) * (abs(v(i)))^2;
    r1z = v(i);
    
    % the derivatives of v and z at time t = t(i) + h/2:
    r2v = -g - alpha * sign(   v(i) + h/2 * r1v    ) * (abs(    v(i) + h/2 * r1v   ))^2;
    r2z = v(i) + h/2 * r1v;
    
    % the derivatives of v and z at time t = t(i) + h:
    r3v = -g - alpha * sign( v(i) - h * r1v + 2*h * r2v ) * (abs( v(i) - h * r1v + 2*h * r2v ))^2;
    r3z = v(i) - h * r1v + 2*h * r2v;

    
    % calculate v(i+1) and z(i+1) using r1v, r2v and r1z, r2z:
    v(i+1) = v(i) + h * (1/6 * r1v + 2/3 * r2v + 1/6 * r3v);
    z(i+1) = z(i) + h * (1/6 * r1z + 2/3 * r2z + 1/6 * r3z);
        
end

% do plotting:
subplot(2,1,1)
plot(t, z, 'r-');
hold on
xlabel('time');
ylabel('altitude');

subplot(2,1,2)
plot(t, v, 'r-');
hold on
xlabel('time');
ylabel('velocity');






