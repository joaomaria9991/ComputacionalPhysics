% Using ode45 to study frequency of forced oscillator
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
tf = 100;

% the physical parameters are defined in the function "func.m"

% define options (same as in Problem 3.3):
reltol = 3E-14;
abstol_1 = 1E-13;
abstol_2 = 1E-13;
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);

% solve:
[t, y] = ode45(@func,[t0 tf],[x0 v0],options);

N = length(t);
x = y(:,1);
v = y(:,2);


% do plotting:
%figure(1)
%plot(t, x, 'r-')

%figure(2)
%plot(t, v, 'b-')

%figure(3)
%plot(x, v, 'g-')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find maxima of x:
% (this is copy-pasted from a previous program)

% I only want to find maxima after t=50, so I need to find the index that
% corresponds to t=50:
for i=2:N-1
    if ( t(i) > 50 )
        break;
    end
end
start_index = i;  % the index of the first time step after t=50

num = 0;
for i=start_index:N-1

    if (( x(i-1) <= x(i) )&&( x(i) >= x(i+1) ))

        num = num + 1;  % found a maximum
        
        % more precise estimate for the time of maximum:
        aux = lagr(t(i-1:i+1), x(i-1:i+1));
        tvals(num) = aux(1);
        xvals(num) = aux(2);
        
    end

end

% determine period:
T = ( tvals(end) - tvals(1) ) / (num - 1);

% determine frequency:
f = 1 / T
% this is more or less the same as "wD/(2*pi)"





