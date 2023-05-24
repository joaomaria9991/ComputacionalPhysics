% Using ode45 to study damping effect in oscillator
% (no forcing here, so I should set FD = 0 in "func.m")
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
num = 0;
for i=2:N-1

    if (( x(i-1) <= x(i) )&&( x(i) >= x(i+1) ))

        num = num + 1;  % found a maximum
        
        % more precise estimate for the time of maximum:
        aux = lagr(t(i-1:i+1), x(i-1:i+1));
        tvals(num) = aux(1);
        xvals(num) = aux(2);
        
    end

end

% get logarithm of "xvals":
log_xvals = log(xvals);

% plot the log of sizes of maxima vs the time when they happened:
% (we should see a straight line, because damping causes an exponential decay of the maxima)
figure(4)
plot(tvals, log_xvals, 'ks-')
hold on

% use polyfit to extract the slope of the straight line:
p = polyfit(tvals, log_xvals, 1);
slope = p(1);
% the slope is given by "q/2" (the rate of decay of the oscillation)
% this can be confirmed by assuming the following form for x:
% x(t) = a * exp(-b * t) * sin(c * t)
% where a,b,c are some constants.

slope=real(slope)









