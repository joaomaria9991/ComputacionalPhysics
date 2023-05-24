close all
clear all
clc

%Constantes

m=1.5;
K=2;
alpha=0.2;

tol=10^-4;
ti=0;
tf=20;
h=0.1;


%Condições iniciais

x(1)=1.9;
v(1)=0;

%Options
reltol = 3E-14;
abstol_1 = 1E-13;
abstol_2 = 1E-13;
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);


%Condições de Shooting

alpha(1) = -0.2;
alpha(2) = -0.25;

j=0;


negamp(1) = 0;
dnegamp = -1.5;


while (abs( negamp(end) - dnegamp ) > tol)
    j=j+1;
    
    %Integrador ODE45
       
    [t, y] = ode45(@f1, [ti tf], [x(1) v(1)],options, K, m, alpha(j));
    x = y(:,1);
    v = y(:,2);
    N=length(t);
    
    
    n = 0;
    clear xmin;
    
    %Shooting mesmo
    for k = 2:N-1
        if ( (x(k+1) > x(k)) && (x(k-1) >= x(k)) )
            n = n + 1;
            aux = lagr(t(k-1:k+1),-x(k-1:k+1));
            xmin(n) = -aux(2);
        end
    end
    
    negamp(j) = mean(xmin);
    
    
    % Proximo palpite
    if (j > 1)
        declive = (negamp(j) - negamp(j-1)) / (alpha(j) - alpha(j-1));
        alpha(j+1) = alpha(j) + ( dnegamp - negamp(j) ) / declive;
    end
end


x(end)


