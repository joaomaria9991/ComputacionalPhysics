%João Maria Machado,NMEC-89132, PL4

close all
clear all
clc

%Constantes

m=1;
K=1;
alpha=0.1;

tol=10^-6;
ti=0;
tf=60;
h=0.1;
t=0:h:tf;

%Condições iniciais

x(1)=0;
v(1)=1.5;

%Options
reltol = 3E-14;
abstol_1 = 1E-13;
abstol_2 = 1E-13;
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);


%Condições de Shooting

niu(1) = 0.35;
niu(2) = 0.4;

j=0;


negamp(1) = 0;  
dnegamp = 2.15;


while (abs( negamp(end) - dnegamp ) > tol)
    j=j+1;
    
    %Integrador ODE45
    
    [t, y] = ode45(@f1,t, [x(1) v(1)],options, K, m, niu(j),alpha);
    x = y(:,1);
    v = y(:,2);
    N=length(t);
    
    
    n = 0;
    
    %Shooting for the amplitude
    clear xmin
    for k = 2:N-1
        if ( (x(k+1) < x(k)) && (x(k-1) < x(k)) )
            n = n + 1;
            aux = lagr(t(k-1:k+1),x(k-1:k+1));
            xmax(n) = aux(2);
        end
    end
    
    negamp(j) = mean(xmax);
    
    
    % Next Guess
    if (j > 1)
        declive = (negamp(j) - negamp(j-1)) / (niu(j) - niu(j-1));
        niu(j+1) = niu(j) + ( dnegamp - negamp(j) ) / declive;
    end
    
    
    
end



niu
xmax


