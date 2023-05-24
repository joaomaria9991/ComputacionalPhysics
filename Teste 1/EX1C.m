%João Maria Machado, NMEC-89132, PL4
close all
clear all
clc

%Constantes

m=1;
K=1;
alpha=0.2;
F0=0;
w0=2;

%Inicializar Vetores

ti=0;
tf=150;
h=0.01;

y0=1.5;
v0=0;


%Parametros da função ODE45


reltol = 3*10^-14;
abstol_1=1*10^-13;
abstol_2=1*10^-13;

h=0.05;
niu_vals=0:h:0.8; %Vetor com os valores de niu


options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);

for i=1:length(niu_vals)
    niu=niu_vals(i);
    
    
    [t,y] = ode45(@f,[ti tf],[y0 v0],options,F0,w0,niu,alpha,K);
    
    x = y(:,1);
    v = y(:,2);
        
    count = 0;
    for n = 2:length(x) - 1
        if x(n - 1) <= x(n) && x(n) >= x(n + 1)  %maximos 
            count = count + 1; %indices dos maximos
            ixd(count)=n;   
        end
    end
    
    aux=lagr(t(ixd(1)-1:ixd(1)+1),x(ixd(1)-1:ixd(1)+1));
    x_max1 = aux(2);
    t_x_max1 = aux(1);
    
    aux=lagr(t(ixd(2)-1:ixd(2)+1),x(ixd(2)-1:ixd(2)+1)); 
    x_max2 = aux(2); %Unecessary...
    t_x_max2 = aux(1);
    
    %In every iteration I retrive 2 max values for y and t, which is all I
    %need to compute the period and range of motion
    
    T (i)= t_x_max2 - t_x_max1;
    A(i)=x_max1;
    
    
    clear x v y t %Clear as recomended
end


subplot(2,1,2)
plot(niu_vals,A,'-o');
title('Variation of Amplitude with {\mu}')

subplot(2,1,1)
plot(niu_vals,T,'-o');
title('Variation of Period with {\mu}')








