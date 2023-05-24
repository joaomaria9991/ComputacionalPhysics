close all
clear all
clc

%Constantes
rmax=50;
h=0.001;
tol=1e-7;

n=3;
l=1;


%Vetores
r=[0:h:rmax];
N=length(r);

u=zeros(1,N);
g=zeros(1,N);

%Valor proprio E

E_an = -1/2 * n^(-2)


%Primeiras guesses do m�todo de Shooting



E(1) = -0.6/9;
E(2) = -0.7/9;


u(N)=0;
u(N-1)=h/1000;

uf(1)=1; %Mais uma vez serve para entrar no ciclo while


iE=0; %variavel de itera��o

while (abs(uf(end)-0)>tol)
    
    iE=iE+1;
    
    %M�todo de Numerov
    for k=2:N   %Atendento aos valores pr�rpios
        g(k) = 2*E(iE) + 2/r(k) - l*(l+1)/(r(k)^2);
        
        
        
    end
    for k=N-1:-1:3
        u(k-1) = (1+h^2*g(k-1)/12)^(-1)*(2*(1-5*h^2*g(k)/12)*u(k)-(1+h^2*g(k+1)/12)*u(k+1));
        
    end
    
    u(1)=interp1(r(2:10), u(2:10), 0, 'spline');
    
    
    %Guardar o primeiro elemnto
    uf(iE) = u(1); %#ok<SAGROW>
    
    if(iE>1)
        m=(uf(iE)-uf(iE-1))/(E(iE)-E(iE-1));
        E(iE+1)=E(iE)+(0-uf(iE))/m;
        
    end
    
end
E_num=E(end);



%Normalizar a fun��o de onda
I = trapz(r, u.^2);
u = u / sqrt(I);


R(2:N) = u(2:N)./r(2:N);
R(1) = interp1(r(2:5), R(2:5), 0, 'spline');

%Plots

figure(2)
plot(r, R, 'r-');
title('normalized radial wave function');
xlim([0,rmax]);
Rmax = max(R);
Rmin = min(R);
dR = Rmax - Rmin;
ylim([Rmin - dR*0.3, Rmax + dR*0.3]);

hold on

% plot also exact form:
R_exact = -4*sqrt(2) / (27*sqrt(3)) * (1-r/6).*r.*exp(-r/3);   % n=3, l=1
plot(r,R_exact,'k--');

