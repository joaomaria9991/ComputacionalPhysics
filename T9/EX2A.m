close all
clear all
clc

%Constantes
rmax=20;
h=0.001;
tol=1e-7;

n=1;
l=0;


%Vetores
r=[0:h:rmax];
N=length(r);

u=zeros(1,N);
g=zeros(1,N);

%Valor proprio E

E_an=-0.5*N^(-2);


%Primeiras guesses do método de Shooting

E(1)=-0.6;
E(2)=-0.7;

u(N)=0;
u(N-1)=h/1000;

uf(1)=1; %Mais uma vez serve para entrar no ciclo while


iE=0; %variavel de iteração

while (abs(uf(end)-0)>tol)
    
    iE=iE+1;
    
    %Método de Numerov
    for k=2:N   %Atendento aos valores prórpios
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



%Normalizar a função de onda
I = trapz(r, u.^2);
u = u / sqrt(I);


R(2:N) = u(2:N)./r(2:N);
R(1) = interp1(r(2:5), R(2:5), 0, 'spline');

%Plots

figure(2)
plot(r, R, 'r-');
title('Função de Onda normalizada');
xlim([0,rmax]);
Rmax = max(R);
Rmin = min(R);
dR = Rmax - Rmin;
ylim([Rmin - dR*0.3, Rmax + dR*0.3]);
hold on
R_exact = 2.*exp(-r);   % for n=1, l=0
plot(r,R_exact,'k--');


