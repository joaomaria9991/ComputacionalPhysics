close all
clear all
clc

%Constantes

V0=20;
a=1;
b=4;    %a FUNÇÃO DE ONDA TENDE PARA ZERO NESTA REGIÃO
h=0.001;
tol=1e-7;

%Vetores
x=[-b:h:b];
N=length(x);

%Posição onde os shootings se encontram
x_match=0.7;

%Indice de x_match no vetor x.
ind_match=round((x_match+b)/h+1);

%Dividindo o vetor em 2 para os 2 shoots
x_left=x(1:ind_match);
x_right=x(ind_match:N);
N_left=length(x_left);
N_right=length(x_right);


%Vetores do potêncial

V_left=zeros(1,N_left);
V_left(x_left<-a)=V0;

V_right=zeros(1,N_right);
V_right(x_right>a)=V0;


psi_extremo=0;
psi_next_left=h;
psi_next_right=psi_next_left;

%Primeiras guesses de 'E'

E(1)=0.5;
E(2)=0.6;


%Metodo de Shooting:
reldiff(1)=1;
iE=0;


while(abs(reldiff(end))>tol)
    
    %Integrador do shoot à esquerda
    iE=iE+1;
    g = 2 * (E(iE) - V_left);
    aux1=(1+h^2/12*g);
    aux2=2*(1-5*h^2/12*g);
    
    psi_left=zeros(1,N_left);
    psi_left(1)=psi_extremo;
    psi_left(2)=psi_next_left;
    
    
    for k=2:N_left-1
        psi_left(k+1) = ( -aux1(k-1) * psi_left(k-1) + aux2(k) * psi_left(k) ) / aux1(k+1);
    end
    
    clear g  %Importante, o g pode variar
    
    g=2*(E(iE)-V_right);
    aux1=(1+h^2/12*g);
    aux2=2*(1-5*h^2/12*g);
    
    
    psi_right=zeros(1,N_right);
    psi_right(N_right)=psi_extremo;
    psi_right(N_right-1)=psi_next_right;
    
    
    for k=N_right-1:-1:2
        psi_right(k-1) = ( -aux1(k+1) * psi_right(k+1) + aux2(k) * psi_right(k) ) / aux1(k-1);
    end
    
    
    ratio = psi_left(N_left) / psi_right(1);
    psi_right = psi_right * ratio;

    
    D_left = ( 25/12*psi_left(N_left) - 4*psi_left(N_left-1) + 3*psi_left(N_left-2) - 4/3*psi_left(N_left-3) + 1/4*psi_left(N_left-4) ) / h;
    D_right = ( -25/12*psi_right(1) + 4*psi_right(2) - 3*psi_right(3) + 4/3*psi_right(4) - 1/4*psi_right(5) ) / h;
    
    reldiff(iE) = (D_left - D_right) / (D_left + D_right);
    
    if(iE>1)
        slope = (reldiff(iE) - reldiff(iE-1)) / (E(iE) - E(iE-1));
        E(iE+1) = E(iE) + ( 0 - reldiff(iE) ) / slope;
    end
    
    
end

E_num = E(end)

psi = [psi_left(1:end-1), psi_right];



I = trapz(x, psi.^2);
psi = psi / sqrt(I);



figure(2)
plot(x, psi, 'r-');
title('normalized wave function');
xlim([-b,b]);
psimax = max(psi);
psimin = min(psi);
dpsi = psimax - psimin;
ylim([psimin - dpsi*0.3, psimax + dpsi*0.3]);

