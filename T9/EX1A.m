close all
clear all
clc

%Constntes

a=1;    %tamanho de meio poço
h=0.01; %paço temporal
tol=1e-7;

%Vetores
x=[-a:h:a];
N=length(x);

%Valor próprio 'E':
E_an = 0.5*(pi^2/8*a^2);

%Primeiras guesses de E:

E(1)=1.5;
E(2)=1.6;


%Primeira dois valores de 'psi'

psi(1)=0;
psi(2)=h; %um valor pequeno


%Metodo de Shooting

psif(1)=1;  %um valor maior que toler para entrar no ciclo  while
%vai ser um vetor que contem os valore de psi(end)

iE=0;   %variavel de iteração


while(abs(psif(end)-0)>tol) %a condição de paragem
    iE=iE+1;
    
    g=2*E(iE);
    aux=(h^2/12)*g;
    
    
    for i=2:N-1
        psi(i+1)=(1+aux)^(-1)*(-(1+aux)*psi(i-1)+2*(1-5*aux)*psi(i));
    end
    
    %Guardar o ultimo elemnto
    psif(iE)=psi(end);
    
    
    %Proximas guesses com método de secante
    
    
    if(iE>1)
        m=(psif(iE)-psif(iE-1))/(E(iE)-E(iE-1));
        E(iE+1)=E(iE)+(0-psif(iE))/m;
        
    end
end


E_num=E(end);


%Normalizar a função de onda 

I=trapz(x,psi.^2);
psi=psi/sqrt(I);




plot(x, psi, 'r-');
title('normalized wave function');
xlim([-a,a]);
psimax = max(psi);
psimin = min(psi);
dpsi = psimax - psimin;
ylim([psimin - dpsi*0.3, psimax + dpsi*0.3]);





    
    
    
    
