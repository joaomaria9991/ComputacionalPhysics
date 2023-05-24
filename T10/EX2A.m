close all
clear all
clc

%Seeds do RNG
rng(89132)

%Constantes
R=1;
a=0;
b=1;
Nvals=[10,100,1000,10000];
n=1000;

exact=R^2*pi/4; %Valor analitico

for i=1:length(Nvals)   %Percorre o array de vaslores de N
    
    N=Nvals(i);
    clear x f diff2
    
    
    for j=1:n
        
        x=rand(1,N);
        f=sqrt(1-x.^2);
        
        med=(b-a)*mean(f);
        diff2(j)=(med-exact)^2;
        
    end
    
    errvals(i)=sqrt(mean(diff2));
    
end

idx=logspace( log10(min(Nvals)), log10(max(Nvals)), 100 );
err_th = (b-a) * std(f) * idx.^(-1/2);



figure(1)
loglog(Nvals, errvals, 'rs');
ylabel('erro')
xlabel('tamanho da amostra')
title('Variação do erro com  o tamanho da amostra')
hold on
loglog(idx, err_th, 'k-');





