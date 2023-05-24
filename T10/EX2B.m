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
d=4;

exact=pi^(d/2)/gamma(d/2+1)/2^d*R^d; %Valor analitico

for i=1:length(Nvals)   %Percorre o array de vaslores de N
    
    N=Nvals(i);
    clear x f diff2
    
    
    for j=1:n
        x=rand(N,d-1);
        x2=x.^2;
        sumx2=sum(x2,2);
        f=zeros(N,1);
        f(sumx2<=1)=sqrt(1-sumx2(sumx2<=1));
        
        med=(b-a)^(d-1)*mean(f);
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





