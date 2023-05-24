close all
clear all
clc

%Constantes


niu=10^-3;
L=1;
T=1e3;

h=0.01;
x=[0:h:L]';
N=length(x);
n=N-2;


%Gerar Matriz A

A=eye(n,n);
A=-2*A;

for i=1:n
    A(i,i+1) = 1;   
    A(i+1,i) = 1;
    
end

[vec,val] = eigs(A, 3, 'sm');
vals = diag(val);

omega_vals = sqrt(-vals * T / niu) / h


m1=vec(:,1)';
m2=vec(:,2)';
m3=vec(:,3)';

subplot(3,1,1)
plot(x(1:end-1),m1)
title('1� Modo de Vibra��o')

subplot(3,1,2)
plot(x(1:end-1),m2)
title('2� Mode de Vibra��o')

subplot(3,1,3)
plot(x(1:end-1),m3)
title('3� Modo de Vibra��o')



