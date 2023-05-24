%alocar memória
close all
clear all
clc

%Constantes

tf=500;
L=50;
k=0.93;
c=0.094;
p=8.98;


%Inicialização de Vetores

dx=0.5;
dt=0.1;

x=0:dx:L;
t=0:dt:tf;

Nx=length(x);
Nt=length(t);


T=zeros(Nx,Nt);


%Condições iniciais

T(1,:) = 0;
T(Nx,:) = 0;
T(2:Nx-1,1) = 100;

%Estabilidade

n=(k*dt)/(c*p*(dx)^2);

%Construir a Matriz A
NM=Nx-2;
D1=((2/n)+2);
A=eye(NM);

A=D1*A;
A(1,2)=-1;

for i=2:NM-1
    A(i,i-1)=-1;
    A(i,i+1)=-1;
end

A(NM,NM-1)=-1;

%Construir o vetor B

b=zeros(Nx-2,1);
D2=((2/n)-2);

%Integrar

for n=1:Nt-1
    for i=1:Nx-2
        b(i) = T(i,n)+D2*T(i+1,n)+T(i+2,n);
    end
    
    b(1)=b(1)+T(1,n+1); 
    b(NM)=b(NM)+T(Nx,n+1);
    T(2:Nx-1,n+1)=linsolve(A,b);
end



% plot results:
figure(11)
contourf(x,t,T')
colorbar
xlabel('x')
ylabel('t')

figure(12)
mesh(x,t,T')
set(gca,'YDir','reverse');
xlim([0,L])
ylim([0,tf])
zlim([0,100])
xlabel('x')
ylabel('t')
zlabel('T')




