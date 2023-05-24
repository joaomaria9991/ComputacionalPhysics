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
eta = k * dt / ( c * p * dx^2 );

if eta>1/2
   'instavel' 
end



for n=1:Nt-1
    for i=2:Nx-1
        T(i,n+1) = T(i,n) + eta * (T(i+1,n) - 2 * T(i,n) + T(i-1,n));
    end
end


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



ind = ceil(Nx/4);    % index of point at 1/4 length
t_ind = interp1(T(ind, ind:end), t(ind:end), 50, 'linear')

% plot to confirm:
figure(13)
plot(t,T(ind,:))
xlabel('t')
ylabel('T (at 1/4 length)')


