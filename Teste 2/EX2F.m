%João Maria Machado,NMEC-89132, PL4

close
clear
clc


%Constantes

fator_res=4;

m=1;
K=1;
alpha=0.1;
niu=0.5; %from last exercise

tol=10^-5;
ti=0;
tf=102.3;
h=0.1*fator_res;
t=0:h:tf;
N=length(t);

%Condições iniciais

x(1)=0;
v(1)=1.5;

%Options
reltol = 3E-14;
abstol_1 = 1E-13;
abstol_2 = 1E-13;
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);

[t, y] = ode45(@f1, t, [x(1) v(1)],options, K, m, niu,alpha);
x = y(:,1);
v = y(:,2);

dw=2*pi/(N*h);
w=dw*[0:N-1];
w=w- N/2*dw;

ft=fft(x);
ft=fftshift(ft);
dens = (h*abs(ft)).^2;

dens_signal1 = dens(dens<10);
w_signal1=w(dens<10);

% subplot(2,1,1)
% plot(w_signal,dens_signal)
% title('FFT with no amplification')
% subplot(2,1,2)
plot(w_signal1,dens_signal1)
title('FFT with 4x amplification')
