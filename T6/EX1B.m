close all
clear all
clc

%Transformadas de Fourier 

dt=0.05;
N=2048;

t=dt*[0:N-1];
dw=2*pi/(N*dt);
w=dw * [0:N-1];

%Ossilações
y3=sin(10*t)+sin(40*t);

%transformadas sem shift 
ft3=fft(y3);

de=(dt*abs(ft3)).^2;

subplot(2,1,1)
plot(w,de)
title('FFT not centered')


%Centrando na origem

w=w- N/2*dw; %N é par
de=fftshift(de);

subplot(2,1,2)
plot(w,de)
title('FFT Centered')

%Temos de aumentar a amostragem quando aumentamos a frequencia da função
%N passa a 0.05 e N=2*N

