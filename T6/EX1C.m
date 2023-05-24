close all
clear all
clc

%Transformadas de Fourier 

dt=0.1;
N=2^15;

t=dt*[0:N-1];
dw=2*pi/(N*dt);
w=dw * [0:N-1];

%Ossila��es
y3=sin(10*t)+sin(10.5*t);

%transformadas sem shift 
ft3=fft(y3);

de=(dt*abs(ft3)).^2;

subplot(2,1,1)
plot(w,de)
title('FFT not centered')


%Centrando na origem

w=w- N/2*dw; %N � par
de=fftshift(de);

subplot(2,1,2)
plot(w,de)
title('FFT Centered')

%Temos de aumentar N,

