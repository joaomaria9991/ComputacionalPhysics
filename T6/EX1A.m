close all
clear all
clc

%Transformadas de Fourier 

dt=0.1;
t=0:dt:102.3;
N=1024;

dw=2*pi/(N*dt);
w=dw * [0:N-1];

%Oscilações

y1=(sin(t));

y2=(sin(10*t));

y3=y1+y2;

%transformadas sem shift 

ft1=fft(y1);
ft2=fft(y2);
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




