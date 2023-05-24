close all
clear all
clc

%Constantes
f=1000;

dt=1/f;

y=load("data.txt");
N=length(y);
ft=fft(y);

dw=2*pi/(N*dt);

w=dw*[0:N-1];

w=w- N/2*dw;
ft=fftshift(ft);
dens = (dt*abs(ft)).^2;

f=w/(2*pi);

plot(f,dens);


f_signal = f(dens>0.03);
dens_signal = dens(dens>0.03);

f_noise = f(dens<0.03);
dens_noise = dens(dens<0.03);


figure(2)
plot(f_signal, dens_signal, 'r*')

figure(3)
plot(f, dens, 'r-')
hold on
plot(f_noise, dens_noise, 'b-')


