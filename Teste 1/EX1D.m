%João Maria Machado, NMEC-89132, PL4
close all
clear all
clc

%Constantes

m=1;
K=1;
alpha=0.2;
niu=0.8;
F0=0;
w0=1;

%Inicializar Vetores

ti=0;
tf=150;
h=0.0001;
t=ti:h:tf;

y=zeros(length(t),1);
v=zeros(length(t),1);

y(1)=1.5;
v(1)=0;

%Funções Implicitas
fy= @(v) v;
fv= @(y,v,t) (niu*cos(v)*v+F0*cos(w0*t)-K*(y+alpha*y^3))/m;


for i=1:length(t)-1
    r1v=fv(y(i),v(i),t(i));
    r1y=fy(v(i));
    
    r2v=fv(y(i)+r1y*(h/2),v(i)+r1y*(h/2),t(i));
    r2y=fy(v(i)+r1v*(h/2));
    
    r3v=fv(y(i)+r2y*(h/2),v(i)+r2v*(h/2),t(i));
    r3y=fy(v(i)+r2v*(h/2));
    
    
    r4v=fv(y(i)+r3y*h,v(i)+r3v*h,t(i));
    r4y=fy(v(i)+r3v*h);
    
    v(i+1)=v(i)+(h/6)*(r1v+2*r2v+2*r3v+r4v);
    y(i+1)=y(i)+(h/6)*(r1y+2*r2y+2*r3y+r4y);
    
end

% figure(1)
% subplot(3,1,1)
% plot(t,y)
% title('Posição do Oscilador em função do Tempo')
% 
% 
% subplot(3,1,2)
% plot(t,v)
% title('Velocidade do Oscilador em função do Tempo')
% 
% 
% subplot(3,1,3)
plot(y,v)
title('Position-Velocity graph')

