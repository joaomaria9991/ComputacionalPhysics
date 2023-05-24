close all
clear all
clc

%Constantes

tol=1*10^-3;
niu=10^-3;
L=1;
T=1000;

%Inicialização de vetores

h=0.01;
x=0:h:L;
N=length(x);


%Condições Iniciais


yf(1) = 1;
j = 0;


w(1)=2000;
w(2) = 3500;



while( abs( yf(end) ) > tol)
    
    j=j+1;
    
    %Condições iniciais
    y(1)=0;
    dy(1)=1;
    
    
    %Integrador Runge-Kutta
    
    
    %Funções Implicitas
    
    fy = @(dy) dy;
    fdy = @(y) - w(j)^2 * niu / T * y;
    
    
    for i=1:N-1
        
        r1dy=fdy(y(i));
        r1y=fy(dy(i));
        
        r2dy=fdy(y(i)+r1y*(h/2));
        r2y=fy(dy(i)+r1dy*(h/2));
        
        r3dy=fdy(y(i)+(h/2)*r2y);
        r3y=fy(dy(i)+(h/2)*r2dy);
        
        r4dy=fdy(y(i)+r3y*h);
        r4y=fy(dy(i)+r3dy*h);
        
        y(i+1) = y(i) + h * (1/6*r1y + 1/3*r2y + 1/3*r3y + 1/6*r4y);
        dy(i+1) = dy(i) + h * (1/6*r1dy + 1/3*r2dy + 1/3*r3dy + 1/6*r4dy);
            
        
    end
    
    %Guardar Valores que interessam
    
    yf(j) = y(end);
    
    if(j>1)
        diff = (yf(j) - yf(j-1)) / (w(j) - w(j-1));
        w(j+1) = w(j) + ( 0 - yf(j) ) / diff;
    end
    
    
end

w1_num = w(end)

%Eu não plotei porque não gosto do plot e sou feio :(

