close all
clear all
clc

%Constantes

k=16;
m=1;
w=sqrt(k/m);


%Vetores e condições iniciais
ti=0;
tf=15;  %500; 12000;

hvals=[0.0001, 0.001, 0.01, 0.1];



x(1)=1;
v(1)=0;




%Funções Implicitas
fv= @(x) -k*x/m;
fx= @(v) v;


%Integrador Runge-Kutta 2ª Ordem


for k=1:length(hvals)
    h=hvals(k);
    clear t v x
    
    
    t=ti:h:tf;
   
    
    
    x=zeros(length(t),1);
    v=zeros(length(t),1);
    
    x(1)=1;
    v(1)=0;
    
    
     
    %solução analitica
    x_an=x(1)*cos(w*t);
    v_an=-x(1)*w*sin(w*t);
    
    
    
    
    for i=1:length(t)-1
        r1v=fv(x(i));
        r1x=fx(v(i));
        
        r2v=fv(x(i)+r1x*(h/2));
        r2x=fx(v(i)+r1v*(h/2));
        
        r3v=fv(x(i)+(1/2)*r2x*h);
        r3x=fx(v(i)+(1/2)*r2v*h);
        
        r4v=fv(x(i)+r3x*h);
        r4x=fx(v(i)+r3v*h);
        
        
        v(i+1)=v(i)+(1/6)*(r1v+2*r2v+2*r3v +r4v)*h;
        x(i+1)=x(i)+(1/6)*(r1x+2*r2x+2*r3x +r4x)*h;
        
        Em(i+1)=0.5*m*abs(v(i))^2+0.5*k*x(i)^2;
    end
    
    errx(k)=abs(x_an(end)-x(end));
    errv(k)=abs(v_an(end)-v(end));
    
    
end

log_hvals=log(hvals);
log_errx=log(errx);
log_errv=log(errv);

subplot(1,2,1)
plot(log_hvals,log_errx)
title('Logaritmo dos Erros-Posição')


subplot(1,2,2)
plot(log_hvals,log_errv)
title('Logaritmo dos Erros-Velocidade')



px=polyfit(log_hvals,log_errx,1);
px(1)


pv=polyfit(log_hvals,log_errv,1);
pv(1)




