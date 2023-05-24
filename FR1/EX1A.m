close all
clear all
clc

%Constantes

m=1;
k=1;
alpha=0.1;

ti=0;
h=0.1;
tf=20;

%Inicializia??o de Vetores

t=ti:h:tf;
const = [h/2, k*h/(2*m), 2*alpha];
options=optimset('Display','off','Tolx', 1e-10,'TolFun',1e-10);
N=length(t);

x=zeros(1,N);
vx=zeros(1,N);

x(1)=1;
vx(1)=1;




for k=1:N-1
func= @(xv) fcr(xv,x(k),vx(k),const);



xv0=[x(k), vx(k)];
aux=fsolve(func,xv0,options);

x(k+1)=aux(1);
vx(k+1)=aux(2);


end


subplot(3,1,1)
plot(t,x)
title('Posi??o em fun??o do tempo')

subplot(3,1,2)
plot(t,vx)
title('Velocidade em fun??o do tempo')


subplot(3,1,3)
plot(x,vx)
title('Velocidade em fun??o da posi??o')


%Amplitude

A=max(x)


%Periodo
num = 0;
for i=2:N-1

    if (( x(i-1) <= x(i) )&&( x(i) >= x(i+1) ))

        num = num + 1;  % found a maximum
        
        aux = lagr(t(i-1:i+1), x(i-1:i+1));
        tvals(num) = aux(1);
        xvals(num) = aux(2);
        
    end

end

T = ( tvals(end) - tvals(1) ) / (num - 1)















