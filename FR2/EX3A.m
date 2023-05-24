close all
clear all
clc

%Constantes
alpha=0.2;


%Condições fronteira


%Inicialização de Vetores

dx=0.1;
dz=0.1;

xi=-20;
xf=20;

zi=0;
zf=16;

x=xi:dx:xf;
z=zi:dz:zf;

Nx=length(x);
Nz=length(z);


% define 2D vector for solution:
Phi = zeros(Nx, Nz);

% set initial and boundary conditions:
Phi(1:Nx,1) = exp(-x.^2/2);
Phi(1,:) = 0; 
Phi(Nx,:) = 0;




%Estabilidade
eta=(1i*dz)/(4*(dx)^2);
qsi=2*alpha*(dx)^2;


aux1 = 1/eta + 2 + qsi * (x(2:Nx-1)).^2;
aux2 = 1/eta - 2 - qsi * (x(2:Nx-1)).^2;
aux2 = transpose(aux2);



%Construir a Matriz A
NM=Nx-2;
A = diag(aux1) - diag(ones(NM-1,1),-1) - diag(ones(NM-1,1),+1);




%Construir o vetor b

b=zeros(Nx-2,1);


for n=1:Nz-1
    
    b = Phi(1:Nx-2,n) + aux2 .* Phi(2:Nx-1,n) + Phi(3:Nx,n);
    b(1) = b(1) + Phi(1,n+1);
    b(NM) = b(NM) + Phi(Nx,n+1);
    Phi(2:Nx-1,n+1) = linsolve(A,b);
   

end

aPhi = abs(Phi);




figure(1)
contourf(x,z,aPhi')
colorbar
xlabel('x')
ylabel('z')

figure(2)
mesh(x,z,aPhi')
set(gca,'YDir','reverse');

xlabel('x')
ylabel('z')
zlabel('|\Phi|')


