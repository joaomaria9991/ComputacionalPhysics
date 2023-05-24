close all
clear all
clc

%variáveis 

dx = 0.05;
Nx = 1024;
x = dx * [0:Nx-1];
x = x - Nx/2 * dx;

%Condições iniciais

z0 = 0;
zf = 4;
dz = 0.02;
z = [z0:dz:zf];
Nz = length(z);

q = zeros(Nz, Nx);

q(1,:)=exp(-x.^2./2);



dk = 2 * pi / (Nx * dx);
k = dk * [0:Nx-1];



k = k - Nx/2 * dk;
qt0 = fft( q(1,:) );
qt = zeros(Nz, Nx);



for n1 = 1:Nz
    for n2 = 1:Nx
        qt(n1,n2) = qt0(n2) * exp(-i/2 * k(n2)^2 * z(n1));
    end
end

for n1 = 1:Nz
    qt(n1,:) = ifftshift( qt(n1,:) );
    q(n1,:) = ifft( qt(n1,:) );
end



figure(1)
q2 = abs(q).^2;
mesh(x,z,q2)
axis([-Inf, Inf, -Inf,Inf,-Inf,Inf])
set(gca,'YDir','reverse');
xlabel('x')
ylabel('z')
zlabel('|q|^2')


