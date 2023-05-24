% Studying heat conduction using Crank-Nicolson
% (inhomogeneous specific heat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all

% parameters:
k = 0.93;
c1 = 0.094;
c2 = 0.188;
%c2 = 10.0;
ro = 8.9;

L = 50;
dx = 0.5;
dt = 0.1;
tf = 500;

eta1 = k * dt / ( c1 * ro * dx^2 );   % always stable for Crank-Nicolson
eta2 = k * dt / ( c2 * ro * dx^2 );

% vectors:
t = [0:dt:tf]';
x = [0:dx:L]';
nx = length(x);
nt = length(t);
T = zeros(nx,nt);

% boundary conditions:
T(1,:) = 0;
T(nx,:) = 0;
%T(1,:) = 50;
%T(nx,:) = 50;

% initial condition:
T(2:nx-1,1) = 50 * sin(2 * pi * x(2:nx-1) / L);
%T(2:nx-1,1) = 0;


% write matrix A:
n_mat = nx - 2;   % rank of matrix
n_mat_mid = ceil(n_mat / 2);

diagonal(1:n_mat_mid) = 2/eta1 + 2;
diagonal(n_mat_mid+1:n_mat) = 2/eta2 + 2;
A = diag(diagonal) - diag(ones(n_mat-1,1),-1) - diag(ones(n_mat-1,1),+1);

% do LU factorization of matrix A (for Problem 7.2c):
[LL,U,P] = lu(A);   % ('LL' because 'L' is already used...)

% define vector b:
b = zeros(nx-2,1);

% define factors I will need for b vector:
aux_b = zeros(nx-2,1); % Para ser coluna
aux_b(1:n_mat_mid) = 2/eta1 - 2;
aux_b(n_mat_mid+1:n_mat) = 2/eta2 - 2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve using Crank-Nicolson:


tic
% loop over 't' values:
for n=1:nt-1
    
    b = T(1:nx-2,n) + aux_b .* T(2:nx-1,n) + T(3:nx,n);
    b(1) = b(1) + T(1,n+1);
    b(n_mat) = b(n_mat) + T(nx,n+1);
    
    % different possible solutions:-------------
    % 7.2a:
    %T(2:nx-1,n+1) = linsolve(A,b);
    
    % 7.2b:
    T(2:nx-1,n+1) = sol_sist_trid(A,b);
    
    % 7.2c:
    %y = LL \ b; 
    %T(2:nx-1,n+1) = U \ y; 
    
    
    % Alt1:
    %T(2:nx-1,n+1) = inv(A) * b;

    % Alt2:
    %T(2:nx-1,n+1) = A \ b;
    %-------------------------------------------
    
    
    
    
end
toc    % prints the time since "toc"...


% plot results:
figure(11)
contourf(x,t,T')
colorbar
xlabel('x')
ylabel('t')

figure(12)
mesh(x,t,T')
%set(gca,'YDir','reverse');
xlim([0,L])
ylim([0,tf])
zlim([-50,50])
xlabel('x')
ylabel('t')
zlabel('T')







