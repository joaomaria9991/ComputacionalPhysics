close all
clear all
clc

%Constantes
m=0.057;
R=0.0067/2;
A=pi*R^2;
p=1.225;
g=9.8;

w=0;%-3000*((2*pi)/60); %100*pi;




%Vetores
h=0.1;
tf=3;
t=0:h:tf;

z=zeros(1,length(t));
x=zeros(1,length(t));

z(1)=0.7;
vi=-20;

v_x(1)=vi*cos(5);
v_z(1)=vi*sin(5);
v(1)=20;

angle=5*pi/180;
%Integrador - EULER

for i=1:length(t)-1
 S=w*R/v(i);
 CD=0.508+(22.503+4.196*S^-2.5)^-0.4;
 CL=(2.022+0.981*S^-1)^-1;


 a_x(i)=(-0.5*CD*p*(v(i)*cos(angle)).^2*A)/m;
 a_z(i)=-g-(0.5*CD*p*(v(i)*sin(angle)).^2*A)/m-(-0.5*CL*p*A*v(i)^2*(w*v(i)*sin(angle))/m);
 
 v_x(i+1)=v_x(i)+a_x(i)*h;
 v_z(i+1)=v_z(i)+a_z(i)*h;
 
 v(i+1)=sqrt((v_x(i).^2+v_z(i).^2));
 
 z(i+1)=z(i)+v_z(i)*h;
 x(i+1)=x(i)+v_x(i)*h;      
end


plot(x,z)





