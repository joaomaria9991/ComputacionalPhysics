close all
clear all
clc



%Constnates 
g=9.8;
vlim=6.8;
m=1;
alpha=m*g/(vlim^2);

%vetores
ti=0;
tf=0.5;
vi_an=-vlim*tanh((g/vlim)*i);
hvals=[0.1 0.01 0.001 0.0001];
vf_vals=zeros(length(hvals),1);
err_vals=zeros(length(hvals),1);


%Integrador

for i=1:length(hvals)

h=hvals(i);
t=ti:h:tf;
v=zeros(length(t),1);

%Euler

for j=1:length(t)-1
     a(j)=-g-((g*v(j)*abs(v(j))/vlim^2));    
     v(j+1)=v(j)+a(j)*h;
     


end
    vi_an=-vlim*tanh((g/vlim)*tf);
    vf_vals=v(end);
    err_vals(i)=abs(v(end) -vi_an(end));

end
loghvals=log(hvals);
logerrvals=log(err_vals);


figure(1)
plot(loghvals,logerrvals);
lsline
hold on

p=polyfit(loghvals',logerrvals,1);

