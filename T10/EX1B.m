close all
clear all
clc

R=1;
a=0;
b=1;
d=5;
N=10000;

x=rand(N,d-1);
x2=x.^2;
sumx2=sum(x2,2);
f=zeros(N,1);
f(sumx2<=1)=sqrt( 1 - sumx2(sumx2 <= 1) );

med=(b-1)^(d-1)*mean(f);
exact=pi^(d/2)/gamma(d/2+1)/2^d*R^d;