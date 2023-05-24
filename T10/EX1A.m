close all 
clear all
clc

rng(13)

R=1;
a=0;
b=1;
N=1000;

x=rand(N,1);

f=sqrt(1-x.^2);

med=(b-a)*mean(f);
exact=R^2*pi/4;