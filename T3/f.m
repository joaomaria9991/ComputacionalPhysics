function derivadas = f(t,y)
k=16;
m=1;

derivadas = zeros(2,1);
derivadas(1) =y(2);
derivadas(2) =-k/m*y(1);
