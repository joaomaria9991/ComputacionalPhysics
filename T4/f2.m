function drdt = f2(t,r)



% physical parameters:

c=5;

drdt = zeros(3, 1); 



drdt(1) = -r(2)-r(3);  
drdt(2) = r(1)+0.2*r(2);
drdt(3) = 0.2+(r(1)-c)*r(3); 