format long    

lambda = 3e-6;
a = 10;
n2 = 1;
% deltan = 1;
% n1 = n2 + deltan;
n1 = sqrt((2.405*lambda/(2*pi*a))^2+1);
deltan = n1-n2;


V = 2*pi*a*sqrt(n1^2-n2^2)/lambda;

mfra = (0.65+1.619/(V^(3/2))+2.879/V^(6));

disp(V)
disp(mfra)
