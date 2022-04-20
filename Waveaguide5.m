re = 2.818e-15;
lambda = 3e-6;
V = 2.405;
k = 2*pi/lambda;
a = 10;

nc = pi/(re*lambda^2);

ne = 2*nc*(-1+sqrt(1+V^2/(k^2*a^2)));

v = 0.2*3*10^8;
e = 1.602e-19;
J = ne*v*e;
disp(J)


%%Stability
eps0 = 8.85418e-12;
me = 9.109e-31;
gamma = 1/sqrt(1-0.2^2);
meprime = me*gamma;
wp = sqrt(ne*e^2/(meprime*eps0));

n = 1+ne/(2*nc);
dnsqd = (n^2-1^2)/n^2;

Dtheta = (pi/12)*(wp/(k*lambda)^2)*(k/k)*(dnsqd);

disp(Dtheta)
