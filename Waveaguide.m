a= 10;
lambda= 3e-6;
n2 = 1;
delta = 1e-18;
n1 = n2+delta;

V = 2*pi*a/lambda*sqrt(n1^2-n2^2);

MFR = a*(0.65+1.619*V^(-3/2)+2.879*V^(-6));

%For single-mode optics we need V as a function of n1,lambda,a
figure;
subplot(1,2,1)



delta = (1e-18)*logspace(1,4,100);
n1 = n2 + delta;
V = (2*pi*a/lambda).*sqrt(n1.^2-n2.^2);


plot(n1-n2,V,"-k");
hold on
a=100;
V = (2*pi*a/lambda).*sqrt(n1.^2-n2.^2);
plot(n1-n2,V,"-b");
yline(2.405,"-r","Single-Moded Cutoff")
legend("V at = 10m","V at a= 100m")
title("V against {\delta n}")


subplot(1,2,2)
a=10;
delta = (1e-18)*logspace(1,4,100);
n1 = n2 + delta;
V = (2*pi*a/lambda).*sqrt(n1.^2-n2.^2);
MFR = (0.65+1.619.*V.^(-3/2)+2.879.*V.^(-6));

loglog(n1-n2,MFR,"-k");
hold on
a=100;
V = (2*pi*a/lambda).*sqrt(n1.^2-n2.^2);
MFR = (0.65+1.619.*V.^(-3/2)+2.879.*V.^(-6));

loglog(n1-n2,MFR,"-b");
yline(1.1,"-r","Single-Moded Cutoff")
legend("MFR/a at = 10m","MFR/a at a= 100m")
title("MFR/a against {\delta n}")

