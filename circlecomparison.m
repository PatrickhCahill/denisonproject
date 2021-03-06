R = 5;
sym r

psteady = @(r) (1./r).*log((R+r).^2./(R-r).^2);
psteadyo = @(x,y,r) 1./(R^2+x.^2 + y.^2 - 2.*sqrt(x.^2+y.^2).*r);
psteady2 = @(x,y) 2*integral(@(r) psteadyo(x,y,r),-R,R);
p2s = [];
r = linspace(0,12,100);
for r=r
    p2s = [p2s, psteady2(r,0)];
end
r = linspace(0,12,100);
ps = psteady(r);
figure
plot(r,ps,"-r")
hold on
plot(r,p2s,"bo")
xlabel("Distance of ray from origin (m)")
ylabel("Number density")
title("Steady State, numerical and anylitical comparison of number density in yz plane")


pboltzmann = @(r) 2*N./(t.*r*beta0)*(beta0/pi)^(3/2).*exp((-beta0/t^2).*(R^2+r.^2)).*sinh((2*beta0*R/t^2).*r);

pboltzo = @(x,y,r) (N/t^3)*(beta0/pi)^(3/2)*exp(-beta0*(R^2+x.^2 + y.^2 - 2.*sqrt(x.^2+y.^2).*r)/t^2);
pboltz2 = @(x,y) 2*integral(@(r) pboltzo(x,y,r),-R,R);

p2s = [];
r = linspace(0,12,100);
for r=r
    p2s = [p2s, pboltz2(r,0)];
end
r = linspace(0,12,100);
ps = pboltzmann(r);
figure
plot(r,ps,"-r")
hold on
plot(r,p2s,"bo")
xlabel("Distance of ray from origin (m)")
ylabel("Number density")
title("Boltzmann, numerical and anylitical comparison of number density in yz plane")
