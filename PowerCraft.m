p = @(x,y) (N/t^3)*(beta0/pi)^(3/2)*exp(-beta0*(x.^2+y.^2)/t^2); %Boltzmann distribution
n = @(x,y) 1+p(x,y)*2*pi*alpha;

ys = 1:1:80;
xs = [];
for y0 = ys
    xs = [xs; xend(y0,n)];
end
semilogy(ys,log(xs))
xlabel("Radial distance from centre (m)")
ylabel("Length to intercept with x-axis (m)")
title("Intercept length vs radial distance of ray")


function xend = xend(y,n)
    x0 = -180;       
    Tx0 = 1;      
    Ty0 = 0;
    pos = [x0; y;Tx0;Ty0];
    output = petertrace(pos,n);
    
    xend = output(2,1);
end