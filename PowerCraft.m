init()

global n




% p = @(x,y) (N/t^3)*(beta0/pi)^(3/2)*exp(-beta0*(x.^2+y.^2)/t^2); %Boltzmann distribution
n = @(x,y) 1+p(x,y)*2*pi*alpha;

au = 1.496e11;
xs = (au/100)*(1:10:100);
powers = [];
for x = xs
    powers = [powers;getpower(au)];
end
plot(xs,powers)
ylim([0,1e-11])
xlabel("Radial Position Of Craft")
ylabel("Power On Craft (J)")
xlabel("Radial Position Of Craft (m)")
title("Power on craft vs Radial position")

function p = p(x,y)
modelledpressure()
r =  sqrt(x.^2+y.^2);

[minValue,closestIndex] = min(abs(rs-r));

y1 = ps(closestIndex);
x1 = rs(closestIndex);
y2 = ps(closestIndex+sign(x1-r));
x2 = rs(closestIndex+sign(x1-r));

p = ((y2-y1)/(x2-x1))*(r-y1)+x1;
end

function power = getpower(xpos)
    P0 = 1e9;
    w0 = 10;
    
    r1 = bestr(xpos,-10);
    r2 = bestr(xpos,10);
    power = -P0*exp(-2*r2^2/w0^2)+P0*exp(-2*r1^2/w0^2);
    
end

function bestr =bestr(x,rcraft)
    global n
    rinits = 1:1:80;
    yends = [];
    for rinit = rinits
        yends = [yends; yend(rinit,n,x)];
    end
    [minValue,closestIndex] = min(abs(yends-rcraft));
    
    y1 = yends(closestIndex);
    x1 = rinits(closestIndex);
    y2 = yends(closestIndex+sign(x1-rcraft));
    x2 = rinits(closestIndex+sign(x1-rcraft));
    
    bestr = ((x2-x1)/(y2-y1))*(rcraft-y1)+x1;
end
function yend = yend(r,n,x)
    x0 = -180;       
    Tx0 = 1;      
    Ty0 = 0;
    pos = [x0; r;Tx0;Ty0];
    output = petertrace(pos,n,180,x);
    
    yend = output(2,2);
end