%% BOLTZMANN LENS MULTIRAY
%Input pressure function and symbolic computes refractive index and other
%functions

init()
p = @(x,y) (N/t^3)*(beta0/pi)^(3/2).*exp(-beta0*(x.^2+y.^2)/t^2); %Boltzmann distribution

n = @(x,y) 1+p(x,y)*2*pi*alpha;
gradn= symfun(gradient(n,[x,y]),[x,y]);

ngradn =@(x,y) double( n(x,y).*gradn(x,y)); % represents n(r)*gradient(n(r)) is equal to acceleration function.



% trace() then computes the path of the ray.
%Initial conditions
x0 = -180;       
y0 = 5;       
Tx0 = 1;      
Ty0 = 0; 
for i = 0:5:45
    pos = [x0;y0+i;Tx0;Ty0];
    w = trace(pos,ngradn);
    final = extrapolate(w);
    plot(final(:,1),final(:,2),'-b')
    hold on
end
xlim([0,20e9])
xlabel("x (m)")
ylabel("y (m)")

function output = extrapolate(w)
    x1 = w(99,1);
    x2 = w(100,1);
    y1 = w(99,2);
    y2 = w(100,2);
        
    grad = (y2-y1)/(x2-x1);
    xend = (-y2)/grad+x2;
    output = [x2,y2;xend,0];
end
