%% BOLTZMANN LENS
%Input pressure function and symbolic computes refractive index and other
%functions


p = @(a,b) (N/t^3)*(beta0/pi)^(3/2)*exp(-beta0*(a^2+b^2)/t^2); %Boltzmann distribution

n = @(a,b) 1+p(a,b)*2*pi*alpha;
gradn= symfun(gradient(n,[a,b]),[a,b]);

ngradn =@(a,b) double( n(a,b).*gradn(a,b)); % represents n(r)*gradient(n(r)) is equal to acceleration function.

% trace() then computes the path of the ray.
%Initial conditions
x0 = -180;       
y0 = 20;       
Tx0 = 1;      
Ty0 = 0; 
pos = [x0;y0;Tx0;Ty0];

w = plot2d(pos,ngradn,n,1e-12);
figure
plot(w(:,4))
figure
final = extrapolate(w);
plot(final(:,1),final(:,2))


function output = extrapolate(w)
    x1 = w(99,1);
    x2 = w(100,1);
    y1 = w(99,2);
    y2 = w(100,2);
        
    grad = (y2-y1)/(x2-x1);
    xend = (-y2)/grad+x2;
    output = [x2,y2;xend,0];
end