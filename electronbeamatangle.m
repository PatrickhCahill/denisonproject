%% Electron Beam LENS
%Input pressure function and symbolic computes refractive index and other
%functions
init()




f = 0.1e11; %approximately 10% of 1 au
R = 10; %Radius of the ideal lens
re = 2.817e-15; %classical electron radius
d = 1; %thickness of ideal lens
lambda = 3e-6;
nc = pi/(2*lambda^2*re);
ne = R^2 *nc/(d*f); %Maximum electron density

w0 = 1/2; %This makes it such that the thickness of the beam is approximately equal to the thickness of an ideal refractive index.

zR = pi*w0^2*1/lambda; %This is the Rayleigh Range. 1 is a place holder for the refractive index. If we add this, then we can get an implicit equation for n. But right now I am interested in an approximate solution.
w = @(y) w0*sqrt(1+(y/zR)^2);


z0 = 0;
r0 = 0;
theta = pi/2;
zprime = @(x,y) (x-z0)*cos(theta) + (y-r0)*sin(theta);
rprime = @(x,y) (x-z0)*sin(theta) + (y-r0)*cos(theta);


p1 = @(x,y) ne*(w0/(w(zprime(x,y))^2))*exp(-2*rprime(x,y)^2/w(zprime(x,y))^2); %Gaussian focussed at (0,0). This is the electron density

z02 = 0;
r02 = -2;
theta2 = -pi/4;
zprime2 = @(x,y) (x-z02)*cos(theta2) + (y-r02)*sin(theta2);
rprime2 = @(x,y) (x-z02)*sin(theta2) + (y-r02)*cos(theta2);


p2 = @(x,y) ne*(w0/(w(zprime2(x,y))^2))*exp(-2*rprime2(x,y)^2/w(zprime2(x,y))^2); %Gaussian focussed at (0,0). This is the electron density

p = @(x,y)  p1(x,y) + p2(x,y);
n = @(x,y) 1+p(x,y)/(2*nc);
gradn= symfun(gradient(n,[x,y]),[x,y]);

ngradn =@(x,y) double( n(x,y).*gradn(x,y)); % represents n(r)*gradient(n(r)) is equal to acceleration function.

% trace() then computes the path of the ray.
%Initial conditions
x0 = -180;       
y0 = 2;       
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