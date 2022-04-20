%% Electron Beam LENS convergence
%Input pressure function and symbolic computes refractive index and other
%functions
init()

f = 1.5e11/1e6; %approximately 10% of 1 au
R = 10; %Radius of the ideal lens
re = 2.817e-15; %classical electron radius
d = 1; %thickness of ideal lens
lambda = 3e-6;
nc = pi/(2*lambda^2*re);
ne = R^2 *nc/(d*f)*1000; %Maximum electron density

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
p3 = @(x,y) ne*(w0/(w(y)^2))*exp(-2*x^2/w(y)^2); %Gaussian focussed at (0,0). This is the electron density

p = @(x,y)  p1(x,y) + p2(x,y)+p3(x,y);
n = @(x,y) 1+p(x,y)/(2*nc);
gradn= symfun(gradient(n,[x,y]),[x,y]);
ngradn =@(x,y) double( n(x,y).*gradn(x,y));


gradx = symfun(diff(n,x),[x,y]);
gradx = @(x,y) double(gradx(x,y));
Tx = @(y) integral(@(x) gradx(x,y),-180,180);
grady = symfun(diff(n,y),[x,y]);
grady = @(x,y) double(grady(x,y));
Ty = @(y) integral(@(x) grady(x,y),-180,180);


x0 = 0;       
y0 = 0;       
for i = 0:1:5
    xend = x0-(y0+i)*(Tx(y0+i)+1)/(Ty(y0+i));
    
    plot([20;xend],[y0+i;0],'-.r')
    hold on

end
name = ("Focal Length = "+f);
p2 = plot([20;xend],[y0+i;0],'-.r','DisplayName',name);
legend([p2])






