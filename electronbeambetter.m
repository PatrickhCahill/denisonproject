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
w = @(z) w0*sqrt(1+(z/zR)^2);



p = @(x,y,z) ne*(w0/(w(z)^2))*exp(-2*(x^2+y^2)/w(z)^2); %Gaussian focussed at (0,0). This is the electron density

n = @(x,y,z) 1+p(x,y,z)/(2*nc);
gradn= symfun(gradient(n,[x,y,z]),[x,y,z]);

ngradn =@(x,y,z) double( n(x,y,z).*gradn(x,y,z)); % represents n(r)*gradient(n(r)) is equal to acceleration function.

% trace() then computes the path of the ray.
%Initial conditions
x0 = -180;       
y0 = 2;
z0 = 2;
Tx0 = 1;      
Ty0 = 0;
Tz0 = 0;
pos = [x0;y0;z0;Tx0;Ty0;Tz0];

w = plot2d(pos,ngradn,n,1e-12);
% figure
% plot(w(:,5))
% figure
final = extrapolate(w);
% plot(final(:,1),final(:,2))


function output = extrapolate(w)
    x1 = w(end-1,1);
    x2 = w(end,1);
    y1 = w(end-1,2);
    y2 = w(end,2);
        
    grad = (y2-y1)/(x2-x1);
    xend = (-y2)/grad+x2;
    output1 = [x2,y2;xend,0];
    
    x1 = w(end-1,1);
    x2 = w(end,1);
    y1 = w(end-1,3);
    y2 = w(end,3);
        
    grad = (y2-y1)/(x2-x1);
    xend = (-y2)/grad+x2;
    output2 = [x2,y2;xend,0];
    
    output = [output1,output2];
    
end