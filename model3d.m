%% 3D BOLTZMANN LENS

%Input pressure function and symbolic computes refractive index and other
%functions


%Boltzmann distribution
p = @(x,y,z)(N/t^3)*(beta0/pi)^(3/2).*exp(-beta0.*(x.^2+y.^2+z.^2)./t^2);



n = @(x,y,z) 1+p(x,y,z)*2*pi*alpha;

%n = @(x,y,z) 2./(1+(x.^2+y.^2+z.^2)/200^2);
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
pos = [x0; y0; z0; Tx0; Ty0; Tz0];


plot3d(pos,ngradn,n,1e-12);