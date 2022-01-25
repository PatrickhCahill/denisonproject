%% FISHEYE LENS
%Input pressure function and symbolic computes refractive index and other
%functions

n = @(x,y) 2./(1+(a.^2+b.^2)); %Newton's fisheye lens refractive index
gradn= symfun(gradient(n,[a,b]),[a,b]);

ngradn =@(x,y) double( n(x,y).*gradn(x,y)); % represents n(r)*gradient(n(r)) is equal to acceleration function.


%Plots contour map of n(r)
fcontour(@(x,y) log(n(x,y)),[-3,3,-3,3])
hold on
colorbar


% trace() then computes the path of the ray.
%Initial conditions
x0 = -1/sqrt(2);       
y0 = 1/sqrt(2);       
Tx0 = 1;      
Ty0 = 0; 
pos = [x0;y0;Tx0;Ty0];
tSpan = linspace(0,25,100);    

w = trace(pos,ngradn,1e-12,tSpan);
plot(w(:,1),w(:,2))