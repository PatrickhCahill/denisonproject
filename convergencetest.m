%% FISHEYE LENS CONVERGENCE TEST
%Input pressure function and symbolic computes refractive index and other
%functions
init()
n = @(x,y) 2/(1+(x.^2+y.^2)); %Newton's fisheye lens refractive index
gradn= symfun(gradient(n,[x,y]),[x,y]);

ngradn =@(x,y) double( n(x,y).*gradn(x,y)); % represents n(r)*gradient(n(r)) is equal to acceleration function.

% trace() then computes the path of the ray.
%Initial conditions
x0 = -1;       
y0 = 0;       
Tx0 = 0;      
Ty0 = 1; 
pos = [x0; y0;Tx0;Ty0];
tSpan = linspace(0,25,100);
errorlist = [];
for i = 1:1:12
    tol = 10^(-1*i);
    w = trace(pos,ngradn,tol,tSpan);
    error = abs(w(end,1).^2+w(end,2).^2-1);
    errorlist = [errorlist;tol,error];
end

figure
loglog(errorlist(:,1),errorlist(:,2))
set(gca, 'xdir','reverse')
xlabel("Relative Error Tolerance")
ylabel("Distance From Convergence")
title("Convergence Test")
