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
ne = R^2 *nc/(d*f); %Maximum electron density



p = @(x,y,z) ne*exp(-2*(x^2+y^2)); %Gaussian focussed at (0,0). This is the electron density

n = @(x,y,z) 1+p(x,y,z)/(2*nc);
gradn= symfun(gradient(n,[x,y,z]),[x,y,z]);

deltaT =int(gradn,x,-180,180);

Tfinal = @(y,z) [1;0;0] + double(deltaT(y,z));

x0 = 180;

rayline = @(x,y0,z0) [x0;y0;z0] + Tfinal(y0,z0)*(x-x0);
r = @(x,y0,z0) sqrt(subsref(rayline(x,y0,z0), struct('type', '()', 'subs', {{2}})).^2 + subsref(rayline(x,y0,z0), struct('type', '()', 'subs', {{3}})).^2);

bestx = @(y0,z0) fminbnd(@(x) abs(r(x,y0,z0)),180,10000);
best = @(y0,z0) r(bestx(y0,z0),y0,z0);

y= linspace(-5,5,10);
z= linspace(-5,5,10);
[Y,Z] = meshgrid(y,z);
finalx = 1e4;
for i = 1:length(y)
    for j = 1:length(z)
        %Y(i,j) = best(y(i),z(j));
        Y(i,j) = r(finalx,y(i),z(j));
    end
end

imagesc([min(y),max(y)],[min(z),max(z)],Y)
colorbar
xlabel("y0 (m)")
ylabel("z0 (m)")
title("Distance to x-axis at x="+finalx+"m")

% x0 = 0;       
% y0 = 0;       
% for i = 0:1:1
%     xend = x0-(y0+i)*(Tx(y0+i)+1)/(Ty(y0+i));
%     
%     plot([20;xend],[y0+i;0],'-.r')
%     hold on
% 
% end
% name = ("Focal Length = "+f);
% p2 = plot([20;xend],[y0+i;0],'-.r','DisplayName',name);
% legend([p2])






