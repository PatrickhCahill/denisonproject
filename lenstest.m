%% Electron Beam LENS convergence
%Input pressure function and symbolic computes refractive index and other
%functions

syms z r
nemax = 1;
nc    = 1;
R     = 100;
z0    = 1;

const = nemax/(2*nc*R^2);

n = @(z,r)  (1+ (R^2-r^2)*const);
gradn= symfun(gradient(n,[z,r]),[z,r]);

deltaT =int(gradn,z,0,z0);

Tfinal = @(z,r) [1;0] + double(deltaT(r));


r0 = 0;       
for i = 0:20:100
    T = Tfinal(z0,r0+i);
    zend = (T(1)/T(2))*(0-r0-i);
    
    plot([0;zend],[r0+i;0],'-.r')
    hold on

end
% name = ("Focal Length = "+f);
% p2 = plot([20,zend],[y0+i,0],'-.r','DisplayName',name);
% legend([p2])




