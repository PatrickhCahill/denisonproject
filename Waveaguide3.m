lambda= 3e-6;
Na = 6.0221409e+23;
eps0 = 8.85418e-12;

M = 1;
alpha = 1.48e-31* (4*pi*eps0);

% mL = 5e-17*(1+linspace(-10,6,100));

% c = (mL).*Na*alpha/(2*pi*M*eps0);
k = 2*pi/(lambda);


a=10;
V = 2.405;

mL = a*(-k*a+sqrt(k^2 * a^2 + V^2))/k * (2*pi*M*eps0/(Na*alpha));

disp(mL)
% a = c.*k*1./sqrt(V.^2-2.*c.*k^2);
% 
% plot(mL,a)
