lambda= 3e-6;
Na = 6.0221409e+23;
eps0 = 8.85418e-12;

M = 1;
alpha = 1.48e-31* (4*pi*eps0);

% mL = 5e-17*(1+linspace(-10,6,100));

% c = (mL).*Na*alpha/(2*pi*M*eps0);
k = 2*pi/(lambda);
mL = 3.6920e-06; %As per Waveaguide3;


a=10;
V = 2.405;
taperconst = sqrt(2)*k*Na*alpha/((2*pi).^2 * M *eps0);

taperangle = taperconst*mL*(1/a);
disp(taperangle) %Very small taper-angle