syms x y z

%Constant and other parameters %% Note that initial ray conditions are in
%trace() function
alpha = 1.642e-30; %polarizability
eps0 = 8.854187817e-12; 
k = 1.38064852e-23; %boltzmann
avagadro =6.02214086e23 ;

t = 0.1;
gasmass = 6000;
temp0 = 90;
molarmass = 39.948;
N = (gasmass/molarmass)*avagadro;
particlemass = molarmass/avagadro/1000;

beta0 = particlemass/(2*k*temp0);