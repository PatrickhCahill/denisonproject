%% Electron Beam LENS convergence

type = "gaussianAlongAxis";
graphIntercepts = true;
getpower = true;




syms z r

f = 1.5e11; %approximately 10% of 1 au
R = 10; %Radius of the ideal lens
re = 2.817e-15; %classical electron radius
d = 1; %thickness of ideal lens
lambda = 3e-6;
nc = pi/(2*lambda^2*re);
nemax = R^2 *nc/(d*f); %Maximum electron density

z0    = 1;




switch type
    case "ideal"
        const = nemax/(2*nc*R^2);
        n = @(z,r)  (1+ (R^2-r^2)*const);
    case "gaussianAlongAxis"
        % This is the density function of an electron beam that is Gaussian
        % in nature.
        % There are three free variables w0 the beam waist radius, p0 the
        % maximum density and zR, the Rayleigh Range.
        
        pgaussian = @(z,r,w0,p0,zR) p0*(1/(1+z^2/zR^2))*exp(-(r^2/w0^2)/(1+z^2/zR^2));
        
        w0 = R;
        p0 = nemax;
        angularspread = 0.01*pi;
        zR = (1/2)*angularspread*w0;

        p = @(z,r) pgaussian(z,r,w0,p0,zR);
        n  = @(z,r)  1 + p(z,r)/(2*nc);
    case "multipleGaussian"
        % This is the density function of an electron beam that is Gaussian
        % in nature.
        % There are three free variables w0 the beam waist radius, p0 the
        % maximum density and zR, the Rayleigh Range.
        
        pgaussian = @(z,r,w0,p0,zR) p0*(1/(1+z^2/zR^2))*exp(-(r^2/w0^2)/(1+z^2/zR^2));
        
        w01 = 2*R;
        p01 = nemax/3;
        angularspread1 = 0.01*pi;
        zR1 = (1/2)*angularspread1*w01;
        
        w02 = 2/sqrt(2)*R;
        p02 = 2*nemax/3;
        angularspread2 = 0.01*pi;
        zR2 = (1/2)*angularspread2*w02;

        pbeam1 = @(z,r) pgaussian(z,r,w01,p01,zR1);
        pbeam2 = @(z,r) pgaussian(z,r,w02,p02,zR2);
        p     = @(z,r) pbeam1(z,r) + pbeam2(z,r);
        n  = @(z,r)  1 + p(z,r)/(2*nc);
end





gradn= symfun(gradient(n,[z,r]),[z,r]);
deltaT =int(gradn,z,0,z0);


Tfinal = @(z,r) [1;0] + double(deltaT(r));
r0 = linspace(0,R,20);
zstars = r0*0;
for i = 1:length(r0)
    T = Tfinal(z0,r0(i));
    zend = (T(1)/T(2))*(0-r0(i));
    zstars(i) = zend;
end


if graphIntercepts == true
    figure(1)
    for i = 1:length(r0)-1
            plot([z0;zstars(i)],[r0(i);0],'-.r')
            hold on
    end
    name = ("Graph Of Intercepts");
    p2 = plot([z0,zstars(end)],[r0(end),0],'-.r','DisplayName',name);
    legend([p2])
end

if getpower == true
    rs = 1:R;
    zstarbetter = rs;
    for i =1:length(rs)
       zstarbetter(i) = zstar_inter(rs(i),r0,zstars);
    end
    figure(2)
    plot(rs,zstarbetter/f)
    xlabel("Initial Radius r0 (m)")
    ylabel("Distance to intercept z-axis (au)")
    
    
    
    
end

function zstarapprox = zstar_inter(r,r0,zstars)
    %----------
    %We first search through the list to find the first r0 such that r0(i)
    %is bigger than r. 
    i=1; 
    while r>r0(i) && i<length(r0)
        i = i+ 1;
        
    %If r0(i) == r we return zstars(i)
    if r0(i)==r
        zstarapprox = zstars(i);
        return
    end
    end

    %-------------------------
    %Now we interpolate between our exact values to approximate zstar
    grad = ((zstars(i-1)-zstars(i))/(r0(i-1)-r0(i)));
    zstarapprox = grad*(r-r0(i-1)) + zstars(i-1);
    
end

function findz(craftRadius,craftPosition)
    

end


function power = getpower2(zpos,r)
    P0 = 1e9;
    w0 = 10;
    
    r1 = bestr(xpos,-10);
    r2 = bestr(xpos,10);
    power = -P0*exp(-2*r2^2/w0^2)+P0*exp(-2*r1^2/w0^2);
    
end