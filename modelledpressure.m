%%Free Expansion
init()

time = t;
global n
n = 1000; %Number of layers approximating


r0 = 1;
rs = ones(n,1);%Creates initial radiuses.
for k = 1:n
   rs(k) = ((k/n).^(1/3)).*r0; 
end
ds = zeros(n,1);


tSpan = [0 time];
pos = [rs;ds];
options = odeset('RelTol',1e-12,'Stats','off');

[t,w] = ode45(@gas0, tSpan, pos,options);
rs = w(length(t),1:n);
ps = pfinals(rs);



%% dydt is the differential equation that is computed.
function dydt = gas0(t,y)
    %Here are the relevant constants and initial conditions
    global n
    R = 8.314; %Ideal gas constant
    gamma = 5/3;
    kb = 1.38064852e-23;
    
    M = 6000;
    MM = 39.948;
    T0 = 90; %Kelvin
    r0 = 1;
    V0 = (4/3)*pi*r0^3;
    const = (4*pi)*(n^(1-gamma)/MM)*V0^(gamma-1)*R*T0;
    
    m = M/n;
    delta =1000;
    w = 1000;
    
    %Wrangling the data
    rs = y(1:end/2); %Gets radii
    ss = y(end/2+1:end); % and velocities
    
    %The following defines a function that computes the (V)^-gamma, where V
    %is a shell around a sphere
    Vkgminus = @(k) 1/((4*pi/3)*(rs(k).^3 - rs(k-1).^3)).^gamma;
    
    %%%%% Adding mean free path evaporation to outter layer
    
    Tn = T0*((4/(3*n))*pi*r0^3)^(gamma-1)/((4*pi/3)*(rs(end)^3-rs(end-1)^3))^(gamma-1);
    
    p0 = (3/(4*pi*r0^3))*(M/MM)*R*T0;
    meanfreepath = kb*T0/(sqrt(2)*pi*((2*rs(end))^2))*((4*pi/3)*(rs(end)^3-rs(end-1)^3))/(p0*(4/(3*n))*pi*r0^3);
    rs(end) = rs(end)-meanfreepath;
    %%%%%
    
    %The first drs term is the is the inner most shell.
    drs = [ss(1)];
    dss = [const.*(rs(1).^2).* (((4*pi/3)*(rs(1).^3)).^(-gamma)-Vkgminus(2))];
    
    %This loop then calculates the next dt step for the layers
    for k = 2:(n-1)
        drs = [drs; ss(k)];
        dss = [dss;const.*(rs(k).^2).* (Vkgminus(k)-Vkgminus(k+1))];%Relevant ODE
    end
    
    %This computes the outter most layer
    drs = [drs; ss(n)];
      dss = [dss;const.*(rs(n).^2).* (Vkgminus(n))]; %Old model
%     dss = [dss;4*pi*rs(n)^2*R*Tn/(MM*((4*pi/3)*(rs(end)^3-rs(end-1)^3)))]; %mean free path model
%     dss = [dss; const.*(rs(n).^2).*Vkgminus(n)-(4*pi/m)*(rs(n).^2)*(p0 + delta*sin(w*t(end)))]; %Wave test

    dydt = [drs; dss];
end
function p = ptk(rs,t,k)
    global n
    R = 8.314;
    gamma = 5/3;
    MM = 39.948;
    T0 = 90;
    M = 6000;
    
    r0 = 1;
    V0 = (4/3)*pi*r0^3;
    
    C = n^(-gamma)* V0^(gamma-1)*M/MM*R*T0;
    
    Tk = @(k) T0*((4/(3*n))*pi*r0^3)^(gamma-1)/((4*pi/3)*(rs(k)^3-rs(k-1)^3))^(gamma-1);
    T1 = T0*((4/(3*n))*pi*r0^3)^(gamma-1)/((4*pi/3)*(rs(1)^3))^(gamma-1);
    
    Vkgminus = @(t,k) ((4*pi/3).*(rs(t,k).^3 - rs(t,k-1).^3)).^(-gamma);
%     ps = [C.*((4*pi/3)*(rs(1).^3)).^(-gamma).*(6.02214086e23)./(R*T1)];
%     for k = 2:length(rs)
%         ps = [ps; C.*Vkgminus(k).*(6.02214086e23)./(R*Tk(k))];
%     end
    if(k==1)
        p = C.*((4*pi/3).*(rs(t,k).^3)).^(-gamma);
    else
        p = C.*Vkgminus(t,k);
    end
end



function ps = pfinals(rs)
    global n
    R = 8.314;
    gamma = 5/3;
    MM = 39.948;
    T0 = 90;
    M = 6000;
    
    r0 = 1;
    V0 = (4/3)*pi*r0^3;
    
    C = n^(-gamma)* V0^(gamma-1)*M/MM*R*T0;
    
    Tk = @(k) T0*((4/(3*n))*pi*r0^3)^(gamma-1)/((4*pi/3)*(rs(k)^3-rs(k-1)^3))^(gamma-1);
    T1 = T0*((4/(3*n))*pi*r0^3)^(gamma-1)/((4*pi/3)*(rs(1)^3))^(gamma-1);
    
    Vkgminus = @(k) ((4*pi/3)*(rs(k).^3 - rs(k-1).^3)).^(-gamma);
    ps = [C.*((4*pi/3)*(rs(1).^3)).^(-gamma).*(6.02214086e23)./(R*T1)];
    for k = 2:length(rs)
        ps = [ps; C.*Vkgminus(k).*(6.02214086e23)./(R*Tk(k))];
    end
%     ps = [C.*((4*pi/3)*(rs(1).^3)).^(-gamma)];
%     for k = 2:length(rs)
%         ps = [ps; C.*Vkgminus(k)];
%     end
%     ps = ps.*(6.02214086e23)./(R*T0);
end
