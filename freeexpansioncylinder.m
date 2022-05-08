%%Free Expansion
global n
n = 5; %Number of layers approximating


r0 = 1;
rs = ones(n,1);%Creates initial radiuses.
for k = 1:n
   rs(k) = ((k/n).^(1/2)).*r0; 
end
ds = zeros(n,1);


tSpan = [0 1e-8];
pos = [rs;ds];


options = odeset('RelTol',1e-15,'Stats','off');
[t,w] = ode45(@gas0, tSpan, pos,options);

% figure
% for time = 1:length(t)
% rfinals = w(time,:)';
% rfinals = rfinals(1:end/2);
% ps = pfinals(rfinals);
% semilogy(rfinals,ps,"-b")
% 
% 
% xlim([0 25])
% ylim([0 10^0])
% xlabel("Radius (m)")
% ylabel("Pressure ({\times p_0})")
% title("Pressure vs Radius at t="+round(t(time)*1000,0)+"ms. n="+n)
% hold off
% pause(0.05)
% end
% figure 
% r1s = w(:,1);
% plot(t,pr(r1s))

vidObj = VideoWriter("electronBeam");
open(vidObj)
axis tight
set(gca,'nextplot','replacechildren')

ytickslist = [1];
for ytick = 1:18
    ytickslist = [yticks, 10^ytick];
end


for time =1:5:length(t)
rfinals = w(time,:)';
rfinals = rfinals(1:end/2);
disp(rfinals)
ps = pfinals(rfinals);

loglog(rfinals,ps,"-b")
hold on

xlim([0 25])
ylim([0 10^18])


xlabel("Radius (m)")
ylabel("Electron Number Density (m^{-3})")
title("Electron Density vs Radius at t="+round(t(time)*1e12,0)+"ps. n="+n)
legend("Onion-Model")
hold off
currframe = getframe(gcf);
writeVideo(vidObj,currframe);


end
close(vidObj)
%% dydt is the differential equation that is computed.
function dydt = gas0(t,y)
    %Here are the relevant constants and initial conditions
    global n
    q = -1.602e-19; %Coulombs
    re= 2.82e-15; %m Classical electron radius
    me = 9.1e-31; %kg Electron Mass
    eps0 = 8.854187817e-12;
    mu = (4*pi)*10^(-7);
    c = 3e8; %m/s speed of light
    hbar = 1e-34;
    
    lambda = 3e-6; %m Light wavelength
    R= 10; %m Radius of lens
    f = 1.5e10; %m 10% of an au.
    d = 1; %Approximate length thickness. Won't actually matter what this is set to.
    ne0 = pi/(re*lambda.^2)*(R.^2/(f*d)); %Initial density of electrons m^-3.
    vd = 0.2*c;
    
    r00 = 1;
    r0 = ones(n,1);%Creates initial radiuses.
    for k = 1:n
       r0(k) = ((k/n).^(1/2)).*r00; 
    end
    
    %
    D =((q.^2) * ne0 * (r0(1).^2)/(2*me))*(1/(eps0)-mu*vd.^2);
    C =((3*pi.^2).^(2/3)*hbar.^2/(5*me))*(ne0*r0(1).^2).^(5/3)*(2/me)/(ne0*r0(1).^2);
    %Wrangling the data
    rs = y(1:end/2); %Gets radii
    ss = y(end/2+1:end); % and velocities

    
    rhandled = @(k) (rs(k).^(2) - rs(k-1).^2).^(-5/3);
    
    %%%%% Adding mean free path evaporation to outter layer
    
%     Tn = T0*((4/(3*n))*pi*r0^3)^(gamma-1)/((4*pi/3)*(rs(end)^3-rs(end-1)^3))^(gamma-1);
%     
%     p0 = (3/(4*pi*r0^3))*(M/MM)*R*T0;
%     meanfreepath = kb*T0/(sqrt(2)*pi*((2*rs(end))^2))*((4*pi/3)*(rs(end)^3-rs(end-1)^3))/(p0*(4/(3*n))*pi*r0^3);
%     rs(end) = rs(end)-meanfreepath;
    %%%%%
    
    %The first drs term is the is the inner most shell.
    drs = [ss(1)];
    dss = [D./rs(1) + C.*rs(1).*rs(1)^(-10/3)];
    %This loop then calculates the next dt step for the layers
    for k = 2:(n-1)
        drs = [drs; ss(k)];
        dss = [dss;D*rs(k)/(rs(k).^2-rs(k-1).^2) +  C.*rs(k).*(rhandled(k)-rhandled(k+1))];%Relevant ODE
    end
    
    %This computes the outter most layer
    drs = [drs; ss(n)];
    dss = [dss;2*D*rs(k)/(rs(k).^2-rs(k-1).^2) + C.*rs(k).*rhandled(k)]; %Old model
    dydt = [drs; dss];
end

function ps = pfinals(rs)
    global n
    q = -1.602e-19; %Coulombs
    re= 2.82e-15; %m Classical electron radius
    me = 9.1e-31; %kg Electron Mass
    eps0 = 8.854187817e-12;
    mu = (4*pi)*10^(-7);
    c = 3e8; %m/s speed of light
    hbar = 1e-34;
    
    lambda = 3e-6; %m Light wavelength
    R= 10; %m Radius of lens
    f = 1.5e10; %m 10% of an au.
    d = 1; %Approximate length thickness. Won't actually matter what this is set to.
    ne0 = pi/(re*lambda.^2)*(R.^2/(f*d)); %Initial density of electrons m^-3.
    vd = 0.2*c;
    
    r00 = 1;
    r0 = ones(n,1);%Creates initial radiuses.
    for k = 1:n
       r0(k) = ((k/n).^(1/2)).*r00; 
    end
    
    
    ps = [ne0*r0(1).^2/(rs(1).^2)];
    for k = 2:(n)
        ps = [ps;ne0*r0(1).^2/(rs(k).^2-rs(k-1).^2)];%Relevant ODE
    end
    
end

