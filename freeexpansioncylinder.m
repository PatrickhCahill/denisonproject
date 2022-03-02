%%Free Expansion
global n
n = 200; %Number of layers approximating


r0 = 1;
rs = ones(n,1);%Creates initial radiuses.
for k = 1:n
   rs(k) = ((k/n).^(1/2)).*r0; 
end
ds = zeros(n,1);


tSpan = [0 1];
pos = [rs;ds];


options = odeset('RelTol',1e-12,'Stats','off');
[t,w] = ode45(@gas0, tSpan, pos);

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

vidObj = VideoWriter("freeexpansioncylinder");
open(vidObj)
axis tight
set(gca,'nextplot','replacechildren')


for time = 1:5:length(t)
rfinals = w(time,:)';
rfinals = rfinals(1:end/2);
ps = pfinals(rfinals);

plot(rfinals,ps,"-b")
hold on

xlim([0 25])
%ylim([0 10^26])
xlabel("Radius (m)")
ylabel("Pressure ({\times p_0})")
title("Pressure vs Radius at t="+round(t(time)*1000,0)+"ms. n="+n)
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
    R = 8.314; %Ideal gas constant
    gamma = 5/3;
    kb = 1.38064852e-23;
    
    M = 6000;
    MM = 39.948;
    T0 = 90; %Kelvin
    r0 = 1;
    V0 = (4/3)*pi*r0^3;
    const = 2*R*T0*(r0.^2/n).^(gamma-1)/MM;
    
    %Wrangling the data
    rs = y(1:end/2); %Gets radii
    ss = y(end/2+1:end); % and velocities

    
    rhandled = @(k) (rs(k).^2 - rs(k-1).^2).^(-gamma);
    
    %%%%% Adding mean free path evaporation to outter layer
    
%     Tn = T0*((4/(3*n))*pi*r0^3)^(gamma-1)/((4*pi/3)*(rs(end)^3-rs(end-1)^3))^(gamma-1);
%     
%     p0 = (3/(4*pi*r0^3))*(M/MM)*R*T0;
%     meanfreepath = kb*T0/(sqrt(2)*pi*((2*rs(end))^2))*((4*pi/3)*(rs(end)^3-rs(end-1)^3))/(p0*(4/(3*n))*pi*r0^3);
%     rs(end) = rs(end)-meanfreepath;
    %%%%%
    
    %The first drs term is the is the inner most shell.
    drs = [ss(1)];
    dss = [const.*(rs(1)).*(rs(1).^(-2*gamma)-rhandled(2))];
    
    %This loop then calculates the next dt step for the layers
    for k = 2:(n-1)
        drs = [drs; ss(k)];
        dss = [dss;const.*(rs(k)).*(rhandled(k)-rhandled(k+1))];%Relevant ODE
    end
    
    %This computes the outter most layer
    drs = [drs; ss(n)];
    dss = [dss;const.*(rs(end)).*rhandled(n)]; %Old model
    %dss = [dss;4*pi*rs(n)^2*R*Tn/(MM*((4*pi/3)*(rs(end)^3-rs(end-1)^3)))];
    dydt = [drs; dss];
end

function ps = pfinals(rs)
    global n
    r0 = 1;
    gamma = 5/3;

    rhandled = @(k) (rs(k).^2 - rs(k-1).^2).^(-gamma);
    
    ps = [(1/n)^(gamma)*r0^(2*gamma)/rs(1)^(2*gamma)];
    for k = 2:(n)
        ps = [ps;(1/n)^(gamma)*r0^(2*gamma)*rhandled(k)];%Relevant ODE
    end
    
end

