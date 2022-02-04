%Zeroth case
tspan = [0 100];
y0 = [1;0.9;0.8;0.7;0;0;0;0];
[t,y] = ode45(@gas0,tspan,y0);

plot(t,y(:,1),'-b')
hold on
plot(t,y(:,2),'-r')
hold on
plot(t,y(:,3),'-g')
hold on
plot(t,y(:,4),'-k')
title('Solution of Zeroth Case');
xlabel('Time t');
ylabel('Solution y');

function dydt = gas0(t,y)
    kn = 1;
    p0 = 0;
    gamma = 5/3;
    %y is a 2n length vector
    rs = y(1:end/2);
    ss = y(end/2+1:end);
    drs = [ss(1)];
    dss = [kn*rs(1)^(2-3*gamma)-p0*rs(1)^2];
    for n = 2:length(rs)
        drs = [drs; ss(n)];
        dss = [dss;kn*rs(n)^(2-3*gamma)-rs(n-1)^(-3*gamma)*rs(n)^2];
    end
    dydt = [drs; dss];
end
