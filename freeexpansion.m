tspan = [0 0.1];
y0 = [1;1];
[t,y] = ode45(@myplot, tspan, y0);


plot(t,1./(y(:,1).^3),'-o')


function output = myplot(tf,wf)
    mm = 39.948;
    R = 8.314;
    T = 90;
    k = 3*mm*R*T/(16*pi^2);

    output = [wf(2);k/(wf(1).^5)];      
end