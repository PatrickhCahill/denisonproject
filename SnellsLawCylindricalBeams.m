ray1 = [-0.7,0];


R=1;
n1 = 1;
n2 = 0.7;
r = n1/n2;

ray2 = [ray1(1),sqrt(R^2-ray1(1)^2)];



theta = (atan(ray2(2)/ray2(1)));

theta2 = asin(r*sin(theta));

m = -tan(theta-theta2);

b = R*sin(theta) + R*cos(theta)*tan(theta-theta2);

ys = roots([(1+m^2),(2*b),(b^2-R^2*m^2)]);
if (ys(1) ~= ys(2))
    if abs(ys(1)-ray1(2))<abs(ys(2)-ray1(2))
        y = ys(2);
    else
        y = ys(1);
    end
else
    y = real(ys(1));
end
x = sqrt(R^2 - y^2);

theta3 = atan(abs((grad-y/x)/(1+grad*y/x)));
theta4 = asin((1/r)*sin(theta3));


plot([-1,ray2(1)],[ray2(2),ray2(2)],"-b");
hold on

plot([ray2(1),x],[ray2(2),y],"-r");
fplot(@(t) R*cos(t),@(t) R*sin(t))
