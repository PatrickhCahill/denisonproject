function output = petertrace(pos,n,x1) %pos = [x0;y0;Tx0;Ty0]
syms x y z

if ~exist('tol','var')
      x1 = pos(1)+2*abs(pos(1));
end
gradx = symfun(diff(n,x),[x,y]);
gradx = @(x,y) double(gradx(x,y));
Tx = @(y) pos(3)+integral(@(x) gradx(x,y),pos(1),x1);
grady = symfun(diff(n,y),[x,y]);
grady = @(x,y) double(grady(x,y));
Ty = @(y) pos(4)+integral(@(x) grady(x,y),pos(1),x1);

xend = x1-pos(2)*Tx(pos(2))/Ty(pos(2));

output = [x1,pos(2);xend,0];
end


