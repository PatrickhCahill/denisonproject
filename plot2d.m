function output = plot2d(pos,ngradn,n,tol)
figure
fcontour(@(x,y) log(n(x,y)),[-200,200,-200,200])%Plots contour map of n(r)

hold on
colorbar

w = trace(pos,ngradn,tol);
plot(w(:,1),w(:,2),'-b');
hold on
plot(w(1,1),w(1,2),'bo','MarkerFaceColor','r');  % Mark the initial position with a red dot

axis equal
xlim([-200 200]);
ylim([-200 200]);

title('Path');
ylabel('y (m)');
xlabel('x (m)');

output = w
end