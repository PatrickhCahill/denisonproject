function plot3d(pos,ngradn,n,tol)
figure
w = trace(pos,ngradn,tol);
ns = log(n(w(:,1)',w(:,2)',w(:,3)')-1);
surf([w(:,1)';w(:,1)'],[w(:,2)';w(:,2)'],[w(:,3)';w(:,3)'],[ns;ns],'facecolor','none','edgecolor','interp','linewidth',2);
cb = colorbar;
cb.Label.String = '{log (n(r)-1)}';
axis equal
xlim([-200 200]);
ylim([-200 200]);
zlim([-200 200]);

title('Path');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
hold on
for i = -2:2
    for j = -2:2
        if i ~= 0 || j ~=0
            x0 = -180;       
            y0 = 2+i*50;
            z0 = 2+j*50;
            Tx0 = 1;      
            Ty0 = 0; 
            Tz0 = 0;
            pos = [x0; y0; z0; Tx0; Ty0; Tz0];
            w = trace(pos,ngradn,tol);
            ns = log(n(w(:,1)',w(:,2)',w(:,3)')-1);
            surf([w(:,1)';w(:,1)'],[w(:,2)';w(:,2)'],[w(:,3)';w(:,3)'],[ns;ns],'facecolor','none','edgecolor','interp','linewidth',2);
        end
    end
end


end