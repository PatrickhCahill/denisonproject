% coulombsConstant = 1;
% electronMass = 1;
% 
% 
% electronPositions = [0,0,0,0;1,1,-1,-1;1,-1,1,-1];
% electronSpeeds = [1,1,1,1;0,0,0,0;0,0,0,0];
% dt = 0.1;
% t=0;
% while t<1
% %Update electronPositions
% electronPositions = electronPositions+electronSpeeds*dt;
% 
% %Loop to find forces
% for index = 1:length(electronPositions)
%     electron = electronPositions(:,index);
%     electronForce = [0;0;0];
%     %electronForce is given by the sum of all the component forces
%     for otherIndex = 1:length(electronPositions)
%         electron2 = electronPositions(:,otherIndex);
%        if index ~= otherIndex
%            rvector = electron-electron2;
%            forceComponent = coulombsConstant.*rvector./norm(rvector).^3;
%            electronForce = electronForce + forceComponent;
%        end
%     end
%     %Update speeds using forces
%     electronSpeeds(:,index) = electronSpeeds(:,index) + electronForce.*dt./electronMass;
% end
% t = t+dt;
% end
phi = (1+sqrt(5))/2;
n=1:100;
thetas = (2*pi/phi^2).*n;
rs = sqrt(n);
plot(rs.*cos(thetas),rs.*sin(thetas),"o")
   


