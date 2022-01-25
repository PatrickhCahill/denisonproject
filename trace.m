    %% Plan to clean up trace
% Simplify trace to have optional inputs and save in external file. Inputs include
% 
% 1. Ray  position
% 2. Error tolerance
% 3. ngradn function
% 4. Create trace2d and trace 3d for different setups

%% Function

function output = trace(pos,ngradn,tol,tSpan) %pos = [x0; y0; Tx0; Ty0] or 3d equivalent
if ~exist('tol','var')
      tol = 1e-12;
end
if ~exist('tSpan','var')
    tSpan = linspace(0,500,100);    
end
if length(pos) == 4
    is3d = 0;
else
    is3d = 1;
end



% ode45 then solves this using the derivatives function. @derivatives calls
% the function derivatives as an object. This allows us just to change
% derivatives.
options = odeset('RelTol',tol,'Stats','off');
if is3d ==0
    [t,w] = ode45(@ray_equation2d, tSpan, pos,options); %increase precision relative tolerance
elseif is3d==1   
    [t,w] = ode45(@ray_equation3d, tSpan, pos,options);
end


    % Ray equation for geometrical optics.
    function output = ray_equation2d(tf,wf)
        rf = [wf(1); wf(2)]; %radial vector as function of t
        Tf = [wf(3);wf(4)]; % optical ray vector as function of r.
        
        drdt = Tf;
        dTdt = ngradn(rf(1),rf(2));

        output = [drdt;dTdt];      

    end
    function output = ray_equation3d(tf,wf)
        rf = [wf(1); wf(2); wf(3)]; %radial vector as function of t
        Tf = [wf(4); wf(5); wf(6)]; % optical ray vector as function of r.
        
        drdt = Tf;
        dTdt = ngradn(rf(1),rf(2),rf(3));

        output = [drdt;dTdt];      

    end

    output = w;
end