function QOI = calculateQOI(p)
%This file saves the QOI quantity of interest, values after running the ODE model
    arguments
        p; % parameters
    end
    
    % ###Step 1
    t0 = 0; tfinal = 10; % simulation time in days
    IC = setInitialConditions(p); % get the initial values of the model
    IC = struct2cell(IC); IC = [IC{:}];

    %simulate the model 
    [~,Y] = ode15s(@(t, y)odefun(t, y, p),...
        [t0 tfinal], IC); 

    %we chose the S(t) as our QOI   
    QOI = Y(end,3); % ###Step 2, now takes into account the third equation in ode.fun, might need to run multiple times to get each equation 
    % results
end
