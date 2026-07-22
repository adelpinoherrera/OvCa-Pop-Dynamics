%%adjusted by Kyle Adams, originally Mahya Aghaee

function plotModelSimulation
    % ###Step 1
    % define the colors of the state variables
    S_color  = 'r'; % red
    R_color  = 'g'; % green
    T_color = 'b'; %total population blue
    line_width = 4;

    %load the parameter names and values
    p = setParameters();

    % ###Step 2
    % set simulation timespan and load the initial conditions
    t0 = 0; tfinal = 10; % initial and final simulation times in days
    IC = setInitialConditions(p); % get the initial values of the model
    IC = struct2cell(IC); IC = [IC{:}];

    % ###Step 3
    % choose integration settings using "options"
    options = odeset('RelTol',1e-12,'AbsTol',1e-12);
    
    % simulate the model 
    [Tf,Xf] = ode45(@(t, y)odefun(t, y, p),...
        [t0 tfinal], IC, options);
    
    % ###Step 4
    % store the solutions for each variable for plotting
    SF   = Xf(:,1);
    RF   = Xf(:,2);
    TF   = Xf(:,3); 

    % ###Step 5
    %% Create plots of the variables in a 2x3 grid %%
    % each plot has the same aesthetic and similar labelings
    
    tiledlayout(1,3) %1 row and 2 columns for the plots 

    % Top left plot
    nexttile
    plot(Tf, SF, 'Color', S_color, 'LineWidth',line_width)
    xlim([0 tfinal])
    xlabel('Time (days)')
    ylabel('Sensitive cells')
    title('Sensitive cells over time')
    ax = gca;
    formatAxes(ax);

    % Top middle plot
    nexttile
    plot(Tf,RF,'Color', R_color,'LineWidth', line_width)
    title('Cell Populations')
    xlim([0 tfinal]) 
    xlabel('Time (days)')
    ylabel('Resistant cells')
    title('Resistant cells over time')
    ax = gca;
    formatAxes(ax);

    % Top right plot
    nexttile
    plot(Tf,TF,'Color', T_color,'LineWidth', line_width)
    title('Cell Populations')
    xlim([0 tfinal]) 
    xlabel('Time (days)')
    ylabel('Total cells')
    title('Total cells over time')
    ax = gca;
    formatAxes(ax);


  
    %helper function for formatting axes
    function formatAxes(ax)
        ax.Title.FontSize = 15;
        ax.XAxis.FontSize = 14;
        ax.YAxis.FontSize = 14;
    end

end 
