%Written by Kyle Adams, adjusted from modelEQ.m
%This file simulates and plots the model equations with 
% specified percentage increases and decreases in parameter 
% values (50% up and down for us)
function plotsVaryingInfluentialParams
    
    % ###Step 1
    %initialize colors
    L_color = [166/255, 107/255, 97/255]; %brown
    grey_color = [150/255, 150/255, 150/255]; %grey

    %loading parameters
    param = setParameters();

    % ###Step 2
    %solve model with the nominal parameter values
    t0 = 0; tfinal = 30; % simulation time in days
    IC = setInitialConditions(); % get the initial values of the model
    IC = struct2cell(IC); IC = [IC{:}];
    tspan = [t0 tfinal];
    options = odeset('RelTol',1e-12,'AbsTol',1e-12) ;

    %simulates the model / state variables
    [T,Y] = ode45(@(t, y)odefun(t, y, param),...
        [t0 tfinal], IC, options);
    L = Y(:, 1);
    % ###Step 3
    %choose parameters to show their influence on QOI
    [T1, Y1, T2, Y2] = solveWithModifiedParam('aCL', param, tspan); L1 = Y1(:,1); L2 = Y2(:,1);
    [T3, Y3, T4, Y4] = solveWithModifiedParam('dL', param, tspan); L3 = Y3(:,1); L4 = Y4(:,1);
    [T5, Y5, T6, Y6] = solveWithModifiedParam('bCL', param, tspan); L5 = Y5(:,1); L6 = Y6(:,1);
    [T7, Y7, T8, Y8] = solveWithModifiedParam('KC', param, tspan); L7 = Y7(:,1); L8 = Y8(:,1);
    [T9, Y9, T10, Y10] = solveWithModifiedParam('aIC', param, tspan); L9 = Y9(:,1); L10 = Y10(:,1);
    [T11, Y11, T12, Y12] = solveWithModifiedParam('gC', param, tspan); L11 = Y11(:,1); L12 = Y12(:,1);
    
    % ###Step 4
    %plot simulations
    figure;
    hold on;
    tiledlayout(2, 3)

    addTile('aCL', T1, L1, T2, L2, L_color, 30);
    addTile('dL', T3, L3, T4, L4, L_color, 30);
    addTile('bCL', T5, L5, T6, L6, L_color, 30);
    addTile('KC', T7, L7, T8, L8, L_color, 30);
    addTile('aIC', T9, L9, T10, L10, L_color, 30);
    addTile('gC', T11, L11, T12, L12, L_color, 30);
    hold off;

    % ###Step 5
    %plot to show all in favor / all against
    good_params = {'bCL'}; %parameters that if increased, will increase QOI
    bad_params = {'aCL', 'dL', 'KC', 'gC', 'aIC'}; %parameters that if decreased, will decrease QOI
    plotAllForAllAgainst(good_params, bad_params);


    function plotAllForAllAgainst(good_param_set, bad_param_set)
    % This function gets new two new param sets, one where influential parameters were adjusted to be
    % all in favor of QOI (1.5*good, 0.5*bad) and one where all are against QOI
    % (0.5*good, 1.5*bad). non-influential parameters are imported from
    % parameters.m at their nominal values. Then it plots the original
    % solution against the solutions where all influential parameters are
    % in favor of the QOI, and where all influential parameters are against
    % the QOI
        p = setParameters();
        p_favor = p;
        p_against = p;
        for i = 1:length(good_param_set) %this loop replaces each good 
            param_i = good_param_set{i};        %parameter in parameters.m 
            p_favor.(param_i) = 1.5*p.(param_i);  %to be 1.5*nominal in 
            p_against.(param_i) = 0.5*p.(param_i); %p_favor, and 0.5*nominal in p_against
        end
        for i = 1:length(bad_param_set)     %this loop does the reverse^
            param_i = bad_param_set{i};
            p_favor.(param_i) = 0.5*p.(param_i);
            p_against.(param_i) = 1.5*p.(param_i);
        end

        %solve the system with new parameter sets
        [T_favor, Y_favor] = ode45(@(t, y) odefun(t, y, p_favor), tspan, IC, options); L_favor = Y_favor(:,1);
        [T_against, Y_against] = ode45(@(t, y) odefun(t, y, p_against), tspan, IC, options); L_against = Y_against(:,1);
        %plot figure
        figure;
        hold on;
        yline(2e11, 'Color', grey_color, 'Linewidth', 3);
        plot(T, L, 'Color', L_color, 'LineWidth', 3);
        plot(T_favor, L_favor, 'Color', L_color, 'LineWidth', 3, 'LineStyle', "--");
        plot(T_against, L_against, 'Color', L_color, 'LineWidth', 3, 'LineStyle', ":");
        
    
        %Axis settings
        xlim([0 30]);
        xlabel('Time (Days)');
        ylabel('L (cells/\mu L)');
        title('Liver hepatocytes (L)');
    
        %Font size adjustments
        ax = gca;
        ax.Title.FontSize = 10;
        ax.XAxis.FontSize = 10;
        ax.YAxis.FontSize = 10;
    
    %Legend
    legend({("Carrying Capacity of Liver Cells"), ("Nominal Parameter Values"),...
        ("All Influential Parameters in Favor of Liver Cells"),...
        ("All Influential Parameters Against of Liver Cells")}, 'Location', 'southwest', 'FontSize', 16, 'Box', 'off');
    hold off;

    end
 
   function addTile(paramName, x_vals1_5, nominal_param1_5, x_vals_half, nominal_param_half, graph_color, endTime)
    nexttile;
    hold on;
    plot(T, L, 'Color', graph_color, 'LineWidth', 3);
    plot(x_vals1_5, nominal_param1_5, 'Color', graph_color, 'LineWidth', 3, 'LineStyle', "--");
    plot(x_vals_half, nominal_param_half, 'Color', graph_color, 'LineWidth', 3, 'LineStyle', ":");
    
    %Axis settings
    xlim([0 endTime]);
    xticks(0:10:30);
    xlabel('Time (Days)');
    ylabel('L (cells/\mu L)');
    title('Liver hepatocytes (L)');
    
    %Font size adjustments
    ax = gca;
    ax.Title.FontSize = 10;
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;
    
    %Legend
    legend({sprintf('%s', paramName),...
        sprintf('1.5*%s', paramName),...
        sprintf('0.5*%s', paramName)}, 'Location', 'southwest', 'FontSize', 10, 'Box', 'off');
    hold off;
end

function [T1, Y1, T2, Y2] = solveWithModifiedParam(paramName, p, tspan)
    % Solve ODE with modified parameter values (1.5x and 0.5x of nominal param).
    % Inputs:
    %   paramName - String of the parameter to modify
    %   p - Struct containing model parameters
    %   tspan - Time span for ODE solver
    % Outputs:
    %   T1, Y1 - Solution for 1.5 * given parameter
    %   T2, Y2 - Solution for 0.5 * given parameter

    %get initial conditions
    IC = setInitialConditions(); % get the initial values of the model
    IC = struct2cell(IC); IC = [IC{:}];
    %options
    options = odeset('RelTol',1e-12,'AbsTol',1e-12) ;

    %modify parameter and solve with 1.5x value
    p1 = p;
    p1.(paramName) = 1.5 * p.(paramName);
    [T1, Y1] = ode45(@(t, y) odefun(t, y, p1), tspan, IC, options);

    %modify parameter and solve with 0.5x value
    p2 = p;
    p2.(paramName) = 0.5 * p.(paramName);
    [T2, Y2] = ode45(@(t, y) odefun(t, y, p2), tspan, IC);
end


end 
