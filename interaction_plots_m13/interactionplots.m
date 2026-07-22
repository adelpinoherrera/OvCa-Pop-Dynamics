
%get parameters from parameters.m and set seed for reproducibility
random_seed = 1; 
rng(random_seed);

numQOIPoints = 10000;
idx = randperm(numel(QOI), numQOIPoints); % gets numQOIpoints random indices from a total size(QOI) (QOI is saved to MatLab environment after running SobolMain.m)
param_sets = pN(idx); % pN is stored in matlab's environment after SobolMain.m is run
QOI_vals= QOI(idx); % gets the QOIs produced by each set in pN
topNames = {'rS', 'S0','rR', 'R0','K'};  % put in which parameters you want plotted (we used top 6 influential)

paramLabels = {'rS', 'S0','rR', 'R0','K'}; 
%paramColors = {[ 1, 0, 0], [1, 135/255, 0], };

%for 0 drug, 100,000 samples, lb 0.5 and ub 1.5
%topNames = {'rS', 'S0','rR', 'R0','K'}; % put in which parameters you want plotted (we used top 6 influential)
%paramLabels = {'rS', 'S0','rR', 'R0','K'}; 


%for 0.66 drug, 100,000 samples, lb 0.5 and ub 1.5
%topNames = {'rR','rS','R0','alphaS','S0','IC50S','K'};% put in which parameters you want plotted (we used top 6 influential)
%paramLabels = {'rR','rS','R0','alphaS','S0','IC50S','K'};



plotInteractionScatters(param_sets, QOI_vals, topNames, paramLabels);

function plotInteractionScatters(pN, QOI, topNames, paramLabels)%, paramColors)
tiledlayout(2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
%set(gcf, 'Units', 'normalized', 'Position', [0.05, 0.1, 0.9, 0.8]); %new line to make figure wider 

N = numel(pN); % number of model evaluations
numTop = numel(topNames); % number of top influential parameters

topVals = zeros(N, numTop); % initializing matrix to store values of top parameters (values from sensitivity analysis)
for i = 1:numTop
    field = topNames{i};
    topVals(:,i) = cellfun(@(s) s.(field), pN); % for each parameter, we are getting all of its perturbed values from each parameter set
end

combs = nchoosek(1:numTop, 2); % getting all possible index pairs

for i = 1:size(combs, 1)
    idx1 = combs(i, 1); % gets index of first member of parameter pair
    idx2 = combs(i, 2); % gets index of second member of parameter pair

    x = topVals(:, idx1); % gets values of parameter 1
    y = topVals(:, idx2); % gets values of parameter 2
    c = QOI; 

    nexttile;
    scatter(x, y, 10, c, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
    clim([min(QOI) max(QOI)])
    xlabel(paramLabels{idx1}, 'FontWeight', 'bold', 'FontSize', 10);
    xRange = max(x) - min(x);
    xpadding = 0.01 * xRange; % used to get rid of weird whitespace in plots
    xlim([min(x) - xpadding, max(x) + xpadding]); % gives a little extra buffer room
    ylabel(paramLabels{idx2}, 'FontWeight', 'bold', 'FontSize', 10);
    yRange = max(y) - min(y); % used to get rid of weird whitespace in plots
    ypadding = 0.01 * yRange;  % gives a little extra buffer room
    ylim([min(y) - ypadding, max(y) + ypadding]);
    colormap(jet);
    set(gca, 'Box', 'on');
    ax = gca;
ax.FontSize = 8;
end


cb = colorbar(ax, 'eastoutside');
cb.Label.String = 'Total population';
cb.Label.FontWeight = 'bold';
cb.Label.FontSize = 10; %color bar label's font size
cb.Position = [0.96, 0.3, 0.008, 0.25]; % hardcoded based on what the plot looks like
    % [x start, y start, width, height]
end

drawnow; % Forces the figure to finish drawing
exportgraphics(gcf, 'interaction_plot_0.0Ce_sharedK_minandmax_5param.png', 'Resolution', 300);