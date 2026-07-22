% This file runs the Sobol sensitivity analysis on a QSP ODE model
% This code originally stems from code by Jaimit Parikh, and was modified
% by Kyle Adams for this project.

%FOR OVARIAN
%Need to run every iteration of the sensitivity with Ce=1 and Ce=0, can be
%changed in odefun.m

%get parameters from parameters.m and set seed for reproducibility
random_seed = 1; %this is the same seed as compareQOIDists, which was purposeful, but not necessary
rng(random_seed);
p = setParameters();
%ce = 0.0; %define the treatment condition here too, not in the setParameters function 


% ###Step 1
%make output folder with a timestamp
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
outdir = fullfile('sensitivityRuns', ['run_1000samples_totalQOI_1Ce_75_25_m13_0.8_1.2_', timestamp]);
mkdir(outdir);

%save particular parameter values used to generate sensitivity indices
save(fullfile(outdir, 'params_used.mat'), '-struct', 'p');

% ###Step 2
%%Set the sampling bounds for the parameters
%%Our lower bounds are 50% of nominal values, and upper bounds are 150%
%can change to 0.5 and 1.5 
lower_percentage = 0.8;
upper_percentage = 1.2;
base_samples = 1000; %100 samples or 100,000 samples ? 
param_dist = {'Uniform'};

%% ---------- %%
% 1. Set parameter names and length
paramNames = fieldnames(p);
numParam = length(paramNames);
%initialize vectors of length numParam with 0s for parameter bounds
lowBounds = zeros(1, numParam); 
upBounds = zeros(1, numParam);

% 2. Set bounds for parameters
%populates lowBounds and upBounds with bounds listed in parameters.m
%allow the initial condition for S and R to vary from 0 to 100
for i = 1:numParam 
    name = paramNames{i};
    switch name
        case 'S0'
            lowBounds(i) = 1;     % S0 from 0 to 100
            upBounds(i)  = 100;
        case 'R0'
            lowBounds(i) = 1;     % R0 from 0 to 100
            upBounds(i)  = 100;
        otherwise
            lowBounds(i) = p.(name) * lower_percentage;
            upBounds(i)  = p.(name) * upper_percentage;
    end
end

%store this info in structure called parsObj; used in getSamplesSobol
parsObj.name = paramNames'; %transpose so the dimensions of each field match
parsObj.lb = num2cell(repmat(-inf, 1, length(parsObj.name)));
parsObj.ub = num2cell(inf(1, length(parsObj.name)));
parsObj.dist = repmat(param_dist, 1, numParam);
parsObj.N = base_samples; %number of desired Sobol base samples
%below is a field that stores each parameter's bounds
parsObj.parameters = arrayfun(@(i) {lowBounds(i), upBounds(i)}, 1:numParam, 'UniformOutput', false); %changed this line
samples = generateSobolSamples(parsObj, false);

%%
% 3. Simulate model for desired parameter set and extract QOI values
parsName = parsObj.name;
QOI = zeros(1, length(samples)); %initialize vector of QOI values with 0's
pN = cell(1,length(samples));
parfor ii = 1:length(samples)
pN{ii} = updatePars(p, parsName, samples(ii, :));
QOI(ii) = calculateQOI(pN{ii});
end
%store QOI distribution for histogram comparison later
QOI_all_varying = QOI;

%% 
% 4. Make scatter plot of QOI vs selected parameters
% we used this for checking the sensitivity analysis results
plotScatter(samples, QOI, parsName, ...
'scatterSobol.png');
ylabel('QOI')

%%
% 5. Estimate Sobol Index values
Y = QOI';
S = calculateSobolIndices(Y, length(parsName), parsObj.N, 1, 0.95);
Si = calculateSobolIndices(Y, length(parsObj.name), parsObj.N);
mytable = cell2table([Si.S1, Si.ST], 'VariableNames', ["S1", "ST"],...
'RowNames',parsObj.name);
sortTable = sortrows(mytable, 1, 'descend');
disp(sortTable)

%prepare for figure in order of descending total sensitivity indices
S1 = cell2mat(S.S1);
ST = cell2mat(S.ST);
[ST_sorted, idx] = sort(ST, 'descend');
S1_sorted = S1(idx); %idx are indices in order of descending ST value
paramNames_sorted = paramNames(idx);

% ###Step 3
%plot sensitivity indices
hold on;
figure('DefaultAxesFontSize', 16);
b = bar([S1_sorted, ST_sorted], 'grouped');
b(1).FaceColor = [1 0.79 0.63]; %light orange - bar color for S1 values
b(2).FaceColor = [0 0.188 0.69]; %dark blue - bar color for ST values
set(gca, 'FontName', 'Times New Roman')
%maximize figure on screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
xticks(1:length(paramNames_sorted));
xticklabels(paramNames_sorted)
xtickangle(45);
ylabel('Sensitivity Index', 'FontSize', 16, 'FontName', 'serif');
legend({'S1', 'ST'}, 'FontSize', 16, 'FontName', 'serif');
hold off;

%create scatterplot of a parameter's values vs QOI values
function plotScatter(samples, QOI, parsName, fname)
f = figure('DefaultAxesFontSize', 14);
tiledlayout('flow')
for ii = 1:size(samples, 2)
nexttile;
plot(samples(:, ii), QOI, 'ko'); xlabel(parsName(ii));
ylabel('QOI');
end
exportgraphics(f, fname, 'resolution', 300);
end
function p = updatePars(p, parsName, parsValue)
for ii = 1:length(parsName)
p.(parsName{ii}) = parsValue(ii);
end
end

%plots top 14 parameters on a bar graph with S1 and ST side by side
plot_top_params(14, S1_sorted, ST_sorted, paramNames_sorted)


%save sensitivity indices
save(fullfile(outdir, ['sensitivity_results_s1' timestamp '.mat']), 'S1');
save(fullfile(outdir, ['sensitivity_results_sT' timestamp '.mat']), 'ST');
T_ST = table(paramNames(:), ST(:), 'VariableNames', {'Parameter', 'ST'});
writetable(T_ST, fullfile(outdir, ['sensitivity_resultsST_' timestamp '.csv']));

T_S1 = table(paramNames(:), S1(:), 'VariableNames', {'Parameter', 'S1'});
writetable(T_S1, fullfile(outdir, ['sensitivity_resultsS1_' timestamp '.csv']));

%save all 13 files that were used to generate the results, with a time stamp
code_files = {'sobolMain.m', 'setParameters.m', 'odefun.m', 'calculateQOI.m', 'compareQOIDistributions.m', 'setInitialConditions.m',...
    'transformSobolSamplesDistributions.m', 'generateSobolSamples.m', 'makeSensitivityTables.m', 'plotModelSimulation.m', 'plotsVaryingInfluentialParams.m', ...
    'calculateSobolIndices.m', 'generateRandomSubset.m'};
for k = 1:length(code_files)
[~, name, ext] = fileparts(code_files{k});
dest_filename = [name '_' timestamp ext];
copyfile(code_files{k}, fullfile(outdir, dest_filename));
end

% ###Step 4
%optional: plot only top influential parameters
%plot_top_params(16, S1_sorted, ST_sorted, paramNames_sorted);

function plot_top_params(num_params_plotted, S1_sorted, ST_sorted, paramNames_sorted)
trunc_S1_sorted = S1_sorted(1:num_params_plotted, :);
trunc_ST_sorted = ST_sorted(1:num_params_plotted, :);
hold on;
figure('DefaultAxesFontSize', 16);
b = bar([trunc_S1_sorted, trunc_ST_sorted], 'grouped');
b(1).FaceColor = [1 0.79 0.63];
b(2).FaceColor = [0 0.188 0.69];
set(gca, 'FontName', 'Times New Roman')
%maximize figure on screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]); 
xticks(1:num_params_plotted);
xticklabels(paramNames_sorted(1:num_params_plotted))
xtickangle(45);
ylabel('Sensitivity Index', 'FontSize', 16, 'FontName', 'serif');
legend({'S1', 'ST'}, 'FontSize', 16, 'FontName', 'serif');
hold off;
end
