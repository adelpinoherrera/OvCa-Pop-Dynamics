% This code was originally written by Jaimit Parikh, and was edited by 
% Kyle Adams to compare QOI distributions when some parameters are fixed.

random_seed = 1;
rng(random_seed);  %Setting seed 

% ###Step 1
%set up bounds and distributions for Sobol sampling
lower_percentage = 0.5;
upper_percentage = 1.5;
base_samples = 175000;
param_dist = {'Uniform'};
p = setParameters();

% ###Step 2
%define parameter sets to compare QOI distributions
param_sets = {
    'All', fieldnames(p)', ...
    'Top 6 Influential Parameters Varied', {'aCL', 'dL', 'bCL', 'KC', 'gC', 'aIC'}, ...
    'Least 29 Influential Parameters Varied', {'lL', 'dA', 'aAH', 'bAH', 'aRA', 'bRA', ...
                                            'aIRA', 'bIRA', 'gH', 'KH', 'aIH', 'bIH', 'dH', ...
                                            'aHC', 'bHC', 'bIC', 'dC', 'sR', 'dR', 'aIR', 'bIR', ...
                                            'aCI', 'bCI', 'aHI', 'bHI', 'lC', 'lH', 'lR', 'dI'} ...
    };

num_param_sets = length(param_sets) / 2;
QOIs = cell(1, num_param_sets);
labels = cell(1, num_param_sets);

%perform sensitivity analysis for each parameter set
for i = 1:num_param_sets
    set_name = param_sets{2*i - 1};
    varying_params = param_sets{2*i};
    labels{i} = set_name;

    %using the parameter names from param_sets, extract all
    %parameter info from parameters.m
    p_subset = struct();
    for j = 1:length(varying_params)
        field = varying_params{j};
        if isfield(p, field)
            p_subset.(field) = p.(field);
        end
    end

    param_names = fieldnames(p_subset);
    num_params = length(param_names);
    lowBounds = zeros(1, num_params);
    upBounds = zeros(1, num_params);

    for j = 1:num_params
        param = param_names{j};
        lowBounds(j) = p_subset.(param) * lower_percentage;
        upBounds(j) = p_subset.(param) * upper_percentage;
    end

    parsObj.name = param_names';
    parsObj.lb = num2cell(repmat(-inf, 1, length(parsObj.name)));
    parsObj.ub = num2cell(inf(1, length(parsObj.name)));
    parsObj.dist = repmat(param_dist, 1, num_params);
    parsObj.N = base_samples;
    parsObj.parameters = arrayfun(@(i) {'lower', lowBounds(i), 'upper', upBounds(i)}, ...
                                  1:num_params, 'UniformOutput', false);

    %get varied parameter value sets by running getSamplesSobol
    samples = generateSobolSamples(parsObj, false);

    %simulate the QOI
    QOI = zeros(1, length(samples));
    pN = cell(1,length(samples));
    parsName = parsObj.name;
    parfor ii = 1:length(samples)
        pN{ii} = updatePars(p, parsName, samples(ii, :));
        QOI(ii) = calculateQOI(pN{ii});
    end
    QOIs{i} = QOI;
end
%take equal-sized subsets of bigger QOI distributions to compare fairly
smallest_num_model_evals = min(cellfun(@length, QOIs));

% ###Step 3
%plot histogram comparisons of QOI distributions
figure;
hold on;
colors = {'red', 'black', 'blue'}; % all params, top params, bottom params
for i = 1:num_param_sets
    histogram(generateRandomSubset(QOIs{i}, smallest_num_model_evals, random_seed), ...
        'EdgeColor', colors{i}, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
end
legend(labels, 'FontSize', 12, 'FontName', 'serif', 'Location', 'northwest');
xlabel('QOI Value', 'FontSize', 12, 'FontName', 'serif');
ylabel('Frequency', 'FontSize', 12, 'FontName', 'serif');
title('Histogram Comparison of QOI Distributions', 'FontSize', 14, 'FontName', 'serif');
xlim([0 2e11]);
hold off;

%update parameters function
function p = updatePars(p, parsName, parsValue)
    for ii = 1:length(parsName)
        p.(parsName{ii}) = parsValue(ii);
    end
end
