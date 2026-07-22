%Try different treatment models for A2780 for 14 days, the other cell lines
%are in other scripts 
%of data 
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Starts with A2780 treated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Import all data
A2780_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC50_1um_30k_cellsinwell.csv');
A2780_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC75_0.62um_30k_cellsinwell.csv');
A2780_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC25_1.47um_30k_cellsinwell.csv');

A2780_t_days = A2780_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
A2780_IC50_y_raw = [A2780_IC50_data.Rep1, A2780_IC50_data.Rep2, A2780_IC50_data.Rep3];
A2780_IC25_y_raw = [A2780_IC25_data.Rep1, A2780_IC25_data.Rep2, A2780_IC25_data.Rep3];
A2780_IC75_y_raw = [A2780_IC75_data.Rep1, A2780_IC75_data.Rep2, A2780_IC75_data.Rep3];

A2780_10t_days = A2780_t_days(1:11);
A2780_IC25_10y_raw = A2780_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
A2780_IC50_10y_raw = A2780_IC50_y_raw(1:11,:);
A2780_IC75_10y_raw = A2780_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 1.00; 
ce_IC25 = 0.62;
ce_IC75 = 1.47;

%Define logistic growth constants 
rS = 0.8547; %growth rate (day^-1)
KS = 1114404.9924;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 0.9960177; %the IC50 for sensitive cells, after 72hr exposure

%MODEL 1 (Equation 3)
%%%Estimate for n and alpha - A2780 (constant D) multiplicative drug kill term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Estimate for parameters
p0 = [2; 0.5];      % [n; alpha]
lb = [0.01; 0.01];   % reasonable bounds
ub = [10; 2.0];

res_fun_nandalpha_constD = @(p) all_datasets_residual_reps_n_and_alpha_constD_2params(p, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

%run globalsearch with 200 starts
problem_nandalpha_constD = createOptimProblem(...
    'lsqnonlin', ... %can't use lsqcurvefit because it expects x and y arguments
    'x0', p0, ...
    'objective', res_fun_nandalpha_constD, ...
    'lb', lb, ...
    'ub', ub);

ms_nandalpha_constD = MultiStart('UseParallel', true, 'Display', 'iter', 'StartPointsToRun','bounds-ineqs');  % Set true if desired
nStarts = 200; %more start to ensure reproducibility 
[d_hat, resnorm_d] = run(ms_nandalpha_constD, problem_nandalpha_constD, nStarts);
fprintf('Fitted n and alpha - mult - constD A2780 (10 days)\n')
fprintf('Fitted Hill coefficient n = %.4f\n', d_hat(1));
fprintf('Fitted alpha term a = %.4f\n', d_hat(2));
fprintf('Sum of squared residuals = %.2f\n', resnorm_d);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted n and alpha - mult - constD A2780 (10 days)
% Fitted Hill coefficient n = 3.0355
% Fitted alpha term a = 0.4771
% Sum of squared residuals = 2642273852290.18

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_2params(d_hat, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

fprintf('Model (n=%.2f, a=%.2f):\n', d_hat(1), d_hat(2));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (n=3.04, a=0.48):
%   BIC IC50  = 789.39
%   BIC IC25  = 786.37
%   BIC IC75  = 814.36
%   BIC total = 2390.12

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {A2780_IC50_10y_raw, A2780_IC25_10y_raw, A2780_IC75_10y_raw};

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(A2780_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstantD(d_hat(1), d_hat(2), A2780_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(A2780_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 1400000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (n=%.2f, a=%.2f):\n', d_hat(1), d_hat(2)));

saveas(gcf, fullfile(outdir,'A2780_treated_nandalpha_mult_constD_10days.png'));

writematrix(S_raw_all{1,1}, fullfile(outdir, 'A2780_IC50treated_nandalpha_mult_constD_10days.csv'));
writematrix(S_raw_all{1,2}, fullfile(outdir, 'A2780_IC25treated_nandalpha_mult_constD_10days.csv'));
writematrix(S_raw_all{1,3}, fullfile(outdir, 'A2780_IC75treated_nandalpha_mult_constD_10days.csv'));
 
%MODEL 2 (Equation 4)
%%%Estimate for n and alpha - A2780 (constant D) additive drug kill term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated';

%Import all data
A2780_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC50_1um_30k_cellsinwell.csv');
A2780_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC75_0.62um_30k_cellsinwell.csv');
A2780_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC25_1.47um_30k_cellsinwell.csv');

A2780_t_days = A2780_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
A2780_IC50_y_raw = [A2780_IC50_data.Rep1, A2780_IC50_data.Rep2, A2780_IC50_data.Rep3];
A2780_IC25_y_raw = [A2780_IC25_data.Rep1, A2780_IC25_data.Rep2, A2780_IC25_data.Rep3];
A2780_IC75_y_raw = [A2780_IC75_data.Rep1, A2780_IC75_data.Rep2, A2780_IC75_data.Rep3];

%downsample data
A2780_10t_days = A2780_t_days(1:11);
A2780_IC25_10y_raw = A2780_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
A2780_IC50_10y_raw = A2780_IC50_y_raw(1:11,:);
A2780_IC75_10y_raw = A2780_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 1.00; 
ce_IC25 = 0.62;
ce_IC75 = 1.47;

%Define logistic growth constants 
rS = 0.8547; %growth rate (day^-1)
KS = 1114404.9924;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 0.9960177; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [2; 0.5];      % [n; alpha]
lb = [0.01; 0.01];   % reasonable bounds
ub = [10; 2.0];

res_fun_nandalpha_constD_add = @(p) all_datasets_residual_reps_n_and_alpha_constD_add(p, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

%run globalsearch with 200 starts
problem_nandalpha_constD_add = createOptimProblem(...
    'lsqnonlin', ... %can't use lsqcurvefit because it expects x and y arguments
    'x0', p0, ...
    'objective', res_fun_nandalpha_constD_add, ...
    'lb', lb, ...
    'ub', ub);

ms_nandalpha_constD_add = MultiStart('UseParallel', true, 'Display', 'iter', 'StartPointsToRun','bounds-ineqs');  % Set true if desired
nStarts = 200; %more start to ensure reproducibility 
[a_hat, resnorm_a] = run(ms_nandalpha_constD_add, problem_nandalpha_constD_add, nStarts);
fprintf('Fitted n and alpha - add - constD A2780 (10 days)\n')
fprintf('Fitted Hill coefficient n = %.4f\n', a_hat(1));
fprintf('Fitted alpha term a = %.4f\n', a_hat(2));
fprintf('Sum of squared residuals = %.2f\n', resnorm_a);
% 123 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted n and alpha - add - constD A2780 (10 days)
% Fitted Hill coefficient n = 3.8278
% Fitted alpha term a = 0.6237
% Sum of squared residuals = 1869851695003.98

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_add(a_hat, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

fprintf('Model (n=%.2f, a=%.2f):\n', a_hat(1), a_hat(2));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (n=3.83, a=0.62):
%   BIC IC50  = 761.67
%   BIC IC25  = 767.72
%   BIC IC75  = 810.77
%   BIC total = 2340.16

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {A2780_IC50_10y_raw, A2780_IC25_10y_raw, A2780_IC75_10y_raw};

model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(A2780_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstantD_add(a_hat(1),a_hat(2), A2780_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    model_fit_all{i} = S_fit_mean;
    plot(A2780_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 1400000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (n=%.2f, a=%.2f):\n', a_hat(1), a_hat(2)));

saveas(gcf, fullfile(outdir,'A2780_treated_nandalpha_add_constD_10days.png'));

writematrix([A2780_10t_days, model_fit_all{1}], fullfile(outdir,'A2780_IC50treated_nandalpha_add_constD_10days.csv'));
writematrix([A2780_10t_days, model_fit_all{2}], fullfile(outdir,'A2780_IC25treated_nandalpha_add_constD_10days.csv'));
writematrix([A2780_10t_days, model_fit_all{3}], fullfile(outdir,'A2780_IC75treated_nandalpha_add_constD_10days.csv'));

writematrix(S_raw_all{1,1}, fullfile(outdir, 'A2780_IC50treated_nandalpha_add_constD_10days.csv'));
writematrix(S_raw_all{1,2}, fullfile(outdir, 'A2780_IC25treated_nandalpha_add_constD_10days.csv'));
writematrix(S_raw_all{1,3}, fullfile(outdir, 'A2780_IC75treated_nandalpha_add_constD_10days.csv'));

%MODEL 3 (Equation 5)
%%%Estimate for n, alpha and d - A2780 (constant D) multiplicative drug
%%% term and additive death term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated';

%Import all data
A2780_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC50_1um_30k_cellsinwell.csv');
A2780_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC75_0.62um_30k_cellsinwell.csv');
A2780_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC25_1.47um_30k_cellsinwell.csv');

A2780_t_days = A2780_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
A2780_IC50_y_raw = [A2780_IC50_data.Rep1, A2780_IC50_data.Rep2, A2780_IC50_data.Rep3];
A2780_IC25_y_raw = [A2780_IC25_data.Rep1, A2780_IC25_data.Rep2, A2780_IC25_data.Rep3];
A2780_IC75_y_raw = [A2780_IC75_data.Rep1, A2780_IC75_data.Rep2, A2780_IC75_data.Rep3];

%downsample data
A2780_10t_days = A2780_t_days(1:11);
A2780_IC25_10y_raw = A2780_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
A2780_IC50_10y_raw = A2780_IC50_y_raw(1:11,:);
A2780_IC75_10y_raw = A2780_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 1.00; 
ce_IC25 = 0.62;
ce_IC75 = 1.47;

%Define logistic growth constants 
rS = 0.8547; %growth rate (day^-1)
KS = 1114404.9924;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 0.9960177; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [2; 0.5; 0.5];      % [n; alpha]
lb = [0.01; 0.01; 0.01];   % reasonable bounds
ub = [20; 2.0; 4.0];

res_fun_nandalpha_constD_death = @(p) all_datasets_residual_reps_n_and_alpha_constD_death(p, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

%run globalsearch with 200 starts
problem_nandalpha_constD_death = createOptimProblem(...
    'lsqnonlin', ... %can't use lsqcurvefit because it expects x and y arguments
    'x0', p0, ...
    'objective', res_fun_nandalpha_constD_death, ...
    'lb', lb, ...
    'ub', ub);

ms_nandalpha_constD_death = MultiStart('UseParallel', true, 'Display', 'iter', 'StartPointsToRun','bounds-ineqs');  % Set true if desired
nStarts = 200; %more start to ensure reproducibility 
[mdeath_hat, resnorm_mdeath] = run(ms_nandalpha_constD_death, problem_nandalpha_constD_death, nStarts);
fprintf('Fitted n, alpha and d - mult drug - add death - constD A2780 (10 days)\n')
fprintf('Fitted Hill coefficient n = %.4f\n', mdeath_hat(1));
fprintf('Fitted alpha term a = %.4f\n', mdeath_hat(2));
fprintf('Fitted death term d = %.4f\n', mdeath_hat(3));
fprintf('Sum of squared residuals = %.2f\n', resnorm_mdeath);
% 199 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted n, alpha and d - mult drug - add death - constD A2780 (10 days)
% Fitted Hill coefficient n = 20.0000
% Fitted alpha term a = 0.2819
% Fitted death term d = 0.1146
% Sum of squared residuals = 2224373561342.70

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_death(mdeath_hat, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

fprintf('Model (n=%.2f, a=%.2f, d=%.2f):\n', mdeath_hat(1), mdeath_hat(2), mdeath_hat(3));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (n=20.00, a=0.28, d=0.11):
%   BIC IC50  = 783.78
%   BIC IC25  = 774.29
%   BIC IC75  = 816.99
%   BIC total = 2375.07


%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {A2780_IC50_10y_raw, A2780_IC25_10y_raw, A2780_IC75_10y_raw};

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(A2780_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstantD_death(mdeath_hat(1), mdeath_hat(2), mdeath_hat(3), A2780_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(A2780_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 1400000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (n=%.2f, a=%.2f, d=%.2f):\n', mdeath_hat(1), mdeath_hat(2), mdeath_hat(3)));

saveas(gcf, fullfile(outdir,'A2780_treated_n-alphaandd_multdrug_adddeath_constD_10days.png'));

writematrix(S_raw_all{1,1}, fullfile(outdir, 'A2780_IC50treated_n-alphaandd_multdrug_adddeath_constD_10days.csv'));
writematrix(S_raw_all{1,2}, fullfile(outdir, 'A2780_IC25treated_n-alphaandd_multdrug_adddeath_constD_10days.csv'));
writematrix(S_raw_all{1,3}, fullfile(outdir, 'A2780_IC75treated_n-alphaandd_multdrug_adddeath_constD_10days.csv'));

%MODEL 4 (Equation 6)
%%%Estimate for n, alpha and d - A2780 (constant D), additive drug kill term
%%%and death term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated';

%Import all data
A2780_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC50_1um_30k_cellsinwell.csv');
A2780_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC75_0.62um_30k_cellsinwell.csv');
A2780_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC25_1.47um_30k_cellsinwell.csv');

A2780_t_days = A2780_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
A2780_IC50_y_raw = [A2780_IC50_data.Rep1, A2780_IC50_data.Rep2, A2780_IC50_data.Rep3];
A2780_IC25_y_raw = [A2780_IC25_data.Rep1, A2780_IC25_data.Rep2, A2780_IC25_data.Rep3];
A2780_IC75_y_raw = [A2780_IC75_data.Rep1, A2780_IC75_data.Rep2, A2780_IC75_data.Rep3];

%downsample data
A2780_10t_days = A2780_t_days(1:11);
A2780_IC25_10y_raw = A2780_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
A2780_IC50_10y_raw = A2780_IC50_y_raw(1:11,:);
A2780_IC75_10y_raw = A2780_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 1.00; 
ce_IC25 = 0.62;
ce_IC75 = 1.47;

%Define logistic growth constants 
rS = 0.8547; %growth rate (day^-1)
KS = 1114404.9924;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 0.9960177; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [2; 0.5; 0.5];      % [n; alpha, d]
lb = [0.01; 0.01; 0.01];   % reasonable bounds
ub = [20; 2.0; 4.0];

res_fun_nandalpha_constD_add_death = @(p) all_datasets_residual_reps_n_and_alpha_constD_add_death(p, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

%run globalsearch with 200 starts
problem_nandalpha_constD_add_death = createOptimProblem(...
    'lsqnonlin', ... %can't use lsqcurvefit because it expects x and y arguments
    'x0', p0, ...
    'objective', res_fun_nandalpha_constD_add_death, ...
    'lb', lb, ...
    'ub', ub);

ms_nandalpha_constD_add_death = MultiStart('UseParallel', true, 'Display', 'iter', 'StartPointsToRun','bounds-ineqs');  % Set true if desired
nStarts = 200; %more start to ensure reproducibility 
[adeath_hat, resnorm_adeath] = run(ms_nandalpha_constD_add_death, problem_nandalpha_constD_add_death, nStarts);
fprintf('Fitted n, alpha and d - add drug - add death - constD A2780 (10 days)\n')
fprintf('Fitted Hill coefficient n = %.4f\n', adeath_hat(1));
fprintf('Fitted alpha term a = %.4f\n', adeath_hat(2));
fprintf('Fitted death term d = %.4f\n', adeath_hat(3));
fprintf('Sum of squared residuals = %.2f\n', resnorm_adeath);
% 194 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted n, alpha and d - add drug - add death - constD A2780 (10 days)
% Fitted Hill coefficient n = 7.3169
% Fitted alpha term a = 0.4634
% Fitted death term d = 0.0752
% Sum of squared residuals = 1868663606022.08

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_add_death(adeath_hat, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

fprintf('Model (n=%.2f, a=%.2f, d=%.2f):\n', adeath_hat(1), adeath_hat(2), adeath_hat(3));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (n=7.32, a=0.46, d=0.08):
%   BIC IC50  = 765.09
%   BIC IC25  = 771.20
%   BIC IC75  = 814.25
%   BIC total = 2350.55

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {A2780_IC50_10y_raw, A2780_IC25_10y_raw, A2780_IC75_10y_raw};

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(A2780_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstantD_add_death(adeath_hat(1),adeath_hat(2),adeath_hat(3), A2780_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(A2780_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 1400000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (n=%.2f, a=%.2f, d=%.2f):\n', adeath_hat(1), adeath_hat(2), adeath_hat(3)));

saveas(gcf, fullfile(outdir,'A2780_treated_n-alphaandd_adddrug_adddeath_constD_10days.png'));

writematrix(S_raw_all{1,1}, fullfile(outdir, 'A2780_IC50treated_n-alphaandd_adddrug_adddeath_constD_10days.csv'));
writematrix(S_raw_all{1,2}, fullfile(outdir, 'A2780_IC25treated_n-alphaandd_adddrug_adddeath_constD_10days.csv'));
writematrix(S_raw_all{1,3}, fullfile(outdir, 'A2780_IC75treated_n-alphaandd_adddrug_adddeath_constD_10days.csv'));

%MODEL 5 (Equation 15)
%%%Estimate for n, alpha and beta - A2780 (constant D) multiplicative hill term,
%%%alpha is time dependent 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated';

%Import all data
A2780_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC50_1um_30k_cellsinwell.csv');
A2780_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC75_0.62um_30k_cellsinwell.csv');
A2780_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC25_1.47um_30k_cellsinwell.csv');

A2780_t_days = A2780_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
A2780_IC50_y_raw = [A2780_IC50_data.Rep1, A2780_IC50_data.Rep2, A2780_IC50_data.Rep3];
A2780_IC25_y_raw = [A2780_IC25_data.Rep1, A2780_IC25_data.Rep2, A2780_IC25_data.Rep3];
A2780_IC75_y_raw = [A2780_IC75_data.Rep1, A2780_IC75_data.Rep2, A2780_IC75_data.Rep3];

%downsample data
A2780_10t_days = A2780_t_days(1:11);
A2780_IC25_10y_raw = A2780_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
A2780_IC50_10y_raw = A2780_IC50_y_raw(1:11,:);
A2780_IC75_10y_raw = A2780_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 1.00; 
ce_IC25 = 0.62;
ce_IC75 = 1.47;

%Define logistic growth constants 
rS = 0.8547; %growth rate (day^-1)
KS = 1114404.9924;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 0.9960177; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [2; 0.5; 5];      % [n; alpha; beta]
lb = [0.01; 0.01; 0.01];   % reasonable bounds
ub = [20; 20.0; 20.0];

res_fun_nandalpha_constD_mult_alpha_t = @(p) all_datasets_residual_reps_n_and_alpha_constD_mult_alpha_t(p, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

%run globalsearch with 200 starts
problem_nandalpha_constD_mult_alpha_t = createOptimProblem(...
    'lsqnonlin', ... %can't use lsqcurvefit because it expects x and y arguments
    'x0', p0, ...
    'objective', res_fun_nandalpha_constD_mult_alpha_t, ...
    'lb', lb, ...
    'ub', ub);

ms_nandalpha_constD_multi_alpha_t = MultiStart('UseParallel', true, 'Display', 'iter', 'StartPointsToRun','bounds-ineqs');  % Set true if desired
nStarts = 200; %more start to ensure reproducibility 
[malphat_hat, resnorm_malphat] = run(ms_nandalpha_constD_multi_alpha_t, problem_nandalpha_constD_mult_alpha_t, nStarts);
fprintf('Fitted n, alpha and beta - mult drug - alpha T - constD A2780 (10 days)\n')
fprintf('Fitted Hill coefficient n = %.4f\n', malphat_hat(1));
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(2));
fprintf('Fitted beta term b = %.4f\n', malphat_hat(3));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% 75 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted n, alpha and beta - mult drug - alpha T - constD A2780 (10 days)
% Fitted Hill coefficient n = 2.3378
% Fitted alpha term a = 0.0136
% Fitted beta term b = 8.5286
% Sum of squared residuals = 870083330718.38

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_mult_alpha_t(malphat_hat, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

fprintf('Model (n=%.2f, a=%.2f, b=%.2f):\n', malphat_hat(1), malphat_hat(2), malphat_hat(3));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (n=2.34, a=0.01, b=8.53):
%   BIC IC50  = 720.58
%   BIC IC25  = 774.40
%   BIC IC75  = 778.61
%   BIC total = 2273.58

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {A2780_IC50_10y_raw, A2780_IC25_10y_raw, A2780_IC75_10y_raw};

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(A2780_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_mult_alpha_t(malphat_hat(1), malphat_hat(2), malphat_hat(3), A2780_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(A2780_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 1400000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (n=%.2f, a=%.2f, b=%.2f):\n', malphat_hat(1), malphat_hat(2), malphat_hat(3)));

saveas(gcf, fullfile(outdir,'A2780_treated_n-alphaandbeta_multdrug_alphaT_constD_10days.png'));

writematrix(S_raw_all{1,1}, fullfile(outdir, 'A2780_IC50treated_n-alphaandbeta_multdrug_alphaT_constD_10days.csv'));
writematrix(S_raw_all{1,2}, fullfile(outdir, 'A2780_IC25treated_n-alphaandbeta_multdrug_alphaT_constD_10days.csv'));
writematrix(S_raw_all{1,3}, fullfile(outdir, 'A2780_IC75treated_n-alphaandbeta_multdrug_alphaT_constD_10days.csv'));

%MODEL 6 (Equation 16)
%%%Estimate for n, alpha and beta - A2780 (constant D) additive hill term,
%%%alpha is time dependent 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated';

%Import all data
A2780_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC50_1um_30k_cellsinwell.csv');
A2780_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC75_0.62um_30k_cellsinwell.csv');
A2780_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC25_1.47um_30k_cellsinwell.csv');

A2780_t_days = A2780_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
A2780_IC50_y_raw = [A2780_IC50_data.Rep1, A2780_IC50_data.Rep2, A2780_IC50_data.Rep3];
A2780_IC25_y_raw = [A2780_IC25_data.Rep1, A2780_IC25_data.Rep2, A2780_IC25_data.Rep3];
A2780_IC75_y_raw = [A2780_IC75_data.Rep1, A2780_IC75_data.Rep2, A2780_IC75_data.Rep3];

%downsample data
A2780_10t_days = A2780_t_days(1:11);
A2780_IC25_10y_raw = A2780_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
A2780_IC50_10y_raw = A2780_IC50_y_raw(1:11,:);
A2780_IC75_10y_raw = A2780_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 1.00; 
ce_IC25 = 0.62;
ce_IC75 = 1.47;

%Define logistic growth constants 
rS = 0.8547; %growth rate (day^-1)
KS = 1114404.9924;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 0.9960177; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [2; 0.5; 0.5];      % [n; alpha; beta]
lb = [0.01; 0.01; 0.01];   % reasonable bounds
ub = [20; 20.0; 20.0]; %change ub for the last one 

res_fun_nandalpha_constD_add_alpha_t = @(p) all_datasets_residual_reps_n_and_alpha_constD_add_alpha_t(p, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

%run globalsearch with 200 starts
problem_nandalpha_constD_add_alpha_t = createOptimProblem(...
    'lsqnonlin', ... %can't use lsqcurvefit because it expects x and y arguments
    'x0', p0, ...
    'objective', res_fun_nandalpha_constD_add_alpha_t, ...
    'lb', lb, ...
    'ub', ub);

ms_nandalpha_constD_add_alpha_t = MultiStart('UseParallel', true, 'Display', 'iter', 'StartPointsToRun','bounds-ineqs');  % Set true if desired
nStarts = 200; %more start to ensure reproducibility 
[aalphat_hat, resnorm_aalphat] = run(ms_nandalpha_constD_add_alpha_t, problem_nandalpha_constD_add_alpha_t, nStarts);
fprintf('Fitted n, alpha and beta - add drug - alpha T - constD A2780 (10 days)\n')
fprintf('Fitted Hill coefficient n = %.4f\n', aalphat_hat(1));
fprintf('Fitted alpha term a = %.4f\n', aalphat_hat(2));
fprintf('Fitted beta term b = %.4f\n', aalphat_hat(3));
fprintf('Sum of squared residuals = %.2f\n', resnorm_aalphat);
% 73 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted n, alpha and beta - add drug - alpha T - constD A2780 (10 days)
% Fitted Hill coefficient n = 3.7704
% Fitted alpha term a = 0.0100
% Fitted beta term b = 12.6975
% Sum of squared residuals = 694781817274.98

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_add_alpha_t(aalphat_hat, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

fprintf('Model (n=%.2f, a=%.2f, b=%.2f):\n', aalphat_hat(1), aalphat_hat(2), aalphat_hat(3));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (n=3.77, a=0.01, b=12.70):
%   BIC IC50  = 733.51
%   BIC IC25  = 756.95
%   BIC IC75  = 773.83
%   BIC total = 2264.29

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {A2780_IC50_10y_raw, A2780_IC25_10y_raw, A2780_IC75_10y_raw};

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(A2780_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_add_alpha_t(aalphat_hat(1), aalphat_hat(2), aalphat_hat(3), A2780_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(A2780_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 1400000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (n=%.2f, a=%.2f, b=%.2f):\n', aalphat_hat(1), aalphat_hat(2), aalphat_hat(3)));

saveas(gcf, fullfile(outdir,'A2780_treated_n-alphaandbeta_adddrug_alphaT_constD_10days.png'));

writematrix(S_raw_all{1,1}, fullfile(outdir, 'A2780_IC50treated_n-alphaandbeta_adddrug_alphaT_constD_10days.csv'));
writematrix(S_raw_all{1,2}, fullfile(outdir, 'A2780_IC25treated_n-alphaandbeta_adddrug_alphaT_constD_10days.csv'));
writematrix(S_raw_all{1,3}, fullfile(outdir, 'A2780_IC75treated_n-alphaandbeta_adddrug_alphaT_constD_10days.csv'));

%MODEL 7 (Equation 7)
%%%Estimate for n and alpha - A2780 (constant D) multiplicative hill term,
%%%alpha is time dependent, drug kill is 0 at time 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated';

%Import all data
A2780_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC50_1um_30k_cellsinwell.csv');
A2780_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC75_0.62um_30k_cellsinwell.csv');
A2780_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC25_1.47um_30k_cellsinwell.csv');

A2780_t_days = A2780_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
A2780_IC50_y_raw = [A2780_IC50_data.Rep1, A2780_IC50_data.Rep2, A2780_IC50_data.Rep3];
A2780_IC25_y_raw = [A2780_IC25_data.Rep1, A2780_IC25_data.Rep2, A2780_IC25_data.Rep3];
A2780_IC75_y_raw = [A2780_IC75_data.Rep1, A2780_IC75_data.Rep2, A2780_IC75_data.Rep3];

%downsample data
A2780_10t_days = A2780_t_days(1:11);
A2780_IC25_10y_raw = A2780_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
A2780_IC50_10y_raw = A2780_IC50_y_raw(1:11,:);
A2780_IC75_10y_raw = A2780_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 1.00; 
ce_IC25 = 0.62;
ce_IC75 = 1.47;

%Define logistic growth constants 
rS = 0.8547; %growth rate (day^-1)
KS = 1114404.9924;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 0.9960177; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [2; 0.5];      % [n; alpha]
lb = [0.01; 0.01];   % reasonable bounds
ub = [20; 20];

res_fun_nandalpha_constD_mult_alpha_t0 = @(p) all_datasets_residual_reps_n_and_alpha_constD_mult_alpha_t0(p, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

%run globalsearch with 200 starts
problem_nandalpha_constD_mult_alpha_t0 = createOptimProblem(...
    'lsqnonlin', ... %can't use lsqcurvefit because it expects x and y arguments
    'x0', p0, ...
    'objective', res_fun_nandalpha_constD_mult_alpha_t0, ...
    'lb', lb, ...
    'ub', ub);

ms_nandalpha_constD_multi_alpha_t0 = MultiStart('UseParallel', true, 'Display', 'iter', 'StartPointsToRun','bounds-ineqs');  % Set true if desired
nStarts = 200; %more start to ensure reproducibility 
[malphat_hat, resnorm_malphat] = run(ms_nandalpha_constD_multi_alpha_t0, problem_nandalpha_constD_mult_alpha_t0, nStarts);
fprintf('Fitted n, and alpha - mult drug - alpha T t0=0 - constD A2780 (14 days)\n')
fprintf('Fitted Hill coefficient n = %.4f\n', malphat_hat(1));
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(2));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% 189 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted n, alpha and beta - mult drug - alpha T t0=0 - constD A2780 (14 days)
% Fitted Hill coefficient n = 2.0622
% Fitted alpha term a = 0.1201
% Sum of squared residuals = 824823053693.37

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_mult_alpha_t0(malphat_hat, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

fprintf('Model (n=%.2f, a=%.2f):\n', malphat_hat(1), malphat_hat(2));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (n=2.06, a=0.12):
%   BIC IC50  = 715.18
%   BIC IC25  = 769.78
%   BIC IC75  = 772.80
%   BIC total = 2257.76

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {A2780_IC50_10y_raw, A2780_IC25_10y_raw, A2780_IC75_10y_raw};

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(A2780_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_mult_alpha_t0(malphat_hat(1), malphat_hat(2), A2780_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(A2780_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 1400000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (n=%.2f, a=%.2f):\n', malphat_hat(1), malphat_hat(2)));

saveas(gcf, fullfile(outdir,'A2780_treated_n-alphaandbeta_multdrug_alphaT0_constD_10days.png'));

writematrix(S_raw_all{1,1}, fullfile(outdir, 'A2780_IC50treated_n-alphaandbeta_multdrug_alphaT0_constD_10days.csv'));
writematrix(S_raw_all{1,2}, fullfile(outdir, 'A2780_IC25treated_n-alphaandbeta_multdrug_alphaT0_constD_10days.csv'));
writematrix(S_raw_all{1,3}, fullfile(outdir, 'A2780_IC75treated_n-alphaandbeta_multdrug_alphaT0_constD_10days.csv'));

%MODEL 8 (Equation 11)
%%%Estimate for n and alpha - A2780 (constant D) multiplicative hill term,
%%%alpha is time dependent, drug kill is alpha at t0 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated';

%Import all data
A2780_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC50_1um_30k_cellsinwell.csv');
A2780_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC75_0.62um_30k_cellsinwell.csv');
A2780_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC25_1.47um_30k_cellsinwell.csv');

A2780_t_days = A2780_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
A2780_IC50_y_raw = [A2780_IC50_data.Rep1, A2780_IC50_data.Rep2, A2780_IC50_data.Rep3];
A2780_IC25_y_raw = [A2780_IC25_data.Rep1, A2780_IC25_data.Rep2, A2780_IC25_data.Rep3];
A2780_IC75_y_raw = [A2780_IC75_data.Rep1, A2780_IC75_data.Rep2, A2780_IC75_data.Rep3];

%downsample data
A2780_10t_days = A2780_t_days(1:11);
A2780_IC25_10y_raw = A2780_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
A2780_IC50_10y_raw = A2780_IC50_y_raw(1:11,:);
A2780_IC75_10y_raw = A2780_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 1.00; 
ce_IC25 = 0.62;
ce_IC75 = 1.47;

%Define logistic growth constants 
rS = 0.8547; %growth rate (day^-1)
KS = 1114404.9924;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 0.9960177; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [2; 0.5];      % [n; alpha]
lb = [0.01; 0.01];   % reasonable bounds
ub = [20; 20];

res_fun_nandalpha_constD_mult_alpha_talpha = @(p) all_datasets_residual_reps_n_and_alpha_constD_mult_alpha_talpha(p, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

%run globalsearch with 200 starts
problem_nandalpha_constD_mult_alpha_talpha = createOptimProblem(...
    'lsqnonlin', ... %can't use lsqcurvefit because it expects x and y arguments
    'x0', p0, ...
    'objective', res_fun_nandalpha_constD_mult_alpha_talpha, ...
    'lb', lb, ...
    'ub', ub);

ms_nandalpha_constD_multi_alpha_talpha = MultiStart('UseParallel', true, 'Display', 'iter', 'StartPointsToRun','bounds-ineqs');  % Set true if desired
nStarts = 200; %more start to ensure reproducibility 
[malphat_hat, resnorm_malphat] = run(ms_nandalpha_constD_multi_alpha_talpha, problem_nandalpha_constD_mult_alpha_talpha, nStarts);
fprintf('Fitted n, and alpha - mult drug - alpha T t0=alpha - constD A2780 (14 days)\n')
fprintf('Fitted Hill coefficient n = %.4f\n', malphat_hat(1));
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(2));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% 43 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted n, and alpha - mult drug - alpha T t0=alpha - constD A2780 (14 days)
% Fitted Hill coefficient n = 19.9992
% Fitted alpha term a = 0.0855
% Sum of squared residuals = 2003929715078.65

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_mult_alpha_talpha(malphat_hat, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

fprintf('Model (n=%.2f, a=%.2f):\n', malphat_hat(1), malphat_hat(2));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (n=20.00, a=0.09):
%   BIC IC50  = 754.84
%   BIC IC25  = 794.07
%   BIC IC75  = 804.14
%   BIC total = 2353.05

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {A2780_IC50_10y_raw, A2780_IC25_10y_raw, A2780_IC75_10y_raw};

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(A2780_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_mult_alpha_talpha(malphat_hat(1), malphat_hat(2), A2780_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(A2780_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 1400000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (n=%.2f, a=%.2f):\n', malphat_hat(1), malphat_hat(2)));

saveas(gcf, fullfile(outdir,'A2780_treated_n-alphaandbeta_multdrug_alphaTalpha_constD_10days.png'));

writematrix(S_raw_all{1,1}, fullfile(outdir, 'A2780_IC50treated_n-alphaandbeta_multdrug_alphaTalpha_constD_10days.csv'));
writematrix(S_raw_all{1,2}, fullfile(outdir, 'A2780_IC25treated_n-alphaandbeta_multdrug_alphaTalpha_constD_10days.csv'));
writematrix(S_raw_all{1,3}, fullfile(outdir, 'A2780_IC75treated_n-alphaandbeta_multdrug_alphaTalpha_constD_10days.csv'));


%MODEL 10 (Equation 12)
%%%Estimate for n, and alpha - A2780 (constant D) additive hill term,
%%%alpha is time dependent, t0=alpha 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated';

%Import all data
A2780_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC50_1um_30k_cellsinwell.csv');
A2780_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC75_0.62um_30k_cellsinwell.csv');
A2780_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC25_1.47um_30k_cellsinwell.csv');

A2780_t_days = A2780_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
A2780_IC50_y_raw = [A2780_IC50_data.Rep1, A2780_IC50_data.Rep2, A2780_IC50_data.Rep3];
A2780_IC25_y_raw = [A2780_IC25_data.Rep1, A2780_IC25_data.Rep2, A2780_IC25_data.Rep3];
A2780_IC75_y_raw = [A2780_IC75_data.Rep1, A2780_IC75_data.Rep2, A2780_IC75_data.Rep3];

%downsample data
A2780_10t_days = A2780_t_days(1:11);
A2780_IC25_10y_raw = A2780_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
A2780_IC50_10y_raw = A2780_IC50_y_raw(1:11,:);
A2780_IC75_10y_raw = A2780_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 1.00; 
ce_IC25 = 0.62;
ce_IC75 = 1.47;

%Define logistic growth constants 
rS = 0.8547; %growth rate (day^-1)
KS = 1114404.9924;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 0.9960177; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [3.0; 0.1];      % [n; alpha; beta]
lb = [0.01; 0.01];   % reasonable bounds
ub = [20.0; 20.0]; %change ub for the last one 

res_fun_nandalpha_constD_add_alpha_talpha = @(p) all_datasets_residual_reps_n_and_alpha_constD_add_alpha_talpha(p, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

%run globalsearch with 200 starts
problem_nandalpha_constD_add_alpha_talpha = createOptimProblem(...
    'lsqnonlin', ... %can't use lsqcurvefit because it expects x and y arguments
    'x0', p0, ...
    'objective', res_fun_nandalpha_constD_add_alpha_talpha, ...
    'lb', lb, ...
    'ub', ub);

ms_nandalpha_constD_add_alpha_talpha = MultiStart('UseParallel', true, 'Display', 'iter', 'StartPointsToRun','bounds-ineqs');  % Set true if desired
nStarts = 200; %more start to ensure reproducibility 
[aalphat_hat, resnorm_aalphat] = run(ms_nandalpha_constD_add_alpha_talpha, problem_nandalpha_constD_add_alpha_talpha, nStarts);
fprintf('Fitted n and alpha - add drug - alpha T - t0=alpha - constD A2780 (10 days)\n')
fprintf('Fitted Hill coefficient n = %.4f\n', aalphat_hat(1));
fprintf('Fitted alpha term a = %.4f\n', aalphat_hat(2));
fprintf('Sum of squared residuals = %.2f\n', resnorm_aalphat);
% 45 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted n and alpha - add drug - alpha T - t0=alpha - constD A2780 (10 days)
% Fitted Hill coefficient n = 3.7142, intial guess 2 and 0.5
% Fitted alpha term a = 0.1067
% Sum of squared residuals = 786875277563.51

% 63 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted n and alpha - add drug - alpha T - t0=alpha - constD A2780 (10 days)
% Fitted Hill coefficient n = 3.7145, initial guesses 3 and 0.1
% Fitted alpha term a = 0.1067
% Sum of squared residuals = 786875177727.52

% 75 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted n and alpha - add drug - alpha T - t0=alpha - constD A2780 (10 days)
% Fitted Hill coefficient n = 3.7156
% Fitted alpha term a = 0.1067
% Sum of squared residuals = 786875112269.50

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_add_alpha_talpha(aalphat_hat, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

fprintf('Model (n=%.2f, a=%.2f):\n', aalphat_hat(1), aalphat_hat(2));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (n=3.71, a=0.11):
%   BIC IC50  = 720.42
%   BIC IC25  = 754.83
%   BIC IC75  = 778.99
%   BIC total = 2254.24

% Model (n=3.71, a=0.11):
%   BIC IC50  = 720.42
%   BIC IC25  = 754.83
%   BIC IC75  = 778.99
%   BIC total = 2254.24

% Model (n=3.72, a=0.11):
%   BIC IC50  = 720.41
%   BIC IC25  = 754.83
%   BIC IC75  = 778.99
%   BIC total = 2254.23

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {A2780_IC50_10y_raw, A2780_IC25_10y_raw, A2780_IC75_10y_raw};

model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(A2780_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_add_alpha_talpha(aalphat_hat(1), aalphat_hat(2), A2780_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    model_fit_all{i} = S_fit_mean;
    plot(A2780_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 1400000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (n=%.2f, a=%.2f):\n', aalphat_hat(1), aalphat_hat(2)));

saveas(gcf, fullfile(outdir,'A2780_treated_nandalpha_adddrug_alphaT_t0alpha_constD_10days.png'));

writematrix([A2780_10t_days, model_fit_all{1}], fullfile(outdir,'A2780_IC50treated_nandalpha_adddrug_alphaT_t0alpha_constD_10days_model10.csv'));
writematrix([A2780_10t_days, model_fit_all{2}], fullfile(outdir,'A2780_IC25treated_nandalpha_adddrug_alphaT_t0alpha_constD_10days_model10.csv'));
writematrix([A2780_10t_days, model_fit_all{3}], fullfile(outdir,'A2780_IC75treated_nandalpha_adddrug_alphaT_t0alpha_constD_10days_model10.csv'));

writematrix(S_raw_all{1,1}, fullfile(outdir, 'A2780_IC50treated_nandalpha_adddrug_alphaT_t0alpha_constD_10days.csv'));
writematrix(S_raw_all{1,2}, fullfile(outdir, 'A2780_IC25treated_nandalpha_adddrug_alphaT_t0alpha_constD_10days.csv'));
writematrix(S_raw_all{1,3}, fullfile(outdir, 'A2780_IC75treated_nandalpha_adddrug_alphaT_t0alpha_constD_10days.csv'));

%MODEL 13 (Equation 8)
%%%Estimate for n and alpha - A2780 (constant D) additive hill term,
%%%alpha is time dependent, drug kill is 0 at time 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated';

%Import all data
A2780_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC50_1um_30k_cellsinwell.csv');
A2780_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC75_0.62um_30k_cellsinwell.csv');
A2780_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC25_1.47um_30k_cellsinwell.csv');

A2780_t_days = A2780_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
A2780_IC50_y_raw = [A2780_IC50_data.Rep1, A2780_IC50_data.Rep2, A2780_IC50_data.Rep3];
A2780_IC25_y_raw = [A2780_IC25_data.Rep1, A2780_IC25_data.Rep2, A2780_IC25_data.Rep3];
A2780_IC75_y_raw = [A2780_IC75_data.Rep1, A2780_IC75_data.Rep2, A2780_IC75_data.Rep3];

%downsample data
A2780_10t_days = A2780_t_days(1:11);
A2780_IC25_10y_raw = A2780_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
A2780_IC50_10y_raw = A2780_IC50_y_raw(1:11,:);
A2780_IC75_10y_raw = A2780_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 1.00; 
ce_IC25 = 0.62;
ce_IC75 = 1.47;

%Define logistic growth constants 
rS = 0.8547; %growth rate (day^-1)
KS = 1114404.9924;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 0.9960177; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [1; 1];      % [n; alpha] %initially 2 and 0.5 same results for those initial guesses
lb = [0.01; 0.01];   % reasonable bounds
ub = [20; 20];

res_fun_nandalpha_constD_add_alpha_t0 = @(p) all_datasets_residual_reps_n_and_alpha_constD_add_alpha_t0(p, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

%run globalsearch with 200 starts
problem_nandalpha_constD_add_alpha_t0 = createOptimProblem(...
    'lsqnonlin', ... %can't use lsqcurvefit because it expects x and y arguments
    'x0', p0, ...
    'objective', res_fun_nandalpha_constD_add_alpha_t0, ...
    'lb', lb, ...
    'ub', ub);

ms_nandalpha_constD_add_alpha_t0 = MultiStart('UseParallel', true, 'Display', 'iter', 'StartPointsToRun','bounds-ineqs');  % Set true if desired
nStarts = 200; %more start to ensure reproducibility 
[malphat_hat, resnorm_malphat] = run(ms_nandalpha_constD_add_alpha_t0, problem_nandalpha_constD_add_alpha_t0, nStarts);
fprintf('Fitted n, and alpha - add drug - alpha T t0=0 - constD A2780 (14 days)\n')
fprintf('Fitted Hill coefficient n = %.4f\n', malphat_hat(1));
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(2));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% 198 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted n, and alpha - add drug - alpha T t0=0 - constD A2780 (14 days)
% Fitted Hill coefficient n = 3.7779
% Fitted alpha term a = 0.1290
% Sum of squared residuals = 687683459491.16

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_log_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_nandalpha_constD_add_alpha_t0, malphat_hat, lb, ub);

% Compute 95% confidence intervals
ci = nlparci(p_log_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95 Confidence Intervals:\n');
fprintf('n: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
fprintf('a: [%.4f, %.4f]\n', ci(2,1), ci(2,2));
% Local minimum possible.
% lsqnonlin stopped because the size of the current step is less than
% the value of the step size tolerance.
% <stopping criteria details>
% 95 Confidence Intervals:
% n: [3.2408, 4.3149] (4.3149-3.2408)/3.7779 = 0.2843
% a: [0.1214, 0.1367] (0.1367-0.1214)/0.1290 = 0.1186

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_add_alpha_t0(malphat_hat, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

fprintf('Model (n=%.2f, a=%.2f):\n', malphat_hat(1), malphat_hat(2));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (n=3.78, a=0.13):
%   BIC IC50  = 731.19
%   BIC IC25  = 753.34
%   BIC IC75  = 769.40
%   BIC total = 2253.92

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {A2780_IC50_10y_raw, A2780_IC25_10y_raw, A2780_IC75_10y_raw};

model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(A2780_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_add_alpha_t0(malphat_hat(1), malphat_hat(2), A2780_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    model_fit_all{i} = S_fit_mean;
    plot(A2780_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 1400000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (n=%.2f, a=%.2f):\n', malphat_hat(1), malphat_hat(2)));

saveas(gcf, fullfile(outdir,'A2780_treated_nandalpha_adddrug_alphaT0_constD_10days_Model13.png'));

writematrix([A2780_10t_days, model_fit_all{1}], fullfile(outdir,'A2780_IC50treated_nandalpha_adddrug_alphaT0_constD_10days_Model13.csv'));
writematrix([A2780_10t_days, model_fit_all{2}], fullfile(outdir,'A2780_IC25treated_nandalpha_adddrug_alphaT0_constD_10days_Model13.csv'));
writematrix([A2780_10t_days, model_fit_all{3}], fullfile(outdir,'A2780_IC75treated_nandalpha_adddrug_alphaT0_constD_10days_Model13.csv'));

writematrix(S_raw_all{1,1}, fullfile(outdir, 'A2780_IC50treated_nandalpha_adddrug_alphaT0_constD_10days_Model13.csv'));
writematrix(S_raw_all{1,2}, fullfile(outdir, 'A2780_IC25treated_nandalpha_adddrug_alphaT0_constD_10days_Model13.csv'));
writematrix(S_raw_all{1,3}, fullfile(outdir, 'A2780_IC75treated_nandalpha_adddrug_alphaT0_constD_10days_Model13.csv'));


%R squared calculation
clear all;
A2780_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC50_1um_30k_cellsinwell.csv');
A2780_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC75_0.62um_30k_cellsinwell.csv');
A2780_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC25_1.47um_30k_cellsinwell.csv');

A2780_t_days = A2780_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
A2780_IC50_y_raw = [A2780_IC50_data.Rep1, A2780_IC50_data.Rep2, A2780_IC50_data.Rep3];
A2780_IC25_y_raw = [A2780_IC25_data.Rep1, A2780_IC25_data.Rep2, A2780_IC25_data.Rep3];
A2780_IC75_y_raw = [A2780_IC75_data.Rep1, A2780_IC75_data.Rep2, A2780_IC75_data.Rep3];

%downsample data
A2780_IC25_10y_raw = A2780_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
A2780_IC50_10y_raw = A2780_IC50_y_raw(1:11,:);
A2780_IC75_10y_raw = A2780_IC75_y_raw(1:11,:);

A2780_IC25_10y_all = A2780_IC25_10y_raw(:);
A2780_IC50_10y_all = A2780_IC50_10y_raw(:);
A2780_IC75_10y_all = A2780_IC75_10y_raw(:);

%IC25
A2780_IC25_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated/A2780_IC25treated_nandalpha_adddrug_alphaT0_constD_10days_Model13.csv');
A2780_IC25_predicted = table2array(A2780_IC25_predicted(:,2)); %only the second column 
A2780_IC25_predicted = repmat(A2780_IC25_predicted,3,1); %repeat 3 times 

SS_res_IC25 = sum((A2780_IC25_10y_all - A2780_IC25_predicted).^2);
SS_tot_IC25 = sum((A2780_IC25_10y_all - mean(A2780_IC25_10y_all)).^2);
R2_IC25 = 1 - SS_res_IC25 / SS_tot_IC25;

fprintf('R-squared for IC25: %.4f\n', R2_IC25);
%R-squared for IC25: 0.9695

%IC50
A2780_IC50_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated/A2780_IC50treated_nandalpha_adddrug_alphaT0_constD_10days_Model13.csv');
A2780_IC50_predicted = table2array(A2780_IC50_predicted(:,2)); %only the second column 
A2780_IC50_predicted = repmat(A2780_IC50_predicted,3,1); %repeat 3 times 

SS_res_IC50 = sum((A2780_IC50_10y_all - A2780_IC50_predicted).^2);
SS_tot_IC50 = sum((A2780_IC50_10y_all - mean(A2780_IC50_10y_all)).^2);
R2_IC50 = 1 - SS_res_IC50 / SS_tot_IC50;

fprintf('R-squared for IC50: %.4f\n', R2_IC50);
%R-squared for IC50: 0.9404

%IC75
A2780_IC75_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated/A2780_IC75treated_nandalpha_adddrug_alphaT0_constD_10days_Model13.csv');
A2780_IC75_predicted = table2array(A2780_IC75_predicted(:,2)); %only the second column 
A2780_IC75_predicted = repmat(A2780_IC75_predicted,3,1); %repeat 3 times 

SS_res_IC75 = sum((A2780_IC75_10y_all - A2780_IC75_predicted).^2);
SS_tot_IC75 = sum((A2780_IC75_10y_all - mean(A2780_IC75_10y_all)).^2);
R2_IC75 = 1 - SS_res_IC75 / SS_tot_IC75;

fprintf('R-squared for IC75: %.4f\n', R2_IC75);
%R-squared for IC75: 0.8959



%MODEL 14 (Equation 10)
%%%Estimate for alpha only - A2780 (constant D) additive hill term,
%%%alpha is time dependent, drug kill is 0 at time 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated';

%Import all data
A2780_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC50_1um_30k_cellsinwell.csv');
A2780_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC75_0.62um_30k_cellsinwell.csv');
A2780_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC25_1.47um_30k_cellsinwell.csv');

A2780_t_days = A2780_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
A2780_IC50_y_raw = [A2780_IC50_data.Rep1, A2780_IC50_data.Rep2, A2780_IC50_data.Rep3];
A2780_IC25_y_raw = [A2780_IC25_data.Rep1, A2780_IC25_data.Rep2, A2780_IC25_data.Rep3];
A2780_IC75_y_raw = [A2780_IC75_data.Rep1, A2780_IC75_data.Rep2, A2780_IC75_data.Rep3];

%downsample data
A2780_10t_days = A2780_t_days(1:11);
A2780_IC25_10y_raw = A2780_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
A2780_IC50_10y_raw = A2780_IC50_y_raw(1:11,:);
A2780_IC75_10y_raw = A2780_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 1.00; 
ce_IC25 = 0.62;
ce_IC75 = 1.47;

%Define logistic growth constants 
rS = 0.8547; %growth rate (day^-1)
KS = 1114404.9924;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 0.9960177; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [0.5];      % [alpha]
lb = [0.0001];   % reasonable bounds
ub = [20];

res_fun_nandalpha_constD_add_alphaonly_t0 = @(p) all_datasets_residual_reps_n_and_alpha_constD_add_alphaonly_t0(p, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

%run globalsearch with 200 starts
problem_nandalpha_constD_add_alphaonly_t0 = createOptimProblem(...
    'lsqnonlin', ... %can't use lsqcurvefit because it expects x and y arguments
    'x0', p0, ...
    'objective', res_fun_nandalpha_constD_add_alphaonly_t0, ...
    'lb', lb, ...
    'ub', ub);

ms_nandalpha_constD_add_alphaonly_t0 = MultiStart('UseParallel', true, 'Display', 'iter', 'StartPointsToRun','bounds-ineqs');  % Set true if desired
nStarts = 200; %more start to ensure reproducibility 
[malphat_hat, resnorm_malphat] = run(ms_nandalpha_constD_add_alphaonly_t0, problem_nandalpha_constD_add_alphaonly_t0, nStarts);
fprintf('Fitted alpha only - add drug - alpha T t0=0 - constD A2780 (14 days)\n')
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted alpha only - add drug - alpha T t0=0 - constD A2780 (14 days)
% Fitted alpha term a = 0.1248
% Sum of squared residuals = 2287132727028.39

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_add_alphaonly_t0(malphat_hat, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

fprintf('Model (a=%.2f):\n', malphat_hat(1));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (a=0.12):
%   BIC IC50  = 726.09
%   BIC IC25  = 806.71
%   BIC IC75  = 798.33
%   BIC total = 2331.13

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {A2780_IC50_10y_raw, A2780_IC25_10y_raw, A2780_IC75_10y_raw};
model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(A2780_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_add_alphaonly_t0(malphat_hat(1), A2780_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    model_fit_all{i} = S_fit_mean;
    plot(A2780_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 1400000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (a=%.2f):\n', malphat_hat(1)));

saveas(gcf, fullfile(outdir,'A2780_treated_alphaonly_adddrug_alphaT0_constD_10days_Model14.png'));

% writematrix(S_raw_all{1,1}, fullfile(outdir, 'A2780_IC50treated_alphaonly_adddrug_alphaT0_constD_10days_Model14.csv'));
% writematrix(S_raw_all{1,2}, fullfile(outdir, 'A2780_IC25treated_alphaonly_adddrug_alphaT0_constD_10days_Model14.csv'));
% writematrix(S_raw_all{1,3}, fullfile(outdir, 'A2780_IC75treated_alphaonly_adddrug_alphaT0_constD_10days_Model14.csv'));

writematrix([A2780_10t_days, model_fit_all{1}], fullfile(outdir,'A2780_IC50treated_alphaonly_adddrug_alphaT0_constD_10days_Model14.csv'));
writematrix([A2780_10t_days, model_fit_all{2}], fullfile(outdir,'A2780_IC25treated_alphaonly_adddrug_alphaT0_constD_10days_Model14.csv'));
writematrix([A2780_10t_days, model_fit_all{3}], fullfile(outdir,'A2780_IC75treated_alphaonly_adddrug_alphaT0_constD_10days_Model14.csv'));

%MODEL 15 (Equation 9)
%%%Estimate for alpha only - A2780 (constant D) multiplicative hill term,
%%%alpha is time dependent, drug kill is 0 at time 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated';

%Import all data
A2780_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC50_1um_30k_cellsinwell.csv');
A2780_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC75_0.62um_30k_cellsinwell.csv');
A2780_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC25_1.47um_30k_cellsinwell.csv');

A2780_t_days = A2780_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
A2780_IC50_y_raw = [A2780_IC50_data.Rep1, A2780_IC50_data.Rep2, A2780_IC50_data.Rep3];
A2780_IC25_y_raw = [A2780_IC25_data.Rep1, A2780_IC25_data.Rep2, A2780_IC25_data.Rep3];
A2780_IC75_y_raw = [A2780_IC75_data.Rep1, A2780_IC75_data.Rep2, A2780_IC75_data.Rep3];

%downsample data
A2780_10t_days = A2780_t_days(1:11);
A2780_IC25_10y_raw = A2780_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
A2780_IC50_10y_raw = A2780_IC50_y_raw(1:11,:);
A2780_IC75_10y_raw = A2780_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 1.00; 
ce_IC25 = 0.62;
ce_IC75 = 1.47;

%Define logistic growth constants 
rS = 0.8547; %growth rate (day^-1)
KS = 1114404.9924;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 0.9960177; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [0.5];      % [alpha]
lb = [0.0001];   % reasonable bounds
ub = [20];

res_fun_nandalpha_constD_mult_alphaonly_t0 = @(p) all_datasets_residual_reps_n_and_alpha_constD_mult_alphaonly_t0(p, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

%run globalsearch with 200 starts
problem_nandalpha_constD_mult_alphaonly_t0 = createOptimProblem(...
    'lsqnonlin', ... %can't use lsqcurvefit because it expects x and y arguments
    'x0', p0, ...
    'objective', res_fun_nandalpha_constD_mult_alphaonly_t0, ...
    'lb', lb, ...
    'ub', ub);

ms_nandalpha_constD_mult_alphaonly_t0 = MultiStart('UseParallel', true, 'Display', 'iter', 'StartPointsToRun','bounds-ineqs');  % Set true if desired
nStarts = 200; %more start to ensure reproducibility 
[malphat_hat, resnorm_malphat] = run(ms_nandalpha_constD_mult_alphaonly_t0, problem_nandalpha_constD_mult_alphaonly_t0, nStarts);
fprintf('Fitted alpha only - mult drug - alpha T t0=0 - constD A2780 (10 days)\n')
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted alpha only - mult drug - alpha T t0=0 - constD A2780 (10 days)
% Fitted alpha term a = 0.2447
% Sum of squared residuals = 1387254643073.12

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_mult_alphaonly_t0(malphat_hat, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

fprintf('Model (a=%.2f):\n', malphat_hat(1));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (a=0.24):
%   BIC IC50  = 712.25
%   BIC IC25  = 789.80
%   BIC IC75  = 782.05
%   BIC total = 2284.10

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {A2780_IC50_10y_raw, A2780_IC25_10y_raw, A2780_IC75_10y_raw};

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(A2780_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_mult_alphaonly_t0(malphat_hat(1), A2780_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(A2780_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 1400000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (a=%.2f):\n', malphat_hat(1)));

saveas(gcf, fullfile(outdir,'A2780_treated_alphaonly_multdrug_alphaT0_constD_10days_Model15.png'));

writematrix(S_raw_all{1,1}, fullfile(outdir, 'A2780_IC50treated_alphaonly_multdrug_alphaT0_constD_10days_Model15.csv'));
writematrix(S_raw_all{1,2}, fullfile(outdir, 'A2780_IC25treated_alphaonly_multdrug_alphaT0_constD_10days_Model15.csv'));
writematrix(S_raw_all{1,3}, fullfile(outdir, 'A2780_IC75treated_alphaonly_multdrug_alphaT0_constD_10days_Model15.csv'));

%MODEL 16 (Equation 14)
%%%Estimate for alpha only - A2780 (constant D) additive hill term,
%%%alpha is time dependent, drug kill is alpha at time 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated';

%Import all data
A2780_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC50_1um_30k_cellsinwell.csv');
A2780_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC75_0.62um_30k_cellsinwell.csv');
A2780_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC25_1.47um_30k_cellsinwell.csv');

A2780_t_days = A2780_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
A2780_IC50_y_raw = [A2780_IC50_data.Rep1, A2780_IC50_data.Rep2, A2780_IC50_data.Rep3];
A2780_IC25_y_raw = [A2780_IC25_data.Rep1, A2780_IC25_data.Rep2, A2780_IC25_data.Rep3];
A2780_IC75_y_raw = [A2780_IC75_data.Rep1, A2780_IC75_data.Rep2, A2780_IC75_data.Rep3];

%downsample data
A2780_10t_days = A2780_t_days(1:11);
A2780_IC25_10y_raw = A2780_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
A2780_IC50_10y_raw = A2780_IC50_y_raw(1:11,:);
A2780_IC75_10y_raw = A2780_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 1.00; 
ce_IC25 = 0.62;
ce_IC75 = 1.47;

%Define logistic growth constants 
rS = 0.8547; %growth rate (day^-1)
KS = 1114404.9924;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 0.9960177; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [0.5];      % [alpha]
lb = [0.0001];   % reasonable bounds
ub = [20];

res_fun_nandalpha_constD_add_alphaonly_talpha = @(p) all_datasets_residual_reps_n_and_alpha_constD_add_alphaonly_talpha(p, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

%run globalsearch with 200 starts
problem_nandalpha_constD_add_alphaonly_talpha = createOptimProblem(...
    'lsqnonlin', ... %can't use lsqcurvefit because it expects x and y arguments
    'x0', p0, ...
    'objective', res_fun_nandalpha_constD_add_alphaonly_talpha, ...
    'lb', lb, ...
    'ub', ub);

ms_nandalpha_constD_add_alphaonly_talpha = MultiStart('UseParallel', true, 'Display', 'iter', 'StartPointsToRun','bounds-ineqs');  % Set true if desired
nStarts = 200; %more start to ensure reproducibility 
[malphat_hat, resnorm_malphat] = run(ms_nandalpha_constD_add_alphaonly_talpha, problem_nandalpha_constD_add_alphaonly_talpha, nStarts);
fprintf('Fitted alpha only - add drug - alpha T t0=alpha - constD A2780 (10 days)\n')
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted alpha only - add drug - alpha T t0=alpha - constD A2780 (10 days)
% Fitted alpha term a = 0.1044
% Sum of squared residuals = 2353616407928.27

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_add_alphaonly_talpha(malphat_hat, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

fprintf('Model (a=%.2f):\n', malphat_hat(1));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (a=0.10):
%   BIC IC50  = 715.87
%   BIC IC25  = 806.95
%   BIC IC75  = 801.18
%   BIC total = 2324.00

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {A2780_IC50_10y_raw, A2780_IC25_10y_raw, A2780_IC75_10y_raw};

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(A2780_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_add_alphaonly_talpha(malphat_hat(1), A2780_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(A2780_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 1400000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (a=%.2f):\n', malphat_hat(1)));

saveas(gcf, fullfile(outdir,'A2780_treated_alphaonly_adddrug_alphaTalpha_constD_10days_Model16.png'));

writematrix(S_raw_all{1,1}, fullfile(outdir, 'A2780_IC50treated_alphaonly_adddrug_alphaTalpha_constD_10days_Model16.csv'));
writematrix(S_raw_all{1,2}, fullfile(outdir, 'A2780_IC25treated_alphaonly_adddrug_alphaTalpha_constD_10days_Model16.csv'));
writematrix(S_raw_all{1,3}, fullfile(outdir, 'A2780_IC75treated_alphaonly_adddrug_alphaTalpha_constD_10days_Model16.csv'));

%MODEL 17 (Equation 13)
%%%Estimate for alpha only - A2780 (constant D) multiplicative hill term,
%%%alpha is time dependent, drug kill is alpha at time 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated';

%Import all data
A2780_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC50_1um_30k_cellsinwell.csv');
A2780_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC75_0.62um_30k_cellsinwell.csv');
A2780_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC25_1.47um_30k_cellsinwell.csv');

A2780_t_days = A2780_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
A2780_IC50_y_raw = [A2780_IC50_data.Rep1, A2780_IC50_data.Rep2, A2780_IC50_data.Rep3];
A2780_IC25_y_raw = [A2780_IC25_data.Rep1, A2780_IC25_data.Rep2, A2780_IC25_data.Rep3];
A2780_IC75_y_raw = [A2780_IC75_data.Rep1, A2780_IC75_data.Rep2, A2780_IC75_data.Rep3];

%downsample data
A2780_10t_days = A2780_t_days(1:11);
A2780_IC25_10y_raw = A2780_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
A2780_IC50_10y_raw = A2780_IC50_y_raw(1:11,:);
A2780_IC75_10y_raw = A2780_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 1.00; 
ce_IC25 = 0.62;
ce_IC75 = 1.47;

%Define logistic growth constants 
rS = 0.8547; %growth rate (day^-1)
KS = 1114404.9924;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 0.9960177; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [0.5];      % [alpha]
lb = [0.0001];   % reasonable bounds
ub = [20];

res_fun_nandalpha_constD_mult_alphaonly_talpha = @(p) all_datasets_residual_reps_n_and_alpha_constD_mult_alphaonly_talpha(p, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

%run globalsearch with 200 starts
problem_nandalpha_constD_mult_alphaonly_talpha = createOptimProblem(...
    'lsqnonlin', ... %can't use lsqcurvefit because it expects x and y arguments
    'x0', p0, ...
    'objective', res_fun_nandalpha_constD_mult_alphaonly_talpha, ...
    'lb', lb, ...
    'ub', ub);

ms_nandalpha_constD_mult_alphaonly_talpha = MultiStart('UseParallel', true, 'Display', 'iter', 'StartPointsToRun','bounds-ineqs');  % Set true if desired
nStarts = 200; %more start to ensure reproducibility 
[malphat_hat, resnorm_malphat] = run(ms_nandalpha_constD_mult_alphaonly_talpha, problem_nandalpha_constD_mult_alphaonly_talpha, nStarts);
fprintf('Fitted alpha only - mult drug - alpha T t0=alpha - constD A2780 (10 days)\n')
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted alpha only - mult drug - alpha T t0=alpha - constD A2780 (10 days)
% Fitted alpha term a = 0.1970
% Sum of squared residuals = 1695405401379.65

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_mult_alphaonly_talpha(malphat_hat, rS, KS, IC50S, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

fprintf('Model (a=%.2f):\n', malphat_hat(1));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (a=0.20):
%   BIC IC50  = 722.65
%   BIC IC25  = 795.57
%   BIC IC75  = 789.24
%   BIC total = 2307.46

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {A2780_IC50_10y_raw, A2780_IC25_10y_raw, A2780_IC75_10y_raw};

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(A2780_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_mult_alphaonly_talpha(malphat_hat(1), A2780_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(A2780_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 1400000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (a=%.2f):\n', malphat_hat(1)));

saveas(gcf, fullfile(outdir,'A2780_treated_alphaonly_multdrug_alphaTalpha_constD_10days_Model17.png'));

writematrix(S_raw_all{1,1}, fullfile(outdir, 'A2780_IC50treated_alphaonly_multdrug_alphaTalpha_constD_10days_Model17.csv'));
writematrix(S_raw_all{1,2}, fullfile(outdir, 'A2780_IC25treated_alphaonly_multdrug_alphaTalpha_constD_10days_Model17.csv'));
writematrix(S_raw_all{1,3}, fullfile(outdir, 'A2780_IC75treated_alphaonly_multdrug_alphaTalpha_constD_10days_Model17.csv'));
%this one looks a little better for IC50 

%MODEL 18 (Equation 17)
%%%Estimate for alpha only - A2780 (constant D) multiplicative hill term max,
%%%alpha is time dependent, drug kill is 0 at time 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated';

%Import all data
A2780_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC50_1um_30k_cellsinwell.csv');
A2780_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC75_0.62um_30k_cellsinwell.csv');
A2780_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC25_1.47um_30k_cellsinwell.csv');

A2780_t_days = A2780_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
A2780_IC50_y_raw = [A2780_IC50_data.Rep1, A2780_IC50_data.Rep2, A2780_IC50_data.Rep3];
A2780_IC25_y_raw = [A2780_IC25_data.Rep1, A2780_IC25_data.Rep2, A2780_IC25_data.Rep3];
A2780_IC75_y_raw = [A2780_IC75_data.Rep1, A2780_IC75_data.Rep2, A2780_IC75_data.Rep3];

%downsample data
A2780_10t_days = A2780_t_days(1:11);
A2780_IC25_10y_raw = A2780_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
A2780_IC50_10y_raw = A2780_IC50_y_raw(1:11,:);
A2780_IC75_10y_raw = A2780_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 1.00; 
ce_IC25 = 0.62;
ce_IC75 = 1.47;

%Define logistic growth constants 
rS = 0.8547; %growth rate (day^-1)
KS = 1114404.9924;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 0.9960177; %the IC50 for sensitive cells, after 72hr exposure
Dmax = 27.15376;

%Estimate for parameters
p0 = [0.5];      % [alpha]
lb = [0.0001];   % reasonable bounds
ub = [20];

res_fun_nandalpha_constD_multmax_alphaonly_t0 = @(p) all_datasets_residual_reps_n_and_alpha_constD_multmax_alphaonly_t0(p, rS, KS, IC50S, Dmax, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

%run globalsearch with 200 starts
problem_nandalpha_constD_multmax_alphaonly_t0 = createOptimProblem(...
    'lsqnonlin', ... %can't use lsqcurvefit because it expects x and y arguments
    'x0', p0, ...
    'objective', res_fun_nandalpha_constD_multmax_alphaonly_t0, ...
    'lb', lb, ...
    'ub', ub);

ms_nandalpha_constD_multmax_alphaonly_t0 = MultiStart('UseParallel', true, 'Display', 'iter', 'StartPointsToRun','bounds-ineqs');  % Set true if desired
nStarts = 200; %more start to ensure reproducibility 
[malphat_hat, resnorm_malphat] = run(ms_nandalpha_constD_multmax_alphaonly_t0, problem_nandalpha_constD_multmax_alphaonly_t0, nStarts);
fprintf('Fitted alpha only - mult drug max - alpha T t0=0 - constD A2780 (10 days)\n')
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted alpha only - mult drug max - alpha T t0=0 - constD A2780 (10 days)
% Fitted alpha term a = 0.0090
% Sum of squared residuals = 1387254856883.53

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_multmax_alphaonly_t0(malphat_hat, rS, KS, IC50S, Dmax, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

fprintf('Model (a=%.2f):\n', malphat_hat(1));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (a=0.01):
%   BIC IC50  = 712.25
%   BIC IC25  = 789.80
%   BIC IC75  = 782.05
%   BIC total = 2284.10

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {A2780_IC50_10y_raw, A2780_IC25_10y_raw, A2780_IC75_10y_raw};

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(A2780_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_multmax_alphaonly_t0(malphat_hat(1), A2780_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S, Dmax);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(A2780_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 1400000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (a=%.2f):\n', malphat_hat(1)));

saveas(gcf, fullfile(outdir,'A2780_treated_alphaonly_multdrugmax_alphaT0_constD_10days_Model18.png'));

writematrix(S_raw_all{1,1}, fullfile(outdir, 'A2780_IC50treated_alphaonly_multdrugmax_alphaT0_constD_10days_Model18.csv'));
writematrix(S_raw_all{1,2}, fullfile(outdir, 'A2780_IC25treated_alphaonly_multdrugmax_alphaT0_constD_10days_Model18.csv'));
writematrix(S_raw_all{1,3}, fullfile(outdir, 'A2780_IC75treated_alphaonly_multdrugmax_alphaT0_constD_10days_Model18.csv'));
 
%MODEL 19 (Equation 18)
%%%Estimate for alpha only - A2780 (constant D) additive hill term max,
%%%alpha is time dependent, drug kill is 0 at time 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated';

%Import all data
A2780_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC50_1um_30k_cellsinwell.csv');
A2780_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC75_0.62um_30k_cellsinwell.csv');
A2780_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC25_1.47um_30k_cellsinwell.csv');

A2780_t_days = A2780_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
A2780_IC50_y_raw = [A2780_IC50_data.Rep1, A2780_IC50_data.Rep2, A2780_IC50_data.Rep3];
A2780_IC25_y_raw = [A2780_IC25_data.Rep1, A2780_IC25_data.Rep2, A2780_IC25_data.Rep3];
A2780_IC75_y_raw = [A2780_IC75_data.Rep1, A2780_IC75_data.Rep2, A2780_IC75_data.Rep3];

%downsample data
A2780_10t_days = A2780_t_days(1:11);
A2780_IC25_10y_raw = A2780_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
A2780_IC50_10y_raw = A2780_IC50_y_raw(1:11,:);
A2780_IC75_10y_raw = A2780_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 1.00; 
ce_IC25 = 0.62;
ce_IC75 = 1.47;

%Define logistic growth constants 
rS = 0.8547; %growth rate (day^-1)
KS = 1114404.9924;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 0.9960177; %the IC50 for sensitive cells, after 72hr exposure
Dmax = 27.15376;

%Estimate for parameters
p0 = [0.5];      % [alpha]
lb = [0.0001];   % reasonable bounds
ub = [20];

res_fun_nandalpha_constD_addmax_alphaonly_t0 = @(p) all_datasets_residual_reps_n_and_alpha_constD_addmax_alphaonly_t0(p, rS, KS, IC50S, Dmax, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

%run globalsearch with 200 starts
problem_nandalpha_constD_addmax_alphaonly_t0 = createOptimProblem(...
    'lsqnonlin', ... %can't use lsqcurvefit because it expects x and y arguments
    'x0', p0, ...
    'objective', res_fun_nandalpha_constD_addmax_alphaonly_t0, ...
    'lb', lb, ...
    'ub', ub);

ms_nandalpha_constD_addmax_alphaonly_t0 = MultiStart('UseParallel', true, 'Display', 'iter', 'StartPointsToRun','bounds-ineqs');  % Set true if desired
nStarts = 200; %more start to ensure reproducibility 
[malphat_hat, resnorm_malphat] = run(ms_nandalpha_constD_addmax_alphaonly_t0, problem_nandalpha_constD_addmax_alphaonly_t0, nStarts);
fprintf('Fitted alpha only - add drug max - alpha T t0=0 - constD A2780 (10 days)\n')
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted alpha only - add drug max - alpha T t0=0 - constD A2780 (10 days)
% Fitted alpha term a = 0.0046
% Sum of squared residuals = 2287133184802.88

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_addmax_alphaonly_t0(malphat_hat, rS, KS, IC50S, Dmax, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

fprintf('Model (a=%.2f):\n', malphat_hat(1));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (a=0.00):
%   BIC IC50  = 726.09
%   BIC IC25  = 806.67
%   BIC IC75  = 798.39
%   BIC total = 2331.14

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {A2780_IC50_10y_raw, A2780_IC25_10y_raw, A2780_IC75_10y_raw};

model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(A2780_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_addmax_alphaonly_t0(malphat_hat(1), A2780_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S, Dmax);
    end
    S_fit_mean = mean(S_fit_all, 2);
    model_fit_all{i} = S_fit_mean;
    plot(A2780_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 1400000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (a=%.2f):\n', malphat_hat(1)));

saveas(gcf, fullfile(outdir,'A2780_treated_alphaonly_adddrugmax_alphaT0_constD_10days_Model19.png'));

writematrix([A2780_10t_days, model_fit_all{1}], fullfile(outdir,'A2780_IC50treated_alphaonly_adddrugmax_alphaT0_constD_10days_Model19.csv'));
writematrix([A2780_10t_days, model_fit_all{2}], fullfile(outdir,'A2780_IC25treated_alphaonly_adddrugmax_alphaT0_constD_10days_Model19.csv'));
writematrix([A2780_10t_days, model_fit_all{3}], fullfile(outdir,'A2780_IC75treated_alphaonly_adddrugmax_alphaT0_constD_10days_Model19.csv'));

writematrix(S_raw_all{1,1}, fullfile(outdir, 'A2780_IC50treated_alphaonly_adddrugmax_alphaT0_constD_10days_Model19.csv'));
writematrix(S_raw_all{1,2}, fullfile(outdir, 'A2780_IC25treated_alphaonly_adddrugmax_alphaT0_constD_10days_Model19.csv'));
writematrix(S_raw_all{1,3}, fullfile(outdir, 'A2780_IC75treated_alphaonly_adddrugmax_alphaT0_constD_10days_Model19.csv'));

%MODEL 20 (Equation 19)
%%%Estimate for alpha only - A2780 (constant D) additive hill,
%%%alpha is time dependent, drug kill is 0 at time 0, LINEAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated';
 
%Import all data
A2780_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC50_1um_30k_cellsinwell.csv');
A2780_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC75_0.62um_30k_cellsinwell.csv');
A2780_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC25_1.47um_30k_cellsinwell.csv');

A2780_t_days = A2780_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
A2780_IC50_y_raw = [A2780_IC50_data.Rep1, A2780_IC50_data.Rep2, A2780_IC50_data.Rep3];
A2780_IC25_y_raw = [A2780_IC25_data.Rep1, A2780_IC25_data.Rep2, A2780_IC25_data.Rep3];
A2780_IC75_y_raw = [A2780_IC75_data.Rep1, A2780_IC75_data.Rep2, A2780_IC75_data.Rep3];

%downsample data
A2780_10t_days = A2780_t_days(1:11);
A2780_IC25_10y_raw = A2780_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
A2780_IC50_10y_raw = A2780_IC50_y_raw(1:11,:);
A2780_IC75_10y_raw = A2780_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 1.00; 
ce_IC25 = 0.62;
ce_IC75 = 1.47;

%Define logistic growth constants 
rS = 0.8547; %growth rate (day^-1)
KS = 1114404.9924;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 0.9960177; %the IC50 for sensitive cells, after 72hr exposure
Dmax = 27.15376;

%Estimate for parameters
p0 = [0.5];      % [alpha]
lb = [0.0001];   % reasonable bounds
ub = [20];

res_fun_nandalpha_constD_linear_alphaonly_t0 = @(p) all_datasets_residual_reps_n_and_alpha_constD_linear_alphaonly_t0(p, rS, KS, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

%run globalsearch with 200 starts
problem_nandalpha_constD_linear_alphaonly_t0 = createOptimProblem(...
    'lsqnonlin', ... %can't use lsqcurvefit because it expects x and y arguments
    'x0', p0, ...
    'objective', res_fun_nandalpha_constD_linear_alphaonly_t0, ...
    'lb', lb, ...
    'ub', ub);

ms_nandalpha_constD_linear_alphaonly_t0 = MultiStart('UseParallel', true, 'Display', 'iter', 'StartPointsToRun','bounds-ineqs');  % Set true if desired
nStarts = 200; %more start to ensure reproducibility 
[malphat_hat, resnorm_malphat] = run(ms_nandalpha_constD_linear_alphaonly_t0, problem_nandalpha_constD_linear_alphaonly_t0, nStarts);
fprintf('Fitted alpha only - add drug max - alpha T t0=0 - constD A2780 (10 days)\n')
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted alpha only - add drug max - alpha T t0=0 - constD A2780 (10 days)
% Fitted alpha term a = 0.0612
% Sum of squared residuals = 1271999329713.16

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_linear_alphaonly_t0(malphat_hat, rS, KS, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

fprintf('Model (a=%.2f):\n', malphat_hat(1));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (a=0.06):
%   BIC IC50  = 726.59
%   BIC IC25  = 786.54
%   BIC IC75  = 776.89
%   BIC total = 2290.02

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {A2780_IC50_10y_raw, A2780_IC25_10y_raw, A2780_IC75_10y_raw};

model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(A2780_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_linear_alphaonly_t0(malphat_hat(1), A2780_10t_days, S0, ...
            ce_all(i), rS, KS);
    end
    S_fit_mean = mean(S_fit_all, 2);
    model_fit_all{i} = S_fit_mean;
    plot(A2780_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 1400000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (a=%.2f):\n', malphat_hat(1)));

saveas(gcf, fullfile(outdir,'A2780_treated_alphaonly_adddruglinear_alphaT0_constD_10days_Model20.png'));

writematrix([A2780_10t_days, model_fit_all{1}], fullfile(outdir,'A2780_IC50treated_alphaonly_adddruglinear_alphaT0_constD_10days_Model20.csv'));
writematrix([A2780_10t_days, model_fit_all{2}], fullfile(outdir,'A2780_IC25treated_alphaonly_adddruglinear_alphaT0_constD_10days_Model20.csv'));
writematrix([A2780_10t_days, model_fit_all{3}], fullfile(outdir,'A2780_IC75treated_alphaonly_adddruglinear_alphaT0_constD_10days_Model20.csv'));

%writematrix(S_raw_all{1,1}, fullfile(outdir, 'A2780_IC50treated_alphaonly_adddruglinear_alphaT0_constD_10days_Model20.csv'));
%writematrix(S_raw_all{1,2}, fullfile(outdir, 'A2780_IC25treated_alphaonly_adddruglinear_alphaT0_constD_10days_Model20.csv'));
%writematrix(S_raw_all{1,3}, fullfile(outdir, 'A2780_IC75treated_alphaonly_adddruglinear_alphaT0_constD_10days_Model20.csv'));

%MODEL 21, equivalent to model 15 and 7 but linear same as 20 (Equation 20)
%%%Estimate for alpha only - A2780 (constant D) additive hill term,
%%%alpha is time dependent, drug kill is 0 at time 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_treated';

%Import all data
A2780_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC50_1um_30k_cellsinwell.csv');
A2780_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC75_0.62um_30k_cellsinwell.csv');
A2780_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC25_1.47um_30k_cellsinwell.csv');

A2780_t_days = A2780_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
A2780_IC50_y_raw = [A2780_IC50_data.Rep1, A2780_IC50_data.Rep2, A2780_IC50_data.Rep3];
A2780_IC25_y_raw = [A2780_IC25_data.Rep1, A2780_IC25_data.Rep2, A2780_IC25_data.Rep3];
A2780_IC75_y_raw = [A2780_IC75_data.Rep1, A2780_IC75_data.Rep2, A2780_IC75_data.Rep3];

%downsample data
A2780_10t_days = A2780_t_days(1:11);
A2780_IC25_10y_raw = A2780_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
A2780_IC50_10y_raw = A2780_IC50_y_raw(1:11,:);
A2780_IC75_10y_raw = A2780_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 1.00; 
ce_IC25 = 0.62;
ce_IC75 = 1.47;

%Define logistic growth constants 
rS = 0.8547; %growth rate (day^-1)
KS = 1114404.9924;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 0.9960177; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [0.5];      % [alpha]
lb = [0.0001];   % reasonable bounds
ub = [20];

res_fun_nandalpha_constD_multlinear_alphaonly_t0 = @(p) all_datasets_residual_reps_n_and_alpha_constD_multlinear_alphaonly_t0(p, rS, KS, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

%run globalsearch with 200 starts
problem_nandalpha_constD_multlinear_alphaonly_t0 = createOptimProblem(...
    'lsqnonlin', ... %can't use lsqcurvefit because it expects x and y arguments
    'x0', p0, ...
    'objective', res_fun_nandalpha_constD_multlinear_alphaonly_t0, ...
    'lb', lb, ...
    'ub', ub);

ms_nandalpha_constD_multlinear_alphaonly_t0 = MultiStart('UseParallel', true, 'Display', 'iter', 'StartPointsToRun','bounds-ineqs');  % Set true if desired
nStarts = 200; %more start to ensure reproducibility 
[malphat_hat, resnorm_malphat] = run(ms_nandalpha_constD_multlinear_alphaonly_t0, problem_nandalpha_constD_multlinear_alphaonly_t0, nStarts);
fprintf('Fitted alpha only - add drug - alpha T t0=0 - constD A2780 (14 days)\n')
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted alpha only - add drug - alpha T t0=0 - constD A2780 (14 days)
% Fitted alpha term a = 0.1168
% Sum of squared residuals = 880169838113.89

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_nandalpha_constD_multlinear_alphaonly_t0, malphat_hat, lb, ub);

% Compute 95% confidence intervals
ci = nlparci(p_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95 Confidence Intervals:\n')
fprintf('a: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
% Local minimum possible.
% lsqnonlin stopped because the size of the current step is less than
% the value of the step size tolerance.
% <stopping criteria details>
% 95 Confidence Intervals:
% a: [0.1454, 0.1877]


[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_multlinear_alphaonly_t0(malphat_hat, rS, KS, ...
                                         A2780_10t_days, A2780_IC50_10y_raw, ce_IC50, ...
                                         A2780_10t_days, A2780_IC25_10y_raw, ce_IC25, ...
                                         A2780_10t_days, A2780_IC75_10y_raw, ce_IC75);

fprintf('Model (a=%.2f):\n', malphat_hat(1));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (a=0.12):
%   BIC IC50  = 720.54
%   BIC IC25  = 767.93
%   BIC IC75  = 770.58
%   BIC total = 2259.05

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {A2780_IC50_10y_raw, A2780_IC25_10y_raw, A2780_IC75_10y_raw};
model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(A2780_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_multlinear_alphaonly_t0(malphat_hat(1), A2780_10t_days, S0, ...
            ce_all(i), rS, KS);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(A2780_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 1400000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (a=%.2f):\n', malphat_hat(1)));

saveas(gcf, fullfile(outdir,'A2780_treated_alphaonly_multdruglinear_alphaT0_constD_10days_Model21.png'));

writematrix([A2780_10t_days, model_fit_all{1}], fullfile(outdir,'A2780_IC50treated_alphaonly_multdruglinear_alphaT0_constD_10days_Model21.csv'));
writematrix([A2780_10t_days, model_fit_all{2}], fullfile(outdir,'A2780_IC25treated_alphaonly_multdruglinear_alphaT0_constD_10days_Model21.csv'));
writematrix([A2780_10t_days, model_fit_all{3}], fullfile(outdir,'A2780_IC75treated_alphaonly_multdruglinear_alphaT0_constD_10days_Model21.csv'));