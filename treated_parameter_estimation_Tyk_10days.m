%Try different treatment models for Tyk for 10 days, the other cell lines
%are in other scripts 
%of data 
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_treated_10days';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Starts with Tyk treated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Import all data
Tyk_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC50_0.66um_40k_cellsinwell.csv');
Tyk_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC75_0.12um_40k_cellsinwell.csv');
Tyk_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC25_0.85um_40k_cellsinwell.csv');

Tyk_t_days = Tyk_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
Tyk_IC50_y_raw = [Tyk_IC50_data.Rep1, Tyk_IC50_data.Rep2, Tyk_IC50_data.Rep3];
Tyk_IC25_y_raw = [Tyk_IC25_data.Rep1, Tyk_IC25_data.Rep2, Tyk_IC25_data.Rep3];
Tyk_IC75_y_raw = [Tyk_IC75_data.Rep1, Tyk_IC75_data.Rep2, Tyk_IC75_data.Rep3];

Tyk_10t_days = Tyk_t_days(1:11);
Tyk_IC25_10y_raw = Tyk_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
Tyk_IC50_10y_raw = Tyk_IC50_y_raw(1:11,:);
Tyk_IC75_10y_raw = Tyk_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 0.66; 
ce_IC25 = 0.12;
ce_IC75 = 0.85;

%Define logistic growth constants 
rS = 0.6433; %growth rate (day^-1)
KS = 545681.7362;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 1.974574; %the IC50 for sensitive cells, after 72hr exposure

%MODEL 1 (Equation 3)
%%%Estimate for n and alpha - Tyk (constant D) multiplicative drug kill term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Estimate for parameters
p0 = [2; 0.5];      % [n; alpha]
lb = [0.01; 0.01];   % reasonable bounds
ub = [10; 2.0];

res_fun_nandalpha_constD = @(p) all_datasets_residual_reps_n_and_alpha_constD_2params(p, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

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
fprintf('Fitted n and alpha - mult - constD Tyk (10 days)\n')
fprintf('Fitted Hill coefficient n = %.4f\n', d_hat(1));
fprintf('Fitted alpha term a = %.4f\n', d_hat(2));
fprintf('Sum of squared residuals = %.2f\n', resnorm_d);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted n and alpha - mult - constD Tyk (10 days)
% Fitted Hill coefficient n = 1.5958
% Fitted alpha term a = 1.4426
% Sum of squared residuals = 1075311714348.74

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_nandalpha_constD, d_hat, lb, ub);

% Compute 95% confidence intervals
ci = nlparci(p_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95 Confidence Intervals:\n')
fprintf('n: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
fprintf('a: [%.4f, %.4f]\n', ci(2,1), ci(2,2));
% Local minimum possible.
% lsqnonlin stopped because the size of the current step is less than
% the value of the step size tolerance.
% <stopping criteria details>
% 95 Confidence Intervals:
% n: [0.7882, 2.4033]
% a: [0.5507, 2.3344]

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_2params(d_hat, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

fprintf('Model (n=%.2f, a=%.2f):\n', d_hat(1), d_hat(2));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (n=1.60, a=1.44):
%   BIC IC50  = 772.90
%   BIC IC25  = 753.50
%   BIC IC75  = 777.71
%   BIC total = 2304.11

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {Tyk_IC50_10y_raw, Tyk_IC25_10y_raw, Tyk_IC75_10y_raw};

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(Tyk_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstantD(d_hat(1), d_hat(2), Tyk_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(Tyk_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 600000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (n=%.2f, a=%.2f):\n', d_hat(1), d_hat(2)));

saveas(gcf, fullfile(outdir,'Tyk_treated_nandalpha_mult_constD_10days_model1.png'));

writematrix(S_raw_all{1,1}, fullfile(outdir, 'Tyk_IC50treated_nandalpha_mult_constD_10days_model1.csv'));
writematrix(S_raw_all{1,2}, fullfile(outdir, 'Tyk_IC25treated_nandalpha_mult_constD_10days_model1.csv'));
writematrix(S_raw_all{1,3}, fullfile(outdir, 'Tyk_IC75treated_nandalpha_mult_constD_10days_model1.csv'));

%MODEL 2 (Equation 4)
%%%Estimate for n and alpha - Tyk (constant D) additive drug kill term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_treated_10days';

%Import all data
Tyk_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC50_0.66um_40k_cellsinwell.csv');
Tyk_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC75_0.12um_40k_cellsinwell.csv');
Tyk_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC25_0.85um_40k_cellsinwell.csv');

Tyk_t_days = Tyk_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
Tyk_IC50_y_raw = [Tyk_IC50_data.Rep1, Tyk_IC50_data.Rep2, Tyk_IC50_data.Rep3];
Tyk_IC25_y_raw = [Tyk_IC25_data.Rep1, Tyk_IC25_data.Rep2, Tyk_IC25_data.Rep3];
Tyk_IC75_y_raw = [Tyk_IC75_data.Rep1, Tyk_IC75_data.Rep2, Tyk_IC75_data.Rep3];

%downsample data
Tyk_10t_days = Tyk_t_days(1:11);
Tyk_IC25_10y_raw = Tyk_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
Tyk_IC50_10y_raw = Tyk_IC50_y_raw(1:11,:);
Tyk_IC75_10y_raw = Tyk_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 0.66; 
ce_IC25 = 0.12;
ce_IC75 = 0.85;

%Define logistic growth constants 
rS = 0.6433; %growth rate (day^-1)
KS = 545681.7362;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 1.974574; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [2; 0.5];      % [n; alpha]
lb = [0.01; 0.01];   % reasonable bounds
ub = [10; 2.0];

res_fun_nandalpha_constD_add = @(p) all_datasets_residual_reps_n_and_alpha_constD_add(p, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

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
fprintf('Fitted n and alpha - add - constD Tyk (10 days)\n')
fprintf('Fitted Hill coefficient n = %.4f\n', a_hat(1));
fprintf('Fitted alpha term a = %.4f\n', a_hat(2));
fprintf('Sum of squared residuals = %.2f\n', resnorm_a);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted n and alpha - add - constD Tyk (10 days)
% Fitted Hill coefficient n = 1.9034
% Fitted alpha term a = 1.5825
% Sum of squared residuals = 823033880369.43

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_nandalpha_constD_add, a_hat, lb, ub);

% Compute 95% confidence intervals
ci = nlparci(p_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95 Confidence Intervals:\n')
fprintf('n: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
fprintf('a: [%.4f, %.4f]\n', ci(2,1), ci(2,2));
% Local minimum possible.
% lsqnonlin stopped because the size of the current step is less than
% the value of the step size tolerance.
% <stopping criteria details>
% 95 Confidence Intervals:
% n: [0.8508, 2.9560]
% a: [0.2318, 2.9332]

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_add(a_hat, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

fprintf('Model (n=%.2f, a=%.2f):\n', a_hat(1), a_hat(2));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (n=1.90, a=1.58):
%   BIC IC50  = 760.35
%   BIC IC25  = 751.37
%   BIC IC75  = 768.36
%   BIC total = 2280.08

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {Tyk_IC50_10y_raw, Tyk_IC25_10y_raw, Tyk_IC75_10y_raw};

model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(Tyk_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstantD_add(a_hat(1),a_hat(2), Tyk_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    model_fit_all{i} = S_fit_mean;
    plot(Tyk_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 600000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (n=%.2f, a=%.2f):\n', a_hat(1), a_hat(2)));

saveas(gcf, fullfile(outdir,'Tyk_treated_nandalpha_add_constD_10days_model2.png'));

writematrix([Tyk_10t_days, model_fit_all{1}], fullfile(outdir,'Tyk_IC50treated_nandalpha_add_constD_10days_model2.csv'));
writematrix([Tyk_10t_days, model_fit_all{2}], fullfile(outdir,'Tyk_IC25treated_nandalpha_add_constD_10days_model2.csv'));
writematrix([Tyk_10t_days, model_fit_all{3}], fullfile(outdir,'Tyk_IC75treated_nandalpha_add_constD_10days_model2.csv'));

% writematrix(S_raw_all{1,1}, fullfile(outdir, 'Tyk_IC50treated_nandalpha_add_constD_10days_model2.csv'));
% writematrix(S_raw_all{1,2}, fullfile(outdir, 'Tyk_IC25treated_nandalpha_add_constD_10days_model2.csv'));
% writematrix(S_raw_all{1,3}, fullfile(outdir, 'Tyk_IC75treated_nandalpha_add_constD_10days_model2.csv'));

%MODEL 3 (Equation 5)
%%%Estimate for n, alpha and d - Tyk (constant D) multiplicative drug
%%% term and additive death term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_treated_10days';

%Import all data
Tyk_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC50_0.66um_40k_cellsinwell.csv');
Tyk_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC75_0.12um_40k_cellsinwell.csv');
Tyk_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC25_0.85um_40k_cellsinwell.csv');

Tyk_t_days = Tyk_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
Tyk_IC50_y_raw = [Tyk_IC50_data.Rep1, Tyk_IC50_data.Rep2, Tyk_IC50_data.Rep3];
Tyk_IC25_y_raw = [Tyk_IC25_data.Rep1, Tyk_IC25_data.Rep2, Tyk_IC25_data.Rep3];
Tyk_IC75_y_raw = [Tyk_IC75_data.Rep1, Tyk_IC75_data.Rep2, Tyk_IC75_data.Rep3];

%downsample data
Tyk_10t_days = Tyk_t_days(1:11);
Tyk_IC25_10y_raw = Tyk_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
Tyk_IC50_10y_raw = Tyk_IC50_y_raw(1:11,:);
Tyk_IC75_10y_raw = Tyk_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 0.66; 
ce_IC25 = 0.12;
ce_IC75 = 0.85;

%Define logistic growth constants 
rS = 0.6433; %growth rate (day^-1)
KS = 545681.7362;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 1.974574; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [2; 0.5; 0.5];      % [n; alpha]
lb = [0.01; 0.01; 0.01];   % reasonable bounds
ub = [20; 2.0; 4.0];

res_fun_nandalpha_constD_death = @(p) all_datasets_residual_reps_n_and_alpha_constD_death(p, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

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
fprintf('Fitted n, alpha and d - mult drug - add death - constD Tyk (10 days)\n')
fprintf('Fitted Hill coefficient n = %.4f\n', mdeath_hat(1));
fprintf('Fitted alpha term a = %.4f\n', mdeath_hat(2));
fprintf('Fitted death term d = %.4f\n', mdeath_hat(3));
fprintf('Sum of squared residuals = %.2f\n', resnorm_mdeath);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted n, alpha and d - mult drug - add death - constD Tyk (10 days)
% Fitted Hill coefficient n = 1.6703
% Fitted alpha term a = 1.4548
% Fitted death term d = 0.0100
% Sum of squared residuals = 1084561987710.42

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_nandalpha_constD_death, mdeath_hat, lb, ub);

% Compute 95% confidence intervals
ci = nlparci(p_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95 Confidence Intervals:\n')
fprintf('n: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
fprintf('a: [%.4f, %.4f]\n', ci(2,1), ci(2,2));
fprintf('d: [%.4f, %.4f]\n', ci(3,1), ci(3,2));
% Local minimum possible.
% lsqnonlin stopped because the final change in the sum of squares relative to 
% its initial value is less than the value of the function tolerance.
% 95 Confidence Intervals:
% n: [0.5248, 2.8159]
% a: [0.3607, 2.5492]
% d: [-0.0451, 0.0651]

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_death(mdeath_hat, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

fprintf('Model (n=%.2f, a=%.2f, d=%.2f):\n', mdeath_hat(1), mdeath_hat(2), mdeath_hat(3));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (n=1.67, a=1.45, d=0.01):
%   BIC IC50  = 775.80
%   BIC IC25  = 760.07
%   BIC IC75  = 780.83
%   BIC total = 2316.70


%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {Tyk_IC50_10y_raw, Tyk_IC25_10y_raw, Tyk_IC75_10y_raw};
model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(Tyk_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstantD_death(mdeath_hat(1), mdeath_hat(2), mdeath_hat(3), Tyk_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(Tyk_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 600000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (n=%.2f, a=%.2f, d=%.2f):\n', mdeath_hat(1), mdeath_hat(2), mdeath_hat(3)));

saveas(gcf, fullfile(outdir,'Tyk_treated_n-alphaandd_multdrug_adddeath_constD_10days_model3.png'));

writematrix([Tyk_10t_days, model_fit_all{1}], fullfile(outdir,'Tyk_IC50treated_n-alphaandd_multdrug_adddeath_constD_10days_model3.csv'));
writematrix([Tyk_10t_days, model_fit_all{2}], fullfile(outdir,'Tyk_IC25treated_n-alphaandd_multdrug_adddeath_constD_10days_model3.csv'));
writematrix([Tyk_10t_days, model_fit_all{3}], fullfile(outdir,'Tyk_IC75treated_n-alphaandd_multdrug_adddeath_constD_10days_model3.csv'));

% writematrix(S_raw_all{1,1}, fullfile(outdir, 'Tyk_IC50treated_n-alphaandd_multdrug_adddeath_constD_10days_model3.csv'));
% writematrix(S_raw_all{1,2}, fullfile(outdir, 'Tyk_IC25treated_n-alphaandd_multdrug_adddeath_constD_10days_model3.csv'));
% writematrix(S_raw_all{1,3}, fullfile(outdir, 'Tyk_IC75treated_n-alphaandd_multdrug_adddeath_constD_10days_model3.csv'));

%MODEL 4 (Equation 6)
%%%Estimate for n, alpha and d - Tyk (constant D), additive drug kill term
%%%and death term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_treated_10days';

%Import all data
Tyk_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC50_0.66um_40k_cellsinwell.csv');
Tyk_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC75_0.12um_40k_cellsinwell.csv');
Tyk_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC25_0.85um_40k_cellsinwell.csv');

Tyk_t_days = Tyk_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
Tyk_IC50_y_raw = [Tyk_IC50_data.Rep1, Tyk_IC50_data.Rep2, Tyk_IC50_data.Rep3];
Tyk_IC25_y_raw = [Tyk_IC25_data.Rep1, Tyk_IC25_data.Rep2, Tyk_IC25_data.Rep3];
Tyk_IC75_y_raw = [Tyk_IC75_data.Rep1, Tyk_IC75_data.Rep2, Tyk_IC75_data.Rep3];

%downsample data
Tyk_10t_days = Tyk_t_days(1:11);
Tyk_IC25_10y_raw = Tyk_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
Tyk_IC50_10y_raw = Tyk_IC50_y_raw(1:11,:);
Tyk_IC75_10y_raw = Tyk_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 0.66; 
ce_IC25 = 0.12;
ce_IC75 = 0.85;

%Define logistic growth constants 
rS = 0.6433; %growth rate (day^-1)
KS = 545681.7362;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 1.974574; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [2; 0.5; 0.5];      % [n; alpha, d]
lb = [0.01; 0.01; 0.01];   % reasonable bounds
ub = [20; 2.0; 4.0];

res_fun_nandalpha_constD_add_death = @(p) all_datasets_residual_reps_n_and_alpha_constD_add_death(p, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

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
fprintf('Fitted n, alpha and d - add drug - add death - constD Tyk (10 days)\n')
fprintf('Fitted Hill coefficient n = %.4f\n', adeath_hat(1));
fprintf('Fitted alpha term a = %.4f\n', adeath_hat(2));
fprintf('Fitted death term d = %.4f\n', adeath_hat(3));
fprintf('Sum of squared residuals = %.2f\n', resnorm_adeath);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted n, alpha and d - add drug - add death - constD Tyk (10 days)
% Fitted Hill coefficient n = 1.9695
% Fitted alpha term a = 1.5925
% Fitted death term d = 0.0100
% Sum of squared residuals = 846487334259.08

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_nandalpha_constD_add_death, adeath_hat, lb, ub);

% Compute 95% confidence intervals
ci = nlparci(p_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95 Confidence Intervals:\n')
fprintf('n: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
fprintf('a: [%.4f, %.4f]\n', ci(2,1), ci(2,2));
fprintf('d: [%.4f, %.4f]\n', ci(3,1), ci(3,2));
% Local minimum possible.
% lsqnonlin stopped because the final change in the sum of squares relative to 
% its initial value is less than the value of the function tolerance.
% <stopping criteria details>
% 95 Confidence Intervals:
% n: [0.5575, 3.3820]
% a: [-0.0386, 3.2243]
% d: [-0.0408, 0.0608]

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_add_death(adeath_hat, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

fprintf('Model (n=%.2f, a=%.2f, d=%.2f):\n', adeath_hat(1), adeath_hat(2), adeath_hat(3));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (n=1.97, a=1.59, d=0.01):
%   BIC IC50  = 763.84
%   BIC IC25  = 758.46
%   BIC IC75  = 771.84
%   BIC total = 2294.14

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {Tyk_IC50_10y_raw, Tyk_IC25_10y_raw, Tyk_IC75_10y_raw};
model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(Tyk_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstantD_add_death(adeath_hat(1),adeath_hat(2),adeath_hat(3), Tyk_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(Tyk_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 600000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (n=%.2f, a=%.2f, d=%.2f):\n', adeath_hat(1), adeath_hat(2), adeath_hat(3)));

saveas(gcf, fullfile(outdir,'Tyk_treated_n-alphaandd_adddrug_adddeath_constD_10days_model4.png'));

writematrix([Tyk_10t_days, model_fit_all{1}], fullfile(outdir,'Tyk_IC50treated_n-alphaandd_adddrug_adddeath_constD_10days_model4.csv'));
writematrix([Tyk_10t_days, model_fit_all{2}], fullfile(outdir,'Tyk_IC25treated_n-alphaandd_adddrug_adddeath_constD_10days_model4.csv'));
writematrix([Tyk_10t_days, model_fit_all{3}], fullfile(outdir,'Tyk_IC75treated_n-alphaandd_adddrug_adddeath_constD_10days_model4.csv'));

% writematrix(S_raw_all{1,1}, fullfile(outdir, 'Tyk_IC50treated_n-alphaandd_adddrug_adddeath_constD_10days_model4.csv'));
% writematrix(S_raw_all{1,2}, fullfile(outdir, 'Tyk_IC25treated_n-alphaandd_adddrug_adddeath_constD_10days_model4.csv'));
% writematrix(S_raw_all{1,3}, fullfile(outdir, 'Tyk_IC75treated_n-alphaandd_adddrug_adddeath_constD_10days_model4.csv'));

%MODEL 7 (Equation 7)
%%%Estimate for n and alpha - Tyk (constant D) multiplicative hill term,
%%%alpha is time dependent, drug kill is 0 at time 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_treated_10days';

%Import all data
Tyk_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC50_0.66um_40k_cellsinwell.csv');
Tyk_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC75_0.12um_40k_cellsinwell.csv');
Tyk_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC25_0.85um_40k_cellsinwell.csv');

Tyk_t_days = Tyk_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
Tyk_IC50_y_raw = [Tyk_IC50_data.Rep1, Tyk_IC50_data.Rep2, Tyk_IC50_data.Rep3];
Tyk_IC25_y_raw = [Tyk_IC25_data.Rep1, Tyk_IC25_data.Rep2, Tyk_IC25_data.Rep3];
Tyk_IC75_y_raw = [Tyk_IC75_data.Rep1, Tyk_IC75_data.Rep2, Tyk_IC75_data.Rep3];

%downsample data
Tyk_10t_days = Tyk_t_days(1:11);
Tyk_IC25_10y_raw = Tyk_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
Tyk_IC50_10y_raw = Tyk_IC50_y_raw(1:11,:);
Tyk_IC75_10y_raw = Tyk_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 0.66; 
ce_IC25 = 0.12;
ce_IC75 = 0.85;

%Define logistic growth constants 
rS = 0.6433; %growth rate (day^-1)
KS = 545681.7362;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 1.974574; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [2; 0.5];      % [n; alpha]
lb = [0.01; 0.01];   % reasonable bounds
ub = [20; 20];

res_fun_nandalpha_constD_mult_alpha_t0 = @(p) all_datasets_residual_reps_n_and_alpha_constD_mult_alpha_t0(p, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

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
fprintf('Fitted n, and alpha - mult drug - alpha T t0=0 - constD Tyk (14 days)\n')
fprintf('Fitted Hill coefficient n = %.4f\n', malphat_hat(1));
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(2));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% 126 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted n, and alpha - mult drug - alpha T t0=0 - constD Tyk (14 days)
% Fitted Hill coefficient n = 1.2892
% Fitted alpha term a = 0.2939
% Sum of squared residuals = 523783466829.76

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_nandalpha_constD_mult_alpha_t0, malphat_hat, lb, ub);

% Compute 95% confidence intervals
ci = nlparci(p_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95 Confidence Intervals:\n')
fprintf('n: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
fprintf('a: [%.4f, %.4f]\n', ci(2,1), ci(2,2));
% Local minimum possible.
% lsqnonlin stopped because the final change in the sum of squares relative to 
% its initial value is less than the value of the function tolerance.
% <stopping criteria details>
% 95 Confidence Intervals:
% n: [0.7913, 1.7871]
% a: [0.1839, 0.4038]

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_mult_alpha_t0(malphat_hat, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

fprintf('Model (n=%.2f, a=%.2f):\n', malphat_hat(1), malphat_hat(2));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (n=1.29, a=0.29):
%   BIC IC50  = 738.84
%   BIC IC25  = 754.13
%   BIC IC75  = 742.56
%   BIC total = 2235.53

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {Tyk_IC50_10y_raw, Tyk_IC25_10y_raw, Tyk_IC75_10y_raw};
model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(Tyk_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_mult_alpha_t0(malphat_hat(1), malphat_hat(2), Tyk_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(Tyk_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 600000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (n=%.2f, a=%.2f):\n', malphat_hat(1), malphat_hat(2)));

saveas(gcf, fullfile(outdir,'Tyk_treated_nandalpha_multdrug_alphaT0_constD_10days_model7.png'));

writematrix([Tyk_10t_days, model_fit_all{1}], fullfile(outdir,'Tyk_IC50treated_nandalpha_multdrug_alphaT0_constD_10days_model7.csv'));
writematrix([Tyk_10t_days, model_fit_all{2}], fullfile(outdir,'Tyk_IC25treated_nandalpha_multdrug_alphaT0_constD_10days_model7.csv'));
writematrix([Tyk_10t_days, model_fit_all{3}], fullfile(outdir,'Tyk_IC75treated_nandalpha_multdrug_alphaT0_constD_10days_model7.csv'));

% writematrix(S_raw_all{1,1}, fullfile(outdir, 'Tyk_IC50treated_n-alphaandbeta_multdrug_alphaT0_constD_10days.csv'));
% writematrix(S_raw_all{1,2}, fullfile(outdir, 'Tyk_IC25treated_n-alphaandbeta_multdrug_alphaT0_constD_10days.csv'));
% writematrix(S_raw_all{1,3}, fullfile(outdir, 'Tyk_IC75treated_n-alphaandbeta_multdrug_alphaT0_constD_10days.csv'));


%MODEL 13 (Equation 8)
%%%Estimate for n and alpha - Tyk (constant D) additive hill term,
%%%alpha is time dependent, drug kill is 0 at time 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_treated_10days';

%Import all data
Tyk_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC50_0.66um_40k_cellsinwell.csv');
Tyk_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC75_0.12um_40k_cellsinwell.csv');
Tyk_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC25_0.85um_40k_cellsinwell.csv');

Tyk_t_days = Tyk_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
Tyk_IC50_y_raw = [Tyk_IC50_data.Rep1, Tyk_IC50_data.Rep2, Tyk_IC50_data.Rep3];
Tyk_IC25_y_raw = [Tyk_IC25_data.Rep1, Tyk_IC25_data.Rep2, Tyk_IC25_data.Rep3];
Tyk_IC75_y_raw = [Tyk_IC75_data.Rep1, Tyk_IC75_data.Rep2, Tyk_IC75_data.Rep3];

%downsample data
Tyk_10t_days = Tyk_t_days(1:11);
Tyk_IC25_10y_raw = Tyk_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
Tyk_IC50_10y_raw = Tyk_IC50_y_raw(1:11,:);
Tyk_IC75_10y_raw = Tyk_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 0.66; 
ce_IC25 = 0.12;
ce_IC75 = 0.85;

%Define logistic growth constants 
rS = 0.6433; %growth rate (day^-1)
KS = 545681.7362;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 1.974574; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [1; 1];      % [n; alpha] %initially 2 and 0.5 same results for those initial guesses
lb = [0.01; 0.01];   % reasonable bounds
ub = [20; 20];

res_fun_nandalpha_constD_add_alpha_t0 = @(p) all_datasets_residual_reps_n_and_alpha_constD_add_alpha_t0(p, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

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
fprintf('Fitted n, and alpha - add drug - alpha T t0=0 - constD Tyk (14 days)\n')
fprintf('Fitted Hill coefficient n = %.4f\n', malphat_hat(1));
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(2));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% 175 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted n, and alpha - add drug - alpha T t0=0 - constD Tyk (14 days)
% Fitted Hill coefficient n = 1.8145
% Fitted alpha term a = 0.3230
% Sum of squared residuals = 446995591076.50

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_nandalpha_constD_add_alpha_t0, malphat_hat, lb, ub);

% Compute 95% confidence intervals
ci = nlparci(p_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95 Confidence Intervals:\n')
fprintf('n: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
fprintf('a: [%.4f, %.4f]\n', ci(2,1), ci(2,2));
% Solver stopped prematurely.
% lsqnonlin stopped because it exceeded the function evaluation limit,
% options.MaxFunctionEvaluations = 2.000000e+02.
% 95 Confidence Intervals:
% n: [1.0607, 2.5683]
% a: [0.1233, 0.5227]

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_add_alpha_t0(malphat_hat, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

fprintf('Model (n=%.2f, a=%.2f):\n', malphat_hat(1), malphat_hat(2));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (n=1.81, a=0.32):
%   BIC IC50  = 730.21
%   BIC IC25  = 750.87
%   BIC IC75  = 737.33
%   BIC total = 2218.41

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {Tyk_IC50_10y_raw, Tyk_IC25_10y_raw, Tyk_IC75_10y_raw};

model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(Tyk_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_add_alpha_t0(malphat_hat(1), malphat_hat(2), Tyk_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    model_fit_all{i} = S_fit_mean;
    plot(Tyk_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 600000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (n=%.2f, a=%.2f):\n', malphat_hat(1), malphat_hat(2)));

saveas(gcf, fullfile(outdir,'Tyk_treated_nandalpha_adddrug_alphaT0_constD_10days_Model13.png'));

writematrix([Tyk_10t_days, model_fit_all{1}], fullfile(outdir,'Tyk_IC50treated_nandalpha_adddrug_alphaT0_constD_10days_Model13.csv'));
writematrix([Tyk_10t_days, model_fit_all{2}], fullfile(outdir,'Tyk_IC25treated_nandalpha_adddrug_alphaT0_constD_10days_Model13.csv'));
writematrix([Tyk_10t_days, model_fit_all{3}], fullfile(outdir,'Tyk_IC75treated_nandalpha_adddrug_alphaT0_constD_10days_Model13.csv'));

% writematrix(S_raw_all{1,1}, fullfile(outdir, 'Tyk_IC50treated_nandalpha_adddrug_alphaT0_constD_10days_Model13.csv'));
% writematrix(S_raw_all{1,2}, fullfile(outdir, 'Tyk_IC25treated_nandalpha_adddrug_alphaT0_constD_10days_Model13.csv'));
% writematrix(S_raw_all{1,3}, fullfile(outdir, 'Tyk_IC75treated_nandalpha_adddrug_alphaT0_constD_10days_Model13.csv'));

%MODEL 15 (Equation 9)
%%%Estimate for alpha only - Tyk (constant D) multiplicative hill term,
%%%alpha is time dependent, drug kill is 0 at time 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_treated_10days';

%Import all data
Tyk_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC50_0.66um_40k_cellsinwell.csv');
Tyk_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC75_0.12um_40k_cellsinwell.csv');
Tyk_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC25_0.85um_40k_cellsinwell.csv');

Tyk_t_days = Tyk_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
Tyk_IC50_y_raw = [Tyk_IC50_data.Rep1, Tyk_IC50_data.Rep2, Tyk_IC50_data.Rep3];
Tyk_IC25_y_raw = [Tyk_IC25_data.Rep1, Tyk_IC25_data.Rep2, Tyk_IC25_data.Rep3];
Tyk_IC75_y_raw = [Tyk_IC75_data.Rep1, Tyk_IC75_data.Rep2, Tyk_IC75_data.Rep3];

%downsample data
Tyk_10t_days = Tyk_t_days(1:11);
Tyk_IC25_10y_raw = Tyk_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
Tyk_IC50_10y_raw = Tyk_IC50_y_raw(1:11,:);
Tyk_IC75_10y_raw = Tyk_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 0.66; 
ce_IC25 = 0.12;
ce_IC75 = 0.85;

%Define logistic growth constants 
rS = 0.6433; %growth rate (day^-1)
KS = 545681.7362;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 1.974574; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [0.5];      % [alpha]
lb = [0.0001];   % reasonable bounds
ub = [20];

res_fun_nandalpha_constD_mult_alphaonly_t0 = @(p) all_datasets_residual_reps_n_and_alpha_constD_mult_alphaonly_t0(p, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

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
fprintf('Fitted alpha only - mult drug - alpha T t0=0 - constD Tyk (10 days)\n')
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted alpha only - mult drug - alpha T t0=0 - constD Tyk (10 days)
% Fitted alpha term a = 0.4766
% Sum of squared residuals = 541825699651.46

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_nandalpha_constD_mult_alphaonly_t0, malphat_hat, lb, ub);

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
% a: [0.4440, 0.5091]

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_mult_alphaonly_t0(malphat_hat, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

fprintf('Model (a=%.2f):\n', malphat_hat(1));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (a=0.48):
%   BIC IC50  = 732.83
%   BIC IC25  = 755.48
%   BIC IC75  = 737.54
%   BIC total = 2225.85

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {Tyk_IC50_10y_raw, Tyk_IC25_10y_raw, Tyk_IC75_10y_raw};
model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(Tyk_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_mult_alphaonly_t0(malphat_hat(1), Tyk_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(Tyk_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 600000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (a=%.2f):\n', malphat_hat(1)));

saveas(gcf, fullfile(outdir,'Tyk_treated_alphaonly_multdrug_alphaT0_constD_10days_Model15.png'));

writematrix([Tyk_10t_days, model_fit_all{1}], fullfile(outdir,'Tyk_IC50treated_alphaonly_multdrug_alphaT0_constD_10days_Model15.csv'));
writematrix([Tyk_10t_days, model_fit_all{2}], fullfile(outdir,'Tyk_IC25treated_alphaonly_multdrug_alphaT0_constD_10days_Model15.csv'));
writematrix([Tyk_10t_days, model_fit_all{3}], fullfile(outdir,'Tyk_IC75treated_alphaonly_multdrug_alphaT0_constD_10days_Model15.csv'));

% writematrix(S_raw_all{1,1}, fullfile(outdir, 'Tyk_IC50treated_alphaonly_multdrug_alphaT0_constD_10days_Model15.csv'));
% writematrix(S_raw_all{1,2}, fullfile(outdir, 'Tyk_IC25treated_alphaonly_multdrug_alphaT0_constD_10days_Model15.csv'));
% writematrix(S_raw_all{1,3}, fullfile(outdir, 'Tyk_IC75treated_alphaonly_multdrug_alphaT0_constD_10days_Model15.csv'));


%MODEL 14 (Equation 10)
%%%Estimate for alpha only - Tyk (constant D) additive hill term,
%%%alpha is time dependent, drug kill is 0 at time 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_treated_10days';

%Import all data
Tyk_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC50_0.66um_40k_cellsinwell.csv');
Tyk_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC75_0.12um_40k_cellsinwell.csv');
Tyk_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC25_0.85um_40k_cellsinwell.csv');

Tyk_t_days = Tyk_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
Tyk_IC50_y_raw = [Tyk_IC50_data.Rep1, Tyk_IC50_data.Rep2, Tyk_IC50_data.Rep3];
Tyk_IC25_y_raw = [Tyk_IC25_data.Rep1, Tyk_IC25_data.Rep2, Tyk_IC25_data.Rep3];
Tyk_IC75_y_raw = [Tyk_IC75_data.Rep1, Tyk_IC75_data.Rep2, Tyk_IC75_data.Rep3];

%downsample data
Tyk_10t_days = Tyk_t_days(1:11);
Tyk_IC25_10y_raw = Tyk_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
Tyk_IC50_10y_raw = Tyk_IC50_y_raw(1:11,:);
Tyk_IC75_10y_raw = Tyk_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 0.66; 
ce_IC25 = 0.12;
ce_IC75 = 0.85;

%Define logistic growth constants 
rS = 0.6433; %growth rate (day^-1)
KS = 545681.7362;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 1.974574; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [0.5];      % [alpha]
lb = [0.0001];   % reasonable bounds
ub = [20];

res_fun_nandalpha_constD_add_alphaonly_t0 = @(p) all_datasets_residual_reps_n_and_alpha_constD_add_alphaonly_t0(p, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

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
fprintf('Fitted alpha only - add drug - alpha T t0=0 - constD Tyk (14 days)\n')
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted alpha only - add drug - alpha T t0=0 - constD Tyk (14 days)
% Fitted alpha term a = 0.1666
% Sum of squared residuals = 533839592190.73 5e11

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_nandalpha_constD_add_alphaonly_t0, malphat_hat, lb, ub);

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


[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_add_alphaonly_t0(malphat_hat, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

fprintf('Model (a=%.2f):\n', malphat_hat(1));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (a=0.17):
%   BIC IC50  = 724.54
%   BIC IC25  = 760.01
%   BIC IC75  = 733.56
%   BIC total = 2218.11

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {Tyk_IC50_10y_raw, Tyk_IC25_10y_raw, Tyk_IC75_10y_raw};
model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(Tyk_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_add_alphaonly_t0(malphat_hat(1), Tyk_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    model_fit_all{i} = S_fit_mean;
    plot(Tyk_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 600000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (a=%.2f):\n', malphat_hat(1)));

saveas(gcf, fullfile(outdir,'Tyk_treated_alphaonly_adddrug_alphaT0_constD_10days_Model14.png'));

writematrix([Tyk_10t_days, model_fit_all{1}], fullfile(outdir,'Tyk_IC50treated_alphaonly_adddrug_alphaT0_constD_10days_Model14.csv'));
writematrix([Tyk_10t_days, model_fit_all{2}], fullfile(outdir,'Tyk_IC25treated_alphaonly_adddrug_alphaT0_constD_10days_Model14.csv'));
writematrix([Tyk_10t_days, model_fit_all{3}], fullfile(outdir,'Tyk_IC75treated_alphaonly_adddrug_alphaT0_constD_10days_Model14.csv'));

% writematrix(S_raw_all{1,1}, fullfile(outdir, 'Tyk_IC50treated_alphaonly_adddrug_alphaT0_constD_10days_Model14.csv'));
% writematrix(S_raw_all{1,2}, fullfile(outdir, 'Tyk_IC25treated_alphaonly_adddrug_alphaT0_constD_10days_Model14.csv'));
% writematrix(S_raw_all{1,3}, fullfile(outdir, 'Tyk_IC75treated_alphaonly_adddrug_alphaT0_constD_10days_Model14.csv'));

%R squared calculation
clear all;
Tyk_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC50_0.66um_40k_cellsinwell.csv');
Tyk_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC75_0.12um_40k_cellsinwell.csv');
Tyk_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC25_0.85um_40k_cellsinwell.csv');

Tyk_t_days = Tyk_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
Tyk_IC50_y_raw = [Tyk_IC50_data.Rep1, Tyk_IC50_data.Rep2, Tyk_IC50_data.Rep3];
Tyk_IC25_y_raw = [Tyk_IC25_data.Rep1, Tyk_IC25_data.Rep2, Tyk_IC25_data.Rep3];
Tyk_IC75_y_raw = [Tyk_IC75_data.Rep1, Tyk_IC75_data.Rep2, Tyk_IC75_data.Rep3];

Tyk_IC25_10y_raw = Tyk_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
Tyk_IC50_10y_raw = Tyk_IC50_y_raw(1:11,:);
Tyk_IC75_10y_raw = Tyk_IC75_y_raw(1:11,:);

Tyk_IC25_10y_all = Tyk_IC25_10y_raw(:);
Tyk_IC50_10y_all = Tyk_IC50_10y_raw(:);
Tyk_IC75_10y_all = Tyk_IC75_10y_raw(:);

%IC25
Tyk_IC25_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_treated_10days/Tyk_IC25treated_alphaonly_adddrug_alphaT0_constD_10days_Model14.csv');
Tyk_IC25_predicted = table2array(Tyk_IC25_predicted(:,2)); %only the second column 
Tyk_IC25_predicted = repmat(Tyk_IC25_predicted,3,1); %repeat 3 times 

SS_res_IC25 = sum((Tyk_IC25_10y_all - Tyk_IC25_predicted).^2);
SS_tot_IC25 = sum((Tyk_IC25_10y_all - mean(Tyk_IC25_10y_all)).^2);
R2_IC25 = 1 - SS_res_IC25 / SS_tot_IC25;

fprintf('R-squared for IC25: %.4f\n', R2_IC25);
%R-squared for IC25: 0.6670

%IC50
Tyk_IC50_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_treated_10days/Tyk_IC50treated_alphaonly_adddrug_alphaT0_constD_10days_Model14.csv');
Tyk_IC50_predicted = table2array(Tyk_IC50_predicted(:,2)); %only the second column 
Tyk_IC50_predicted = repmat(Tyk_IC50_predicted,3,1); %repeat 3 times 

SS_res_IC50 = sum((Tyk_IC50_10y_all - Tyk_IC50_predicted).^2);
SS_tot_IC50 = sum((Tyk_IC50_10y_all - mean(Tyk_IC50_10y_all)).^2);
R2_IC50 = 1 - SS_res_IC50 / SS_tot_IC50;

fprintf('R-squared for IC50: %.4f\n', R2_IC50);
%R-squared for IC50: 0.3850

%IC75
Tyk_IC75_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_treated_10days/Tyk_IC75treated_alphaonly_adddrug_alphaT0_constD_10days_Model14.csv');
Tyk_IC75_predicted = table2array(Tyk_IC75_predicted(:,2)); %only the second column 
Tyk_IC75_predicted = repmat(Tyk_IC75_predicted,3,1); %repeat 3 times 

SS_res_IC75 = sum((Tyk_IC75_10y_all - Tyk_IC75_predicted).^2);
SS_tot_IC75 = sum((Tyk_IC75_10y_all - mean(Tyk_IC75_10y_all)).^2);
R2_IC75 = 1 - SS_res_IC75 / SS_tot_IC75;

fprintf('R-squared for IC75: %.4f\n', R2_IC75);
%R-squared for IC75: 0.3259

%average of the 3 R2 squared average: mean([0.667,0.385,0.3259]) = 0.4593


%MODEL 8 (Equation 11)
%%%Estimate for n and alpha - Tyk (constant D) multiplicative hill term,
%%%alpha is time dependent, drug kill is alpha at t0 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_treated_10days';

%Import all data
Tyk_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC50_0.66um_40k_cellsinwell.csv');
Tyk_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC75_0.12um_40k_cellsinwell.csv');
Tyk_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC25_0.85um_40k_cellsinwell.csv');

Tyk_t_days = Tyk_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
Tyk_IC50_y_raw = [Tyk_IC50_data.Rep1, Tyk_IC50_data.Rep2, Tyk_IC50_data.Rep3];
Tyk_IC25_y_raw = [Tyk_IC25_data.Rep1, Tyk_IC25_data.Rep2, Tyk_IC25_data.Rep3];
Tyk_IC75_y_raw = [Tyk_IC75_data.Rep1, Tyk_IC75_data.Rep2, Tyk_IC75_data.Rep3];

%downsample data
Tyk_10t_days = Tyk_t_days(1:11);
Tyk_IC25_10y_raw = Tyk_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
Tyk_IC50_10y_raw = Tyk_IC50_y_raw(1:11,:);
Tyk_IC75_10y_raw = Tyk_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 0.66; 
ce_IC25 = 0.12;
ce_IC75 = 0.85;

%Define logistic growth constants 
rS = 0.6433; %growth rate (day^-1)
KS = 545681.7362;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 1.974574; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [2; 0.5];      % [n; alpha]
lb = [0.01; 0.01];   % reasonable bounds
ub = [20; 20];

res_fun_nandalpha_constD_mult_alpha_talpha = @(p) all_datasets_residual_reps_n_and_alpha_constD_mult_alpha_talpha(p, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

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
fprintf('Fitted n, and alpha - mult drug - alpha T t0=alpha - constD Tyk (14 days)\n')
fprintf('Fitted Hill coefficient n = %.4f\n', malphat_hat(1));
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(2));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% 25 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted n, and alpha - mult drug - alpha T t0=alpha - constD Tyk (14 days)
% Fitted Hill coefficient n = 1.3303
% Fitted alpha term a = 0.2422
% Sum of squared residuals = 608903530845.97

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_nandalpha_constD_mult_alpha_talpha, malphat_hat, lb, ub);

% Compute 95% confidence intervals
ci = nlparci(p_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95 Confidence Intervals:\n')
fprintf('n: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
fprintf('a: [%.4f, %.4f]\n', ci(2,1), ci(2,2));
% Solver stopped prematurely.
% lsqnonlin stopped because it exceeded the function evaluation limit,
% options.MaxFunctionEvaluations = 2.000000e+02.
% 95 Confidence Intervals:
% n: [0.7551, 1.9055]
% a: [0.1395, 0.3450]

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_mult_alpha_talpha(malphat_hat, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

fprintf('Model (n=%.2f, a=%.2f):\n', malphat_hat(1), malphat_hat(2));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (n=1.33, a=0.24):
%   BIC IC50  = 743.45
%   BIC IC25  = 750.90
%   BIC IC75  = 747.16
%   BIC total = 2241.51

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {Tyk_IC50_10y_raw, Tyk_IC25_10y_raw, Tyk_IC75_10y_raw};
model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(Tyk_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_mult_alpha_talpha(malphat_hat(1), malphat_hat(2), Tyk_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(Tyk_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 600000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (n=%.2f, a=%.2f):\n', malphat_hat(1), malphat_hat(2)));

saveas(gcf, fullfile(outdir,'Tyk_treated_nandalpha_multdrug_alphaTalpha_constD_10days_model8.png'));

writematrix([Tyk_10t_days, model_fit_all{1}], fullfile(outdir,'Tyk_IC50treated_nandalpha_multdrug_alphaTalpha_constD_10days_model8.csv'));
writematrix([Tyk_10t_days, model_fit_all{2}], fullfile(outdir,'Tyk_IC25treated_nandalpha_multdrug_alphaTalpha_constD_10days_model8.csv'));
writematrix([Tyk_10t_days, model_fit_all{3}], fullfile(outdir,'Tyk_IC75treated_nandalpha_multdrug_alphaTalpha_constD_10days_model8.csv'));

% writematrix(S_raw_all{1,1}, fullfile(outdir, 'Tyk_IC50treated_n-alphaandbeta_multdrug_alphaTalpha_constD_10days.csv'));
% writematrix(S_raw_all{1,2}, fullfile(outdir, 'Tyk_IC25treated_n-alphaandbeta_multdrug_alphaTalpha_constD_10days.csv'));
% writematrix(S_raw_all{1,3}, fullfile(outdir, 'Tyk_IC75treated_n-alphaandbeta_multdrug_alphaTalpha_constD_10days.csv'));


%MODEL 10 (Equation 12)
%%%Estimate for n, and alpha - Tyk (constant D) additive hill term,
%%%alpha is time dependent, t0=alpha 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_treated_10days';

%Import all data
Tyk_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC50_0.66um_40k_cellsinwell.csv');
Tyk_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC75_0.12um_40k_cellsinwell.csv');
Tyk_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC25_0.85um_40k_cellsinwell.csv');

Tyk_t_days = Tyk_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
Tyk_IC50_y_raw = [Tyk_IC50_data.Rep1, Tyk_IC50_data.Rep2, Tyk_IC50_data.Rep3];
Tyk_IC25_y_raw = [Tyk_IC25_data.Rep1, Tyk_IC25_data.Rep2, Tyk_IC25_data.Rep3];
Tyk_IC75_y_raw = [Tyk_IC75_data.Rep1, Tyk_IC75_data.Rep2, Tyk_IC75_data.Rep3];

%downsample data
Tyk_10t_days = Tyk_t_days(1:11);
Tyk_IC25_10y_raw = Tyk_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
Tyk_IC50_10y_raw = Tyk_IC50_y_raw(1:11,:);
Tyk_IC75_10y_raw = Tyk_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 0.66; 
ce_IC25 = 0.12;
ce_IC75 = 0.85;

%Define logistic growth constants 
rS = 0.6433; %growth rate (day^-1)
KS = 545681.7362;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 1.974574; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [3.0; 0.1];      % [n; alpha; beta]
lb = [0.01; 0.01];   % reasonable bounds
ub = [20.0; 20.0]; %change ub for the last one 

res_fun_nandalpha_constD_add_alpha_talpha = @(p) all_datasets_residual_reps_n_and_alpha_constD_add_alpha_talpha(p, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

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
fprintf('Fitted n and alpha - add drug - alpha T - t0=alpha - constD Tyk (10 days)\n')
fprintf('Fitted Hill coefficient n = %.4f\n', aalphat_hat(1));
fprintf('Fitted alpha term a = %.4f\n', aalphat_hat(2));
fprintf('Sum of squared residuals = %.2f\n', resnorm_aalphat);
% 4 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted n and alpha - add drug - alpha T - t0=alpha - constD Tyk (10 days)
% Fitted Hill coefficient n = 2.5255
% Fitted alpha term a = 0.4815
% Sum of squared residuals = 517604599262.18


% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_nandalpha_constD_add_alpha_talpha, aalphat_hat, lb, ub);

% Compute 95% confidence intervals
ci = nlparci(p_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95 Confidence Intervals:\n')
fprintf('n: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
fprintf('a: [%.4f, %.4f]\n', ci(2,1), ci(2,2));
% Solver stopped prematurely.
% lsqnonlin stopped because it exceeded the function evaluation limit,
% options.MaxFunctionEvaluations = 2.000000e+02.
% 95 Confidence Intervals:
% n: [2.2551, 2.4185]
% a: [0.3768, 0.4331]

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_add_alpha_talpha(aalphat_hat, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

fprintf('Model (n=%.2f, a=%.2f):\n', aalphat_hat(1), aalphat_hat(2));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (n=2.53, a=0.48):
%   BIC IC50  = 740.99
%   BIC IC25  = 748.95
%   BIC IC75  = 745.91
%   BIC total = 2235.85

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {Tyk_IC50_10y_raw, Tyk_IC25_10y_raw, Tyk_IC75_10y_raw};

model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(Tyk_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_add_alpha_talpha(aalphat_hat(1), aalphat_hat(2), Tyk_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    model_fit_all{i} = S_fit_mean;
    plot(Tyk_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 600000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (n=%.2f, a=%.2f):\n', aalphat_hat(1), aalphat_hat(2)));

saveas(gcf, fullfile(outdir,'Tyk_treated_nandalpha_adddrug_alphaT_t0alpha_constD_10days_model10.png'));

writematrix([Tyk_10t_days, model_fit_all{1}], fullfile(outdir,'Tyk_IC50treated_nandalpha_adddrug_alphaT_t0alpha_constD_10days_model10.csv'));
writematrix([Tyk_10t_days, model_fit_all{2}], fullfile(outdir,'Tyk_IC25treated_nandalpha_adddrug_alphaT_t0alpha_constD_10days_model10.csv'));
writematrix([Tyk_10t_days, model_fit_all{3}], fullfile(outdir,'Tyk_IC75treated_nandalpha_adddrug_alphaT_t0alpha_constD_10days_model10.csv'));

% writematrix(S_raw_all{1,1}, fullfile(outdir, 'Tyk_IC50treated_nandalpha_adddrug_alphaT_t0alpha_constD_10days.csv'));
% writematrix(S_raw_all{1,2}, fullfile(outdir, 'Tyk_IC25treated_nandalpha_adddrug_alphaT_t0alpha_constD_10days.csv'));
% writematrix(S_raw_all{1,3}, fullfile(outdir, 'Tyk_IC75treated_nandalpha_adddrug_alphaT_t0alpha_constD_10days.csv'));

%MODEL 17 (Equation 13)
%%%Estimate for alpha only - Tyk (constant D) multiplicative hill term,
%%%alpha is time dependent, drug kill is alpha at time 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_treated_10days';

%Import all data
Tyk_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC50_0.66um_40k_cellsinwell.csv');
Tyk_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC75_0.12um_40k_cellsinwell.csv');
Tyk_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC25_0.85um_40k_cellsinwell.csv');

Tyk_t_days = Tyk_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
Tyk_IC50_y_raw = [Tyk_IC50_data.Rep1, Tyk_IC50_data.Rep2, Tyk_IC50_data.Rep3];
Tyk_IC25_y_raw = [Tyk_IC25_data.Rep1, Tyk_IC25_data.Rep2, Tyk_IC25_data.Rep3];
Tyk_IC75_y_raw = [Tyk_IC75_data.Rep1, Tyk_IC75_data.Rep2, Tyk_IC75_data.Rep3];

%downsample data
Tyk_10t_days = Tyk_t_days(1:11);
Tyk_IC25_10y_raw = Tyk_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
Tyk_IC50_10y_raw = Tyk_IC50_y_raw(1:11,:);
Tyk_IC75_10y_raw = Tyk_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 0.66; 
ce_IC25 = 0.12;
ce_IC75 = 0.85;

%Define logistic growth constants 
rS = 0.6433; %growth rate (day^-1)
KS = 545681.7362;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 1.974574; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [0.5];      % [alpha]
lb = [0.0001];   % reasonable bounds
ub = [20];

res_fun_nandalpha_constD_mult_alphaonly_talpha = @(p) all_datasets_residual_reps_n_and_alpha_constD_mult_alphaonly_talpha(p, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

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
fprintf('Fitted alpha only - mult drug - alpha T t0=alpha - constD Tyk (10 days)\n')
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted alpha only - mult drug - alpha T t0=alpha - constD Tyk (10 days)
% Fitted alpha term a = 0.3816
% Sum of squared residuals = 633145375636.36

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_nandalpha_constD_mult_alphaonly_talpha, malphat_hat, lb, ub);

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
% a: [0.3530, 0.4102]


[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_mult_alphaonly_talpha(malphat_hat, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

fprintf('Model (a=%.2f):\n', malphat_hat(1));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (a=0.38):
%   BIC IC50  = 741.04
%   BIC IC25  = 756.84
%   BIC IC75  = 745.89
%   BIC total = 2243.77

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {Tyk_IC50_10y_raw, Tyk_IC25_10y_raw, Tyk_IC75_10y_raw};
model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(Tyk_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_mult_alphaonly_talpha(malphat_hat(1), Tyk_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(Tyk_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 600000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (a=%.2f):\n', malphat_hat(1)));

saveas(gcf, fullfile(outdir,'Tyk_treated_alphaonly_multdrug_alphaTalpha_constD_10days_Model17.png'));

writematrix([Tyk_10t_days, model_fit_all{1}], fullfile(outdir,'Tyk_IC50treated_alphaonly_multdrug_alphaTalpha_constD_10days_Model17.csv'));
writematrix([Tyk_10t_days, model_fit_all{2}], fullfile(outdir,'Tyk_IC25treated_alphaonly_multdrug_alphaTalpha_constD_10days_Model17.csv'));
writematrix([Tyk_10t_days, model_fit_all{3}], fullfile(outdir,'Tyk_IC75treated_alphaonly_multdrug_alphaTalpha_constD_10days_Model17.csv'));

% writematrix(S_raw_all{1,1}, fullfile(outdir, 'Tyk_IC50treated_alphaonly_multdrug_alphaTalpha_constD_10days_Model17.csv'));
% writematrix(S_raw_all{1,2}, fullfile(outdir, 'Tyk_IC25treated_alphaonly_multdrug_alphaTalpha_constD_10days_Model17.csv'));
% writematrix(S_raw_all{1,3}, fullfile(outdir, 'Tyk_IC75treated_alphaonly_multdrug_alphaTalpha_constD_10days_Model17.csv'));
%this one looks a little better for IC50 

%MODEL 16 (Equation 14)
%%%Estimate for alpha only - Tyk (constant D) additive hill term,
%%%alpha is time dependent, drug kill is alpha at time 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_treated_10days';

%Import all data
Tyk_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC50_0.66um_40k_cellsinwell.csv');
Tyk_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC75_0.12um_40k_cellsinwell.csv');
Tyk_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC25_0.85um_40k_cellsinwell.csv');

Tyk_t_days = Tyk_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
Tyk_IC50_y_raw = [Tyk_IC50_data.Rep1, Tyk_IC50_data.Rep2, Tyk_IC50_data.Rep3];
Tyk_IC25_y_raw = [Tyk_IC25_data.Rep1, Tyk_IC25_data.Rep2, Tyk_IC25_data.Rep3];
Tyk_IC75_y_raw = [Tyk_IC75_data.Rep1, Tyk_IC75_data.Rep2, Tyk_IC75_data.Rep3];

%downsample data
Tyk_10t_days = Tyk_t_days(1:11);
Tyk_IC25_10y_raw = Tyk_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
Tyk_IC50_10y_raw = Tyk_IC50_y_raw(1:11,:);
Tyk_IC75_10y_raw = Tyk_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 0.66; 
ce_IC25 = 0.12;
ce_IC75 = 0.85;

%Define logistic growth constants 
rS = 0.6433; %growth rate (day^-1)
KS = 545681.7362;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 1.974574; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [0.5];      % [alpha]
lb = [0.0001];   % reasonable bounds
ub = [20];

res_fun_nandalpha_constD_add_alphaonly_talpha = @(p) all_datasets_residual_reps_n_and_alpha_constD_add_alphaonly_talpha(p, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

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
fprintf('Fitted alpha only - add drug - alpha T t0=alpha - constD Tyk (10 days)\n')
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted alpha only - add drug - alpha T t0=alpha - constD Tyk (10 days)
% Fitted alpha term a = 0.1375
% Sum of squared residuals = 582020101806.49

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_nandalpha_constD_add_alphaonly_talpha, malphat_hat, lb, ub);

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
% a: [0.1236, 0.1514]

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_add_alphaonly_talpha(malphat_hat, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

fprintf('Model (a=%.2f):\n', malphat_hat(1));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (a=0.14):
%   BIC IC50  = 730.12
%   BIC IC25  = 760.49
%   BIC IC75  = 739.22
%   BIC total = 2229.83

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {Tyk_IC50_10y_raw, Tyk_IC25_10y_raw, Tyk_IC75_10y_raw};
model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(Tyk_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_add_alphaonly_talpha(malphat_hat(1), Tyk_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(Tyk_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 600000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (a=%.2f):\n', malphat_hat(1)));

saveas(gcf, fullfile(outdir,'Tyk_treated_alphaonly_adddrug_alphaTalpha_constD_10days_Model16.png'));

writematrix([Tyk_10t_days, model_fit_all{1}], fullfile(outdir,'Tyk_IC50treated_alphaonly_adddrug_alphaTalpha_constD_10days_Model16.csv'));
writematrix([Tyk_10t_days, model_fit_all{2}], fullfile(outdir,'Tyk_IC25treated_alphaonly_adddrug_alphaTalpha_constD_10days_Model16.csv'));
writematrix([Tyk_10t_days, model_fit_all{3}], fullfile(outdir,'Tyk_IC75treated_alphaonly_adddrug_alphaTalpha_constD_10days_Model16.csv'));

% writematrix(S_raw_all{1,1}, fullfile(outdir, 'Tyk_IC50treated_alphaonly_adddrug_alphaTalpha_constD_10days_Model16.csv'));
% writematrix(S_raw_all{1,2}, fullfile(outdir, 'Tyk_IC25treated_alphaonly_adddrug_alphaTalpha_constD_10days_Model16.csv'));
% writematrix(S_raw_all{1,3}, fullfile(outdir, 'Tyk_IC75treated_alphaonly_adddrug_alphaTalpha_constD_10days_Model16.csv'));


%MODEL 5 (Equation 15)
%%%Estimate for n, alpha and beta - Tyk (constant D) multiplicative hill term,
%%%alpha is time dependent 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_treated_10days';

%Import all data
Tyk_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC50_0.66um_40k_cellsinwell.csv');
Tyk_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC75_0.12um_40k_cellsinwell.csv');
Tyk_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC25_0.85um_40k_cellsinwell.csv');

Tyk_t_days = Tyk_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
Tyk_IC50_y_raw = [Tyk_IC50_data.Rep1, Tyk_IC50_data.Rep2, Tyk_IC50_data.Rep3];
Tyk_IC25_y_raw = [Tyk_IC25_data.Rep1, Tyk_IC25_data.Rep2, Tyk_IC25_data.Rep3];
Tyk_IC75_y_raw = [Tyk_IC75_data.Rep1, Tyk_IC75_data.Rep2, Tyk_IC75_data.Rep3];

%downsample data
Tyk_10t_days = Tyk_t_days(1:11);
Tyk_IC25_10y_raw = Tyk_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
Tyk_IC50_10y_raw = Tyk_IC50_y_raw(1:11,:);
Tyk_IC75_10y_raw = Tyk_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 0.66; 
ce_IC25 = 0.12;
ce_IC75 = 0.85;

%Define logistic growth constants 
rS = 0.6433; %growth rate (day^-1)
KS = 545681.7362;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 1.974574; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [2; 0.5; 5];      % [n; alpha; beta]
lb = [0.01; 0.01; 0.01];   % reasonable bounds
ub = [20; 20.0; 20.0];

res_fun_nandalpha_constD_mult_alpha_t = @(p) all_datasets_residual_reps_n_and_alpha_constD_mult_alpha_t(p, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

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
fprintf('Fitted n, alpha and beta - mult drug - alpha T - constD Tyk (10 days)\n')
fprintf('Fitted Hill coefficient n = %.4f\n', malphat_hat(1));
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(2));
fprintf('Fitted beta term b = %.4f\n', malphat_hat(3));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% 8 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted n, alpha and beta - mult drug - alpha T - constD Tyk (10 days)
% Fitted Hill coefficient n = 7.2450
% Fitted alpha term a = 7.6892
% Fitted beta term b = 3.9512
% Sum of squared residuals = 1181384494063.95

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_nandalpha_constD_mult_alpha_t, malphat_hat, lb, ub);

% Compute 95% confidence intervals
ci = nlparci(p_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95 Confidence Intervals:\n')
fprintf('n: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
fprintf('a: [%.4f, %.4f]\n', ci(2,1), ci(2,2));
fprintf('b: [%.4f, %.4f]\n', ci(3,1), ci(3,2));
% Local minimum possible.
% lsqnonlin stopped because the size of the current step is less than
% the value of the step size tolerance.
% 
% <stopping criteria details>
% 95 Confidence Intervals:
% n: [-0.1792, 14.6692]
% a: [1.1188, 14.2596]
% b: [-25.5564, 33.4588]

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_mult_alpha_t(malphat_hat, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

fprintf('Model (n=%.2f, a=%.2f, b=%.2f):\n', malphat_hat(1), malphat_hat(2), malphat_hat(3));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (n=7.25, a=7.69, b=3.95):
%   BIC IC50  = 800.79
%   BIC IC25  = 751.80
%   BIC IC75  = 747.07
%   BIC total = 2299.66

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {Tyk_IC50_10y_raw, Tyk_IC25_10y_raw, Tyk_IC75_10y_raw};
model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(Tyk_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_mult_alpha_t(malphat_hat(1), malphat_hat(2), malphat_hat(3), Tyk_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(Tyk_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 600000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (n=%.2f, a=%.2f, b=%.2f):\n', malphat_hat(1), malphat_hat(2), malphat_hat(3)));

saveas(gcf, fullfile(outdir,'Tyk_treated_n-alphaandbeta_multdrug_alphaT_constD_10days_model5.png'));

writematrix([Tyk_10t_days, model_fit_all{1}], fullfile(outdir,'Tyk_IC50treated_n-alphaandbeta_multdrug_alphaT_constD_10days_model5.csv'));
writematrix([Tyk_10t_days, model_fit_all{2}], fullfile(outdir,'Tyk_IC25treated_n-alphaandbeta_multdrug_alphaT_constD_10days_model5.csv'));
writematrix([Tyk_10t_days, model_fit_all{3}], fullfile(outdir,'Tyk_IC75treated_n-alphaandbeta_multdrug_alphaT_constD_10days_model5.csv'));

% writematrix(S_raw_all{1,1}, fullfile(outdir, 'Tyk_IC50treated_n-alphaandbeta_multdrug_alphaT_constD_10days.csv'));
% writematrix(S_raw_all{1,2}, fullfile(outdir, 'Tyk_IC25treated_n-alphaandbeta_multdrug_alphaT_constD_10days.csv'));
% writematrix(S_raw_all{1,3}, fullfile(outdir, 'Tyk_IC75treated_n-alphaandbeta_multdrug_alphaT_constD_10days.csv'));

%MODEL 6 (Equation 16)
%%%Estimate for n, alpha and beta - Tyk (constant D) additive hill term,
%%%alpha is time dependent 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_treated_10days';

%Import all data
Tyk_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC50_0.66um_40k_cellsinwell.csv');
Tyk_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC75_0.12um_40k_cellsinwell.csv');
Tyk_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC25_0.85um_40k_cellsinwell.csv');

Tyk_t_days = Tyk_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
Tyk_IC50_y_raw = [Tyk_IC50_data.Rep1, Tyk_IC50_data.Rep2, Tyk_IC50_data.Rep3];
Tyk_IC25_y_raw = [Tyk_IC25_data.Rep1, Tyk_IC25_data.Rep2, Tyk_IC25_data.Rep3];
Tyk_IC75_y_raw = [Tyk_IC75_data.Rep1, Tyk_IC75_data.Rep2, Tyk_IC75_data.Rep3];

%downsample data
Tyk_10t_days = Tyk_t_days(1:11);
Tyk_IC25_10y_raw = Tyk_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
Tyk_IC50_10y_raw = Tyk_IC50_y_raw(1:11,:);
Tyk_IC75_10y_raw = Tyk_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 0.66; 
ce_IC25 = 0.12;
ce_IC75 = 0.85;

%Define logistic growth constants 
rS = 0.6433; %growth rate (day^-1)
KS = 545681.7362;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 1.974574; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [2; 0.5; 0.5];      % [n; alpha; beta]
lb = [0.01; 0.01; 0.01];   % reasonable bounds
ub = [20; 20.0; 20.0]; %change ub for the last one 

res_fun_nandalpha_constD_add_alpha_t = @(p) all_datasets_residual_reps_n_and_alpha_constD_add_alpha_t(p, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

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
fprintf('Fitted n, alpha and beta - add drug - alpha T - constD Tyk (10 days)\n')
fprintf('Fitted Hill coefficient n = %.4f\n', aalphat_hat(1));
fprintf('Fitted alpha term a = %.4f\n', aalphat_hat(2));
fprintf('Fitted beta term b = %.4f\n', aalphat_hat(3));
fprintf('Sum of squared residuals = %.2f\n', resnorm_aalphat);
% 1 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted n, alpha and beta - add drug - alpha T - constD Tyk (10 days)
% Fitted Hill coefficient n = 1.7837
% Fitted alpha term a = 0.0215
% Fitted beta term b = 14.4750
% Sum of squared residuals = 450489628862.92

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_nandalpha_constD_add_alpha_t, aalphat_hat, lb, ub);

% Compute 95% confidence intervals
ci = nlparci(p_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95 Confidence Intervals:\n')
fprintf('n: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
fprintf('a: [%.4f, %.4f]\n', ci(2,1), ci(2,2));
fprintf('b: [%.4f, %.4f]\n', ci(3,1), ci(3,2));
% Solver stopped prematurely.
% 
% lsqnonlin stopped because it exceeded the function evaluation limit,
% options.MaxFunctionEvaluations = 3.000000e+02.
% 
% 95 Confidence Intervals:
% n: [1.2890, 2.3373]
% a: [-0.4023, 0.4394]
% b: [-377.6943, 412.1183]

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_add_alpha_t(aalphat_hat, rS, KS, IC50S, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

fprintf('Model (n=%.2f, a=%.2f, b=%.2f):\n', aalphat_hat(1), aalphat_hat(2), aalphat_hat(3));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (n=1.78, a=0.02, b=14.48):
%   BIC IC50  = 733.94
%   BIC IC25  = 754.56
%   BIC IC75  = 741.19
%   BIC total = 2229.69

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {Tyk_IC50_10y_raw, Tyk_IC25_10y_raw, Tyk_IC75_10y_raw};
model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(Tyk_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_add_alpha_t(aalphat_hat(1), aalphat_hat(2), aalphat_hat(3), Tyk_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(Tyk_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 600000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (n=%.2f, a=%.2f, b=%.2f):\n', aalphat_hat(1), aalphat_hat(2), aalphat_hat(3)));

saveas(gcf, fullfile(outdir,'Tyk_treated_n-alphaandbeta_adddrug_alphaT_constD_10days.png'));

writematrix([Tyk_10t_days, model_fit_all{1}], fullfile(outdir,'Tyk_IC50treated_n-alphaandbeta_adddrug_alphaT_constD_10days_model6.csv'));
writematrix([Tyk_10t_days, model_fit_all{2}], fullfile(outdir,'Tyk_IC25treated_n-alphaandbeta_adddrug_alphaT_constD_10days_model6.csv'));
writematrix([Tyk_10t_days, model_fit_all{3}], fullfile(outdir,'Tyk_IC75treated_n-alphaandbeta_adddrug_alphaT_constD_10days_model6.csv'));


% writematrix(S_raw_all{1,1}, fullfile(outdir, 'Tyk_IC50treated_n-alphaandbeta_adddrug_alphaT_constD_10days.csv'));
% writematrix(S_raw_all{1,2}, fullfile(outdir, 'Tyk_IC25treated_n-alphaandbeta_adddrug_alphaT_constD_10days.csv'));
% writematrix(S_raw_all{1,3}, fullfile(outdir, 'Tyk_IC75treated_n-alphaandbeta_adddrug_alphaT_constD_10days.csv'));


%MODEL 18 (Equation 17)
%%%Estimate for alpha only - Tyk (constant D) multiplicative hill term max,
%%%alpha is time dependent, drug kill is 0 at time 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_treated_10days';

%Import all data
Tyk_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC50_0.66um_40k_cellsinwell.csv');
Tyk_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC75_0.12um_40k_cellsinwell.csv');
Tyk_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC25_0.85um_40k_cellsinwell.csv');

Tyk_t_days = Tyk_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
Tyk_IC50_y_raw = [Tyk_IC50_data.Rep1, Tyk_IC50_data.Rep2, Tyk_IC50_data.Rep3];
Tyk_IC25_y_raw = [Tyk_IC25_data.Rep1, Tyk_IC25_data.Rep2, Tyk_IC25_data.Rep3];
Tyk_IC75_y_raw = [Tyk_IC75_data.Rep1, Tyk_IC75_data.Rep2, Tyk_IC75_data.Rep3];

%downsample data
Tyk_10t_days = Tyk_t_days(1:11);
Tyk_IC25_10y_raw = Tyk_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
Tyk_IC50_10y_raw = Tyk_IC50_y_raw(1:11,:);
Tyk_IC75_10y_raw = Tyk_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 0.66; 
ce_IC25 = 0.12;
ce_IC75 = 0.85;

%Define logistic growth constants 
rS = 0.6433; %growth rate (day^-1)
KS = 545681.7362;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 1.974574; %the IC50 for sensitive cells, after 72hr exposure
Dmax = 50;

%Estimate for parameters
p0 = [0.5];      % [alpha]
lb = [0.0001];   % reasonable bounds
ub = [20];

res_fun_nandalpha_constD_multmax_alphaonly_t0 = @(p) all_datasets_residual_reps_n_and_alpha_constD_multmax_alphaonly_t0(p, rS, KS, IC50S, Dmax, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

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
fprintf('Fitted alpha only - mult drug max - alpha T t0=0 - constD Tyk (10 days)\n')
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted alpha only - mult drug max - alpha T t0=0 - constD Tyk (10 days)
% Fitted alpha term a = 0.0095
% Sum of squared residuals = 541825699651.26

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_nandalpha_constD_multmax_alphaonly_t0, malphat_hat, lb, ub);

% Compute 95% confidence intervals
ci = nlparci(p_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95 Confidence Intervals:\n')
fprintf('a: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
% 95 Confidence Intervals:
% a: [0.0089, 0.0102]

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_multmax_alphaonly_t0(malphat_hat, rS, KS, IC50S, Dmax, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

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
S_raw_all = {Tyk_IC50_10y_raw, Tyk_IC25_10y_raw, Tyk_IC75_10y_raw};
model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(Tyk_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_multmax_alphaonly_t0(malphat_hat(1), Tyk_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S, Dmax);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(Tyk_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 600000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (a=%.2f):\n', malphat_hat(1)));

saveas(gcf, fullfile(outdir,'Tyk_treated_alphaonly_multdrugmax_alphaT0_constD_10days_Model18.png'));

writematrix([Tyk_10t_days, model_fit_all{1}], fullfile(outdir,'Tyk_IC50treated_alphaonly_multdrugmax_alphaT0_constD_10days_Model18.csv'));
writematrix([Tyk_10t_days, model_fit_all{2}], fullfile(outdir,'Tyk_IC25treated_alphaonly_multdrugmax_alphaT0_constD_10days_Model18.csv'));
writematrix([Tyk_10t_days, model_fit_all{3}], fullfile(outdir,'Tyk_IC75treated_alphaonly_multdrugmax_alphaT0_constD_10days_Model18.csv'));

% writematrix(S_raw_all{1,1}, fullfile(outdir, 'Tyk_IC50treated_alphaonly_multdrugmax_alphaT0_constD_10days_Model18.csv'));
% writematrix(S_raw_all{1,2}, fullfile(outdir, 'Tyk_IC25treated_alphaonly_multdrugmax_alphaT0_constD_10days_Model18.csv'));
% writematrix(S_raw_all{1,3}, fullfile(outdir, 'Tyk_IC75treated_alphaonly_multdrugmax_alphaT0_constD_10days_Model18.csv'));

%MODEL 19 (Equation 18)
%%%Estimate for alpha only - Tyk (constant D) additive hill term max,
%%%alpha is time dependent, drug kill is 0 at time 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_treated_10days';

%Import all data
Tyk_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC50_0.66um_40k_cellsinwell.csv');
Tyk_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC75_0.12um_40k_cellsinwell.csv');
Tyk_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC25_0.85um_40k_cellsinwell.csv');

Tyk_t_days = Tyk_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
Tyk_IC50_y_raw = [Tyk_IC50_data.Rep1, Tyk_IC50_data.Rep2, Tyk_IC50_data.Rep3];
Tyk_IC25_y_raw = [Tyk_IC25_data.Rep1, Tyk_IC25_data.Rep2, Tyk_IC25_data.Rep3];
Tyk_IC75_y_raw = [Tyk_IC75_data.Rep1, Tyk_IC75_data.Rep2, Tyk_IC75_data.Rep3];

%downsample data
Tyk_10t_days = Tyk_t_days(1:11);
Tyk_IC25_10y_raw = Tyk_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
Tyk_IC50_10y_raw = Tyk_IC50_y_raw(1:11,:);
Tyk_IC75_10y_raw = Tyk_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 0.66; 
ce_IC25 = 0.12;
ce_IC75 = 0.85;

%Define logistic growth constants 
rS = 0.6433; %growth rate (day^-1)
KS = 545681.7362;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 1.974574; %the IC50 for sensitive cells, after 72hr exposure
Dmax = 50;

%Estimate for parameters
p0 = [0.5];      % [alpha]
lb = [0.0001];   % reasonable bounds
ub = [20];

res_fun_nandalpha_constD_addmax_alphaonly_t0 = @(p) all_datasets_residual_reps_n_and_alpha_constD_addmax_alphaonly_t0(p, rS, KS, IC50S, Dmax, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

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
fprintf('Fitted alpha only - add drug max - alpha T t0=0 - constD Tyk (10 days)\n')
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted alpha only - add drug max - alpha T t0=0 - constD Tyk (10 days)
% Fitted alpha term a = 0.0033
% Sum of squared residuals = 533839592190.80

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_nandalpha_constD_addmax_alphaonly_t0, malphat_hat, lb, ub);

% Compute 95% confidence intervals
ci = nlparci(p_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95 Confidence Intervals:\n')
fprintf('a: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
% Local minimum possible.
% lsqnonlin stopped because the final change in the sum of squares relative to 
% its initial value is less than the value of the function tolerance.
% <stopping criteria details>
% 95 Confidence Intervals:
% a: [0.0029, 0.0038]


[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_addmax_alphaonly_t0(malphat_hat, rS, KS, IC50S, Dmax, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

fprintf('Model (a=%.2f):\n', malphat_hat(1));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (a=0.00):
%   BIC IC50  = 724.54
%   BIC IC25  = 760.01
%   BIC IC75  = 733.56
%   BIC total = 2218.11

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {Tyk_IC50_10y_raw, Tyk_IC25_10y_raw, Tyk_IC75_10y_raw};

model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(Tyk_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_addmax_alphaonly_t0(malphat_hat(1), Tyk_10t_days, S0, ...
            ce_all(i), rS, KS, IC50S, Dmax);
    end
    S_fit_mean = mean(S_fit_all, 2);
    model_fit_all{i} = S_fit_mean;
    plot(Tyk_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 600000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (a=%.2f):\n', malphat_hat(1)));

saveas(gcf, fullfile(outdir,'Tyk_treated_alphaonly_adddrugmax_alphaT0_constD_10days_Model19.png'));

writematrix([Tyk_10t_days, model_fit_all{1}], fullfile(outdir,'Tyk_IC50treated_alphaonly_adddrugmax_alphaT0_constD_10days_Model19.csv'));
writematrix([Tyk_10t_days, model_fit_all{2}], fullfile(outdir,'Tyk_IC25treated_alphaonly_adddrugmax_alphaT0_constD_10days_Model19.csv'));
writematrix([Tyk_10t_days, model_fit_all{3}], fullfile(outdir,'Tyk_IC75treated_alphaonly_adddrugmax_alphaT0_constD_10days_Model19.csv'));

% writematrix(S_raw_all{1,1}, fullfile(outdir, 'Tyk_IC50treated_alphaonly_adddrugmax_alphaT0_constD_10days_Model19.csv'));
% writematrix(S_raw_all{1,2}, fullfile(outdir, 'Tyk_IC25treated_alphaonly_adddrugmax_alphaT0_constD_10days_Model19.csv'));
% writematrix(S_raw_all{1,3}, fullfile(outdir, 'Tyk_IC75treated_alphaonly_adddrugmax_alphaT0_constD_10days_Model19.csv'));

%MODEL 20, equivalent to model 13 but linear (Equation 19)
%%%Estimate for alpha only - Tyk (constant D) additive hill,
%%%alpha is time dependent, drug kill is 0 at time 0, LINEAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_treated_10days';

%Import all data
Tyk_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC50_0.66um_40k_cellsinwell.csv');
Tyk_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC75_0.12um_40k_cellsinwell.csv');
Tyk_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC25_0.85um_40k_cellsinwell.csv');

Tyk_t_days = Tyk_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
Tyk_IC50_y_raw = [Tyk_IC50_data.Rep1, Tyk_IC50_data.Rep2, Tyk_IC50_data.Rep3];
Tyk_IC25_y_raw = [Tyk_IC25_data.Rep1, Tyk_IC25_data.Rep2, Tyk_IC25_data.Rep3];
Tyk_IC75_y_raw = [Tyk_IC75_data.Rep1, Tyk_IC75_data.Rep2, Tyk_IC75_data.Rep3];

%downsample data
Tyk_10t_days = Tyk_t_days(1:11);
Tyk_IC25_10y_raw = Tyk_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
Tyk_IC50_10y_raw = Tyk_IC50_y_raw(1:11,:);
Tyk_IC75_10y_raw = Tyk_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 0.66; 
ce_IC25 = 0.12;
ce_IC75 = 0.85;

%Define logistic growth constants 
rS = 0.6433; %growth rate (day^-1)
KS = 545681.7362;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 1.974574; %the IC50 for sensitive cells, after 72hr exposure
%Dmax = 27.15376;

%Estimate for parameters
p0 = [0.5];      % [alpha]
lb = [0.0001];   % reasonable bounds
ub = [20];

res_fun_nandalpha_constD_linear_alphaonly_t0 = @(p) all_datasets_residual_reps_n_and_alpha_constD_linear_alphaonly_t0(p, rS, KS, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

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
% Fitted alpha term a = 0.0620
% Sum of squared residuals = 500123298738.34

[~, BIC_IC50, BIC_IC25, BIC_IC75] = all_datasets_residual_reps_n_and_alpha_constD_linear_alphaonly_t0(malphat_hat, rS, KS, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

fprintf('Model (a=%.2f):\n', malphat_hat(1));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (a=0.06):
%   BIC IC50  = 725.01
%   BIC IC25  = 756.24
%   BIC IC75  = 732.83
%   BIC total = 2214.07

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {Tyk_IC50_10y_raw, Tyk_IC25_10y_raw, Tyk_IC75_10y_raw};

model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(Tyk_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_linear_alphaonly_t0(malphat_hat(1), Tyk_10t_days, S0, ...
            ce_all(i), rS, KS);
    end
    S_fit_mean = mean(S_fit_all, 2);
    model_fit_all{i} = S_fit_mean;
    plot(Tyk_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 600000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (a=%.2f):\n', malphat_hat(1)));

saveas(gcf, fullfile(outdir,'Tyk_treated_alphaonly_adddruglinear_alphaT0_constD_10days_Model20.png'));

writematrix([Tyk_10t_days, model_fit_all{1}], fullfile(outdir,'Tyk_IC50treated_alphaonly_adddruglinear_alphaT0_constD_10days_Model20.csv'));
writematrix([Tyk_10t_days, model_fit_all{2}], fullfile(outdir,'Tyk_IC25treated_alphaonly_adddruglinear_alphaT0_constD_10days_Model20.csv'));
writematrix([Tyk_10t_days, model_fit_all{3}], fullfile(outdir,'Tyk_IC75treated_alphaonly_adddruglinear_alphaT0_constD_10days_Model20.csv'));

% writematrix(S_raw_all{1,1}, fullfile(outdir, 'Tyk_IC50treated_alphaonly_adddruglinear_alphaT0_constD_10days_Model20.csv'));
% writematrix(S_raw_all{1,2}, fullfile(outdir, 'Tyk_IC25treated_alphaonly_adddruglinear_alphaT0_constD_10days_Model20.csv'));
% writematrix(S_raw_all{1,3}, fullfile(outdir, 'Tyk_IC75treated_alphaonly_adddruglinear_alphaT0_constD_10days_Model20.csv'));


%MODEL 21, equivalent to model 15 and 7 but linear same as 20 (Equation 20)
%%%Estimate for alpha only - Tyk (constant D) additive hill term,
%%%alpha is time dependent, drug kill is 0 at time 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results';

%Import all data
Tyk_IC50_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC50_0.66um_40k_cellsinwell.csv');
Tyk_IC25_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC75_0.12um_40k_cellsinwell.csv');
Tyk_IC75_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_IC25_0.85um_40k_cellsinwell.csv');

Tyk_t_days = Tyk_IC50_data.Day; %no need to define multiple times, time is all the same for all concentrations
Tyk_IC50_y_raw = [Tyk_IC50_data.Rep1, Tyk_IC50_data.Rep2, Tyk_IC50_data.Rep3];
Tyk_IC25_y_raw = [Tyk_IC25_data.Rep1, Tyk_IC25_data.Rep2, Tyk_IC25_data.Rep3];
Tyk_IC75_y_raw = [Tyk_IC75_data.Rep1, Tyk_IC75_data.Rep2, Tyk_IC75_data.Rep3];

%downsample data
Tyk_10t_days = Tyk_t_days(1:11);
Tyk_IC25_10y_raw = Tyk_IC25_y_raw(1:11,:); %rows from 1-9 and all columns
Tyk_IC50_10y_raw = Tyk_IC50_y_raw(1:11,:);
Tyk_IC75_10y_raw = Tyk_IC75_y_raw(1:11,:);

%Define extracellular concentration (uM) values for each condition 
ce_IC50 = 0.66; 
ce_IC25 = 0.12;
ce_IC75 = 0.85;

%Define logistic growth constants 
rS = 0.6433; %growth rate (day^-1)
KS = 545681.7362;  %carrying capacity (cells)

%Define drug parameters 
IC50S = 1.974574; %the IC50 for sensitive cells, after 72hr exposure

%Estimate for parameters
p0 = [0.5];      % [alpha]
lb = [0.0001];   % reasonable bounds
ub = [20];

res_fun_nandalpha_constD_multlinear_alphaonly_t0 = @(p) all_datasets_residual_reps_n_and_alpha_constD_multlinear_alphaonly_t0(p, rS, KS, ...
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

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
fprintf('Fitted alpha only - add drug - alpha T t0=0 - constD Tyk (14 days)\n')
fprintf('Fitted alpha term a = %.4f\n', malphat_hat(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_malphat);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted alpha only - add drug - alpha T t0=0 - constD Tyk (14 days)
% Fitted alpha term a = 0.1740
% Sum of squared residuals = 539037681760.95

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
                                         Tyk_10t_days, Tyk_IC50_10y_raw, ce_IC50, ...
                                         Tyk_10t_days, Tyk_IC25_10y_raw, ce_IC25, ...
                                         Tyk_10t_days, Tyk_IC75_10y_raw, ce_IC75);

fprintf('Model (a=%.2f):\n', malphat_hat(1));
fprintf('  BIC IC50  = %.2f\n', BIC_IC50);
fprintf('  BIC IC25  = %.2f\n', BIC_IC25);
fprintf('  BIC IC75  = %.2f\n', BIC_IC75);
fprintf('  BIC total = %.2f\n', BIC_IC50 + BIC_IC25 + BIC_IC75);
% Model (a=0.17):
%   BIC IC50  = 735.53
%   BIC IC25  = 752.80
%   BIC IC75  = 738.90
%   BIC total = 2227.23

%Put datasets together 
ce_all = [ce_IC50, ce_IC25, ce_IC75];
S_raw_all = {Tyk_IC50_10y_raw, Tyk_IC25_10y_raw, Tyk_IC75_10y_raw};
model_fit_all = cell(1,3);

figure; 
for i = 1:3
    subplot(1,3,i); hold on;
    
    % All 3 replicates (gray circles)
    S_raw = S_raw_all{i};
    [nTimes, nReps] = size(S_raw);
    plot(Tyk_10t_days, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);
    
    % Fitted mean curve (red line)
    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_treat_for_dataset_withConstD_multlinear_alphaonly_t0(malphat_hat(1), Tyk_10t_days, S0, ...
            ce_all(i), rS, KS);
    end
    S_fit_mean = mean(S_fit_all, 2);
    plot(Tyk_10t_days, S_fit_mean, 'r-', 'LineWidth', 3);
    
    title(sprintf('Ce = %.2f μM', ce_all(i)));
    xlabel('Time (days)'); ylabel('Sensitive cells');
    ylim([0 600000]); grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Model (a=%.2f):\n', malphat_hat(1)));

saveas(gcf, fullfile(outdir,'Tyk_treated_alphaonly_multdruglinear_alphaT0_constD_10days_Model21.png'));

writematrix([Tyk_10t_days, model_fit_all{1}], fullfile(outdir,'Tyk_IC50treated_alphaonly_multdruglinear_alphaT0_constD_10days_Model21.csv'));
writematrix([Tyk_10t_days, model_fit_all{2}], fullfile(outdir,'Tyk_IC25treated_alphaonly_multdruglinear_alphaT0_constD_10days_Model21.csv'));
writematrix([Tyk_10t_days, model_fit_all{3}], fullfile(outdir,'Tyk_IC75treated_alphaonly_multdruglinear_alphaT0_constD_10days_Model21.csv'));