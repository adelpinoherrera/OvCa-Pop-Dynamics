%%%Attempting to find parameters that would fit the co-culture data at:
%%att treatment conditions (0,0.62,1,1.47) and all sensitivities (75/25,
%%50/50, 25/75)

%Estimate for only shared K -all treatments, all sensitivities, 10 days 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Import data and format it 
clear all;
outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture';

coculture_files = {'/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/ce0_75S25R_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/ce0_50S50R_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/ce0_25S75R_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/ce1_75S25R_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/ce1_50S50R_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/ce1_25S75R_cellsinwell.csv'}

coculture_data = cell(6,1);
S_raw_all = cell(6,1); R_raw_all = cell(6,1);
ce_coculture = [0, 0, 0, 1.00, 1.00, 1.00];  % 2 treatments × 3 ratios
ratio_labels = {'75/25 S/R', '50/50 S/R', '25/75 S/R'};

%organize data with their different replicates
for i = 1:6
    coculture_data{i} = readtable(coculture_files{i});
    S_raw_all{i} = [coculture_data{i}.S_Rep1, coculture_data{i}.S_Rep2, coculture_data{i}.S_Rep3];
    R_raw_all{i} = [coculture_data{i}.R_Rep1, coculture_data{i}.R_Rep2, coculture_data{i}.R_Rep3];
end
t_coculture = coculture_data{1}.Day;

% Logical mask for days 0–10 (inclusive)
idx_0_10 = (t_coculture >= 0) & (t_coculture <= 10);

% Downsampled time vector
t_coculture_10 = t_coculture(idx_0_10);

% Downsample S and R for each condition
S_raw_10_all = cell(6,1);
R_raw_10_all = cell(6,1);
for i = 1:6
    S_raw_10_all{i} = S_raw_all{i}(idx_0_10, :);
    R_raw_10_all{i} = R_raw_all{i}(idx_0_10, :);
end

%define parameters
ce_IC50 = 1.00;

%Sensitive data
rS = 0.8547; 
KS = 1114404.9924; 
IC50S = 0.9960177; 
nS = 3.7779;
alphaS = 0.129;

%Resistant data
rR = 0.5553; 
KR = 1356851.059; 
IC50R = 8.338535; 
nR = 2.5166;
alphaR = 1.7982;

%run the same code as before, with t_coculture_10, S_raw_10_all and R_raw_10_all
p0 = [1250000]; 
lb = [0.001]; 
ub = [2000000];

res_fun_10 = @(p) coculture_residuals_all_sharedK_m13(p, rS, alphaS, nS, IC50S, rR, alphaR, nR, IC50R, ...
    t_coculture_10, S_raw_10_all, R_raw_10_all, ce_coculture);

problem_gd_10 = createOptimProblem('lsqnonlin', 'x0', p0, 'objective', res_fun_10, 'lb', lb, 'ub', ub);
ms_gd_10 = MultiStart('UseParallel', true, 'Display', 'iter');
nStarts = 200;
[gd_hat_10, resnorm_gd_10] = run(ms_gd_10, problem_gd_10, nStarts);
fprintf('Fitted K (10 days)\n');
fprintf('Fitted K = %.4f\n', gd_hat_10(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_gd_10);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted K (10 days)
% Fitted K = 990795.3111
% Sum of squared residuals = 8741994233637.63

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_10, gd_hat_10, lb, ub);

% Compute 95% confidence intervals
ci = nlparci(p_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95 Confidence Intervals:\n')
fprintf('K: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
% Local minimum possible.
% lsqnonlin stopped because the final change in the sum of squares relative to 
% its initial value is less than the value of the function tolerance.
% <stopping criteria details>
% 95 Confidence Intervals:
% K: [936177.0443, 1045413.6792] (1045413.6792-936177.0443)/990795.3111 =
% 0.1103

%calculate BICs
[~, BIC_conditions_10] = coculture_residuals_all_sharedK_m13( ...
    gd_hat_10, rS, alphaS, nS, IC50S, ...
    rR, alphaR, nR, IC50R, ...
    t_coculture_10, S_raw_10_all, R_raw_10_all, ce_coculture);

%output BIC for all the conditions
fprintf('\nBIC per coculture condition:\n');
for i = 1:6
    if i <= 3
        ce_val = 0;
    else
        ce_val = 1.00;
    end
    ratio_labels = {'75/25 S/R','50/50 S/R','25/75 S/R'};
    ratio_idx = mod(i-1,3) + 1;
    fprintf('  Ce=%.2f, %s: BIC = %.2f\n', ce_val, ratio_labels{ratio_idx}, BIC_conditions_10(i));
end

% BIC per coculture condition:
%   Ce=0.00, 75/25 S/R: BIC = 1563.67
%   Ce=0.00, 50/50 S/R: BIC = 1514.06
%   Ce=0.00, 25/75 S/R: BIC = 1551.60
%   Ce=1.00, 75/25 S/R: BIC = 1586.45
%   Ce=1.00, 50/50 S/R: BIC = 1572.04
%   Ce=1.00, 25/75 S/R: BIC = 1621.40

%plot, 6 plots, 3 sensitivities and 2 treatments 
figure;
for i = 1:6
    S_raw = S_raw_10_all{i};
    R_raw = R_raw_10_all{i};
    [nTimes, nReps] = size(S_raw);

    % Simulate each replicate
    S_mod_all = zeros(nTimes, nReps);
    R_mod_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        R0 = R_raw(1,j);
        [S_mod, R_mod] = simulate_coculture_sharedK_m13( ...
            t_coculture_10, S0, R0, ...
            rS, alphaS, nS, IC50S, ...
            rR, gd_hat_10(1), alphaR, nR, IC50R, ce_coculture(i));
        S_mod_all(:,j) = S_mod;
        R_mod_all(:,j) = R_mod;
    end

    % Mean curves (optional, cleaner than plotting all replicates)
    S_mod_mean = mean(S_mod_all, 2);
    R_mod_mean = mean(R_mod_all, 2);

    subplot(2,3,i); hold on;

    % Data: red dots (S), green dots (R)
    Sdata = plot(t_coculture_10, S_raw, 'ro', 'MarkerFaceColor','r', 'MarkerSize',4);
    Rdata = plot(t_coculture_10, R_raw, 'go', 'MarkerFaceColor','g', 'MarkerSize',4);

    % Model: red line (S), green line (R)
    Smodel = plot(t_coculture_10, S_mod_mean, 'r-', 'LineWidth',2);
    Rmodel = plot(t_coculture_10, R_mod_mean, 'g-', 'LineWidth',2);

    title(sprintf('Ce=%.2f, %s', ...
        ce_coculture(i), ratio_labels{mod(i-1,3)+1}));
    xlabel('Time (days)');
    ylabel('Cell count');
    ylim([0 1400000]);
    grid on;
    legend( [Sdata(1), Rdata(1), Smodel(1), Rmodel(1)],...
        {'S data','R data','S model','R model'}, 'Location','best');
end

sgtitle(sprintf('sharedK = %.3f', gd_hat_10(1)));
saveas(gcf, fullfile(outdir,'co-culture_A2780_ALLtreat_ALLsensitivities_sharedK_m13_constD_10days.png'));

figure;
for i = 1:6
    S_raw = S_raw_10_all{i};
    R_raw = R_raw_10_all{i};
    [nTimes, nReps] = size(S_raw);

    % Simulate each replicate
    S_mod_all = zeros(nTimes, nReps);
    R_mod_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        R0 = R_raw(1,j);
        [S_mod, R_mod] = simulate_coculture_sharedK_m13( ...
            t_coculture_10, S0, R0, ...
            rS, alphaS, nS, IC50S,  ...
            rR, gd_hat_10(1), alphaR, nR, IC50R, ce_coculture(i));
        S_mod_all(:,j) = S_mod;
        R_mod_all(:,j) = R_mod;
    end

    % Means and SDs across replicates, mean data
    S_mean = mean(S_raw,  2);
    S_sd   = std(S_raw, 0, 2);
    R_mean = mean(R_raw,  2);
    R_sd   = std(R_raw, 0, 2);

    % Mean curves (optional, cleaner than plotting all replicates)
    S_mod_mean = mean(S_mod_all, 2);
    R_mod_mean = mean(R_mod_all, 2);

    % Totals (data and model)
    Tot_mean      = S_mean + R_mean;          % per replicate
    Tot_mod_mean = S_mod_mean + R_mod_mean;

writematrix([t_coculture_10, S_mod_all], ...
    fullfile(outdir, sprintf('S_sim_i%d_Ce%.2f_sharedK_m13_10days.csv', i, ce_coculture(i))));

writematrix([t_coculture_10, R_mod_all], ...
    fullfile(outdir, sprintf('R_sim_i%d_Ce%.2f_sharedK_m13_10days.csv', i, ce_coculture(i))));
    
subplot(2,3,i); hold on;

    % Data: red dots (S), green dots (R)
    Sdata = errorbar(t_coculture_10, S_mean, S_sd, 'r.', 'MarkerSize', 12, 'LineStyle','none');
    Rdata = errorbar(t_coculture_10, R_mean, R_sd, 'g.', 'MarkerSize',12, 'LineStyle','none');
    Tdata  = plot(t_coculture_10, Tot_mean,'b.', 'MarkerSize',12);

    % Model: red line (S), green line (R)
    Smodel = plot(t_coculture_10, S_mod_mean, 'r-', 'LineWidth',2);
    Rmodel = plot(t_coculture_10, R_mod_mean, 'g-', 'LineWidth',2);
    Tmodel = plot(t_coculture_10, Tot_mod_mean, 'b-', 'LineWidth',2);

    title(sprintf('Ce=%.2f, %s', ...
        ce_coculture(i), ratio_labels{mod(i-1,3)+1}));
    xlabel('Time (days)');
    ylabel('Cell count');
    ylim([0 1400000]);
    grid on;
    legend( [Sdata(1), Rdata(1), Tdata(1), Smodel(1), Rmodel(1), Tmodel(1)],...
        {'S data','R data','Total data','S model','R model','Total model'}, 'Location','best');
end
sgtitle(sprintf('sharedK = %.3f', gd_hat_10(1)));

saveas(gcf, fullfile(outdir,'co-culture_A2780withTotal_ALLtreat_ALLsensitivities_sharedK_m13_constD_10days.png'));

%calculate R squared values 
clear all;
outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results';

coculture_files = {'/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/ce0_75S25R_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/ce0_50S50R_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/ce0_25S75R_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/ce1_75S25R_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/ce1_50S50R_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/ce1_25S75R_cellsinwell.csv'}

coculture_data = cell(6,1);
S_raw_all = cell(6,1); R_raw_all = cell(6,1);
ce_coculture = [0, 0, 0, 1.00, 1.00, 1.00];  % 2 treatments × 3 ratios
ratio_labels = {'75/25 S/R', '50/50 S/R', '25/75 S/R'};

%organize data with their different replicates
for i = 1:6
    coculture_data{i} = readtable(coculture_files{i});
    S_raw_all{i} = [coculture_data{i}.S_Rep1, coculture_data{i}.S_Rep2, coculture_data{i}.S_Rep3];
    R_raw_all{i} = [coculture_data{i}.R_Rep1, coculture_data{i}.R_Rep2, coculture_data{i}.R_Rep3];
end
t_coculture = coculture_data{1}.Day;

% Logical mask for days 0–10 (inclusive)
idx_0_10 = (t_coculture >= 0) & (t_coculture <= 10);

% Downsampled time vector
t_coculture_10 = t_coculture(idx_0_10);

% Downsample S and R for each condition
S_raw_10_all = cell(6,1);
R_raw_10_all = cell(6,1);
for i = 1:6
    S_raw_10_all{i} = S_raw_all{i}(idx_0_10, :);
    R_raw_10_all{i} = R_raw_all{i}(idx_0_10, :);
end


%sensitive data
A2780_75_Ce0_all = S_raw_10_all{1};
A2780_75_Ce0_all = A2780_75_Ce0_all(:);

A2780_50_Ce0_all = S_raw_10_all{2};
A2780_50_Ce0_all = A2780_50_Ce0_all(:);

A2780_25_Ce0_all = S_raw_10_all{3};
A2780_25_Ce0_all = A2780_25_Ce0_all(:);

A2780_75_Ce1_all = S_raw_10_all{4};
A2780_75_Ce1_all = A2780_75_Ce1_all(:);

A2780_50_Ce1_all = S_raw_10_all{5};
A2780_50_Ce1_all = A2780_50_Ce1_all(:);

A2780_25_Ce1_all = S_raw_10_all{6};
A2780_25_Ce1_all = A2780_25_Ce1_all(:);

%resistant data
A2780cis_75_Ce0_all = R_raw_10_all{1};
A2780cis_75_Ce0_all = A2780cis_75_Ce0_all(:);

A2780cis_50_Ce0_all = R_raw_10_all{2};
A2780cis_50_Ce0_all = A2780cis_50_Ce0_all(:);

A2780cis_25_Ce0_all = R_raw_10_all{3};
A2780cis_25_Ce0_all = A2780cis_25_Ce0_all(:);

A2780cis_75_Ce1_all = R_raw_10_all{4};
A2780cis_75_Ce1_all = A2780cis_75_Ce1_all(:);

A2780cis_50_Ce1_all = R_raw_10_all{5};
A2780cis_50_Ce1_all = A2780cis_50_Ce1_all(:);

A2780cis_25_Ce1_all = R_raw_10_all{6};
A2780cis_25_Ce1_all = A2780cis_25_Ce1_all(:);

%read the predicted data
%Sensitive, 75%, Ce=0.0
A2780_75_Ce0_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture/S_sim_i1_Ce0.00_sharedK_m13_10days.csv');
A2780_75_Ce0_predicted = table2array(A2780_75_Ce0_predicted(:,2)); %only the second column 
A2780_75_Ce0_predicted = repmat(A2780_75_Ce0_predicted,3,1); %repeat 3 times 

SS_res_A2780_75_Ce0 = sum((A2780_75_Ce0_all - A2780_75_Ce0_predicted).^2);
SS_tot_A2780_75_Ce0 = sum((A2780_75_Ce0_all - mean(A2780_75_Ce0_all)).^2);
R2_A2780_75_Ce0 = 1 - SS_res_A2780_75_Ce0 / SS_tot_A2780_75_Ce0;

fprintf('R-squared for A2780, Ce=0, 75: %.4f\n', R2_A2780_75_Ce0);
%R-squared for A2780, Ce=0, 75: 0.5716

%Sensitive, 50%, Ce=0.0
A2780_50_Ce0_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture/S_sim_i2_Ce0.00_sharedK_m13_10days.csv');
A2780_50_Ce0_predicted = table2array(A2780_50_Ce0_predicted(:,2)); %only the second column 
A2780_50_Ce0_predicted = repmat(A2780_50_Ce0_predicted,3,1); %repeat 3 times 

SS_res_A2780_50_Ce0 = sum((A2780_50_Ce0_all - A2780_50_Ce0_predicted).^2);
SS_tot_A2780_50_Ce0 = sum((A2780_50_Ce0_all - mean(A2780_50_Ce0_all)).^2);
R2_A2780_50_Ce0 = 1 - SS_res_A2780_50_Ce0 / SS_tot_A2780_50_Ce0;

fprintf('R-squared for A2780, Ce=0, 50: %.4f\n', R2_A2780_50_Ce0);
%R-squared for A2780, Ce=0, 50: 0.7689

%Sensitive, 25%, Ce=0.0
A2780_25_Ce0_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture/S_sim_i3_Ce0.00_sharedK_m13_10days.csv');
A2780_25_Ce0_predicted = table2array(A2780_25_Ce0_predicted(:,2)); %only the second column 
A2780_25_Ce0_predicted = repmat(A2780_25_Ce0_predicted,3,1); %repeat 3 times 

%delete certain rows
% rowsToRemove = [4,15,26];
% A2780_25_Ce0_all(rowsToRemove) = [];
% A2780_25_Ce0_predicted(rowsToRemove) = [];

SS_res_A2780_25_Ce0 = sum((A2780_25_Ce0_all - A2780_25_Ce0_predicted).^2);
SS_tot_A2780_25_Ce0 = sum((A2780_25_Ce0_all - mean(A2780_25_Ce0_all)).^2);
R2_A2780_25_Ce0 = 1 - SS_res_A2780_25_Ce0 / SS_tot_A2780_25_Ce0;

fprintf('R-squared for A2780, Ce=0, 25: %.4f\n', R2_A2780_25_Ce0);
%R-squared for A2780, Ce=0, 25: 0.7933


%Sensitive, 75%, Ce=1.00
A2780_75_Ce1_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture/S_sim_i4_Ce1.00_sharedK_m13_10days.csv');
A2780_75_Ce1_predicted = table2array(A2780_75_Ce1_predicted(:,2)); %only the second column 
A2780_75_Ce1_predicted = repmat(A2780_75_Ce1_predicted,3,1); %repeat 3 times 

SS_res_A2780_75_Ce1 = sum((A2780_75_Ce1_all - A2780_75_Ce1_predicted).^2);
SS_tot_A2780_75_Ce1 = sum((A2780_75_Ce1_all - mean(A2780_75_Ce1_all)).^2);
R2_A2780_75_Ce1 = 1 - SS_res_A2780_75_Ce1 / SS_tot_A2780_75_Ce1;

fprintf('R-squared for A2780, Ce=1.00, 75: %.4f\n', R2_A2780_75_Ce1);
%R-squared for A2780, Ce=1.00, 75: 0.2549

%Sensitive, 50%, Ce=1.00
A2780_50_Ce1_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture/S_sim_i5_Ce1.00_sharedK_m13_10days.csv');
A2780_50_Ce1_predicted = table2array(A2780_50_Ce1_predicted(:,2)); %only the second column 
A2780_50_Ce1_predicted = repmat(A2780_50_Ce1_predicted,3,1); %repeat 3 times 

SS_res_A2780_50_Ce1 = sum((A2780_50_Ce1_all - A2780_50_Ce1_predicted).^2);
SS_tot_A2780_50_Ce1 = sum((A2780_50_Ce1_all - mean(A2780_50_Ce1_all)).^2);
R2_A2780_50_Ce1 = 1 - SS_res_A2780_50_Ce1 / SS_tot_A2780_50_Ce1;

fprintf('R-squared for A2780, Ce=1.00, 50: %.4f\n', R2_A2780_50_Ce1);
%R-squared for A2780, Ce=0.66, 50: 0.4415

%Sensitive, 25%, Ce=1.00
A2780_25_Ce1_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture/S_sim_i6_Ce1.00_sharedK_m13_10days.csv');
A2780_25_Ce1_predicted = table2array(A2780_25_Ce1_predicted(:,2)); %only the second column 
A2780_25_Ce1_predicted = repmat(A2780_25_Ce1_predicted,3,1); %repeat 3 times 

SS_res_A2780_25_Ce1 = sum((A2780_25_Ce1_all - A2780_25_Ce1_predicted).^2);
SS_tot_A2780_25_Ce1 = sum((A2780_25_Ce1_all - mean(A2780_25_Ce1_all)).^2);
R2_A2780_25_Ce1 = 1 - SS_res_A2780_25_Ce1 / SS_tot_A2780_25_Ce1;

fprintf('R-squared for A2780, Ce=1.00, 25: %.4f\n', R2_A2780_25_Ce1);
%R-squared for Tyk, Ce=0.66, 25: 0.5159


%Resistant, 75%, Ce=0.0
A2780cis_75_Ce0_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture/R_sim_i1_Ce0.00_sharedK_m13_10days.csv');
A2780cis_75_Ce0_predicted = table2array(A2780cis_75_Ce0_predicted(:,2)); %only the second column 
A2780cis_75_Ce0_predicted = repmat(A2780cis_75_Ce0_predicted,3,1); %repeat 3 times 

SS_res_A2780cis_75_Ce0 = sum((A2780cis_75_Ce0_all - A2780cis_75_Ce0_predicted).^2);
SS_tot_A2780cis_75_Ce0 = sum((A2780cis_75_Ce0_all - mean(A2780cis_75_Ce0_all)).^2);
R2_A2780cis_75_Ce0 = 1 - SS_res_A2780cis_75_Ce0 / SS_tot_A2780cis_75_Ce0;

fprintf('R-squared for A2780cis, Ce=0, 75: %.4f\n', R2_A2780cis_75_Ce0);
%R-squared for A2780cis, Ce=0, 75: 0.8187

%Resistant, 50%, Ce=0.0
A2780cis_50_Ce0_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture/R_sim_i2_Ce0.00_sharedK_m13_10days.csv');
A2780cis_50_Ce0_predicted = table2array(A2780cis_50_Ce0_predicted(:,2)); %only the second column 
A2780cis_50_Ce0_predicted = repmat(A2780cis_50_Ce0_predicted,3,1); %repeat 3 times 

SS_res_A2780cis_50_Ce0 = sum((A2780cis_50_Ce0_all - A2780cis_50_Ce0_predicted).^2);
SS_tot_A2780cis_50_Ce0 = sum((A2780cis_50_Ce0_all - mean(A2780cis_50_Ce0_all)).^2);
R2_A2780cis_50_Ce0 = 1 - SS_res_A2780cis_50_Ce0 / SS_tot_A2780cis_50_Ce0;

fprintf('R-squared for A2780cis, Ce=0, 50: %.4f\n', R2_A2780cis_50_Ce0);
%R-squared for A2780cis, Ce=0, 50: 0.8365

%Resistant, 25%, Ce=0.0
A2780cis_25_Ce0_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture/R_sim_i3_Ce0.00_sharedK_m13_10days.csv');
A2780cis_25_Ce0_predicted = table2array(A2780cis_25_Ce0_predicted(:,2)); %only the second column 
A2780cis_25_Ce0_predicted = repmat(A2780cis_25_Ce0_predicted,3,1); %repeat 3 times 

SS_res_A2780cis_25_Ce0 = sum((A2780cis_25_Ce0_all - A2780cis_25_Ce0_predicted).^2);
SS_tot_A2780cis_25_Ce0 = sum((A2780cis_25_Ce0_all - mean(A2780cis_25_Ce0_all)).^2);
R2_A2780cis_25_Ce0 = 1 - SS_res_A2780cis_25_Ce0 / SS_tot_A2780cis_25_Ce0;

fprintf('R-squared for A2780cis, Ce=0, 25: %.4f\n', R2_A2780cis_25_Ce0);
%R-squared for A2780cis, Ce=0, 25: 0.6429


%Resistant, 75%, Ce=1.00
A2780cis_75_Ce1_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture/R_sim_i4_Ce1.00_sharedK_m13_10days.csv');
A2780cis_75_Ce1_predicted = table2array(A2780cis_75_Ce1_predicted(:,2)); %only the second column 
A2780cis_75_Ce1_predicted = repmat(A2780cis_75_Ce1_predicted,3,1); %repeat 3 times 

SS_res_A2780cis_75_Ce1 = sum((A2780cis_75_Ce1_all - A2780cis_75_Ce1_predicted).^2);
SS_tot_A2780cis_75_Ce1 = sum((A2780cis_75_Ce1_all - mean(A2780cis_75_Ce1_all)).^2);
R2_A2780cis_75_Ce1 = 1 - SS_res_A2780cis_75_Ce1 / SS_tot_A2780cis_75_Ce1;

fprintf('R-squared for A2780cis, Ce=1.00, 75: %.4f\n', R2_A2780cis_75_Ce1);
%R-squared for A2780cis, Ce=0.66, 75: 0.5320

%Resistant, 50%, Ce=1.00
A2780cis_50_Ce1_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture/R_sim_i5_Ce1.00_sharedK_m13_10days.csv');
A2780cis_50_Ce1_predicted = table2array(A2780cis_50_Ce1_predicted(:,2)); %only the second column 
A2780cis_50_Ce1_predicted = repmat(A2780cis_50_Ce1_predicted,3,1); %repeat 3 times 

SS_res_A2780cis_50_Ce1 = sum((A2780cis_50_Ce1_all - A2780cis_50_Ce1_predicted).^2);
SS_tot_A2780cis_50_Ce1 = sum((A2780cis_50_Ce1_all - mean(A2780cis_50_Ce1_all)).^2);
R2_A2780cis_50_Ce1 = 1 - SS_res_A2780cis_50_Ce1 / SS_tot_A2780cis_50_Ce1;

fprintf('R-squared for A2780cis, Ce=1.00, 50: %.4f\n', R2_A2780cis_50_Ce1);
%R-squared for A2780cis, Ce=1.00, 50: 0.5591

%Resistant, 25%, Ce=1.00
A2780cis_25_Ce1_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture/R_sim_i6_Ce1.00_sharedK_m13_10days.csv');
A2780cis_25_Ce1_predicted = table2array(A2780cis_25_Ce1_predicted(:,2)); %only the second column 
A2780cis_25_Ce1_predicted = repmat(A2780cis_25_Ce1_predicted,3,1); %repeat 3 times 

SS_res_A2780cis_25_Ce1 = sum((A2780cis_25_Ce1_all - A2780cis_25_Ce1_predicted).^2);
SS_tot_A2780cis_25_Ce1 = sum((A2780cis_25_Ce1_all - mean(A2780cis_25_Ce1_all)).^2);
R2_A2780cis_25_Ce1 = 1 - SS_res_A2780cis_25_Ce1 / SS_tot_A2780cis_25_Ce1;

fprintf('R-squared for A2780cis, Ce=1.00, 25: %.4f\n', R2_A2780cis_25_Ce1);
%R-squared for Tykcpr, Ce=0.66, 25: 0.4504




