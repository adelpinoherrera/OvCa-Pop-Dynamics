%%%Attempting to find parameters that would fit the co-culture data at:
%%att treatment conditions (0,0.12,0.66,0.85) and all sensitivities (75/25,
%%50/50, 25/75)

%Estimate for only shared K -all treatments, all sensitivities, 10 days 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Import data and format it 
clear all;
outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_coculture';

coculture_files = {'/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/ce0_75S25R_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/ce0_50S50R_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/ce0_25S75R_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/ce1_75S25R_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/ce1_50S50R_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/ce1_25S75R_cellsinwell.csv'}

coculture_data = cell(6,1);
S_raw_all = cell(6,1); R_raw_all = cell(6,1);
ce_coculture = [0, 0, 0, 0.66, 0.66, 0.66];  % 2 treatments × 3 ratios
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
ce_IC50 = 0.66;

%Sensitive data
rS = 0.6433; 
KS = 545681.7362; 
IC50S = 1.974574;
alphaS = 0.1666;

%Resistant data
rR = 0.5339; 
KR = 421668.8619; 
IC50R = 4.210567;
alphaR = 0.0006;

%run the same code as before, with t_coculture_10, S_raw_10_all and R_raw_10_all
p0 = [500000]; 
lb = [0.001]; 
ub = [2000000];

res_fun_10 = @(p) coculture_residuals_all_sharedK_m14(p, rS, alphaS, IC50S, rR, alphaR, IC50R, ...
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
% Fitted K = 557834.9748
% Sum of squared residuals = 2788621515540.01

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_10, gd_hat_10, lb, ub);

% Compute 95% confidence intervals
ci = nlparci(p_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95 Confidence Intervals:\n')
fprintf('K: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
% 95 Confidence Intervals:
% K: [526483.9485, 589186.0278]

%calculate BICs
[~, BIC_conditions_10] = coculture_residuals_all_sharedK_m14( ...
    gd_hat_10, rS, alphaS, IC50S, ...
    rR, alphaR, IC50R, ...
    t_coculture_10, S_raw_10_all, R_raw_10_all, ce_coculture);

%output BIC for all the conditions
fprintf('\nBIC per coculture condition:\n');
for i = 1:6
    if i <= 3
        ce_val = 0;
    else
        ce_val = 0.66;
    end
    ratio_labels = {'75/25 S/R','50/50 S/R','25/75 S/R'};
    ratio_idx = mod(i-1,3) + 1;
    fprintf('  Ce=%.2f, %s: BIC = %.2f\n', ce_val, ratio_labels{ratio_idx}, BIC_conditions_10(i));
end

% BIC per coculture condition:
%   Ce=0.00, 75/25 S/R: BIC = 1476.35
%   Ce=0.00, 50/50 S/R: BIC = 1480.93
%   Ce=0.00, 25/75 S/R: BIC = 1496.18
%   Ce=0.66, 75/25 S/R: BIC = 1545.12
%   Ce=0.66, 50/50 S/R: BIC = 1482.25
%   Ce=0.66, 25/75 S/R: BIC = 1495.32

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
        [S_mod, R_mod] = simulate_coculture_sharedK_m14( ...
            t_coculture_10, S0, R0, ...
            rS, alphaS, IC50S, ...
            rR, gd_hat_10(1), alphaR, IC50R, ce_coculture(i));
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
    ylim([0 600000]);
    grid on;
    legend( [Sdata(1), Rdata(1), Smodel(1), Rmodel(1)],...
        {'S data','R data','S model','R model'}, 'Location','best');
end

sgtitle(sprintf('sharedK = %.3f', gd_hat_10(1)));
saveas(gcf, fullfile(outdir,'co-culture_Tyk_ALLtreat_ALLsensitivities_sharedK_m14_constD_10days_model6.png'));

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
        [S_mod, R_mod] = simulate_coculture_sharedK_m14( ...
            t_coculture_10, S0, R0, ...
            rS, alphaS, IC50S,  ...
            rR, gd_hat_10(1), alphaR, IC50R, ce_coculture(i));
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
    fullfile(outdir, sprintf('S_sim_i%d_Ce%.2f_sharedK_m14_10days.csv', i, ce_coculture(i))));

writematrix([t_coculture_10, R_mod_all], ...
    fullfile(outdir, sprintf('R_sim_i%d_Ce%.2f_sharedK_m14_10days.csv', i, ce_coculture(i))));
    
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
    ylim([0 600000]);
    grid on;
    legend( [Sdata(1), Rdata(1), Tdata(1), Smodel(1), Rmodel(1), Tmodel(1)],...
        {'S data','R data','Total data','S model','R model','Total model'}, 'Location','best');
end
sgtitle(sprintf('sharedK = %.3f', gd_hat_10(1)));

saveas(gcf, fullfile(outdir,'co-culture_TykwithTotal_ALLtreat_ALLsensitivities_sharedK_m14_constD_10days_model6.png'));

%%%%R squared values calculation 
clear all;

coculture_files = {'/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/ce0_75S25R_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/ce0_50S50R_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/ce0_25S75R_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/ce1_75S25R_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/ce1_50S50R_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/ce1_25S75R_cellsinwell.csv'}

coculture_data = cell(6,1);
S_raw_all = cell(6,1); R_raw_all = cell(6,1);
ce_coculture = [0, 0, 0, 0.66, 0.66, 0.66];  % 2 treatments × 3 ratios
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
Tyk_75_Ce0_all = S_raw_10_all{1};
Tyk_75_Ce0_all = Tyk_75_Ce0_all(:);

Tyk_50_Ce0_all = S_raw_10_all{2};
Tyk_50_Ce0_all = Tyk_50_Ce0_all(:);

Tyk_25_Ce0_all = S_raw_10_all{3};
Tyk_25_Ce0_all = Tyk_25_Ce0_all(:);

Tyk_75_Ce066_all = S_raw_10_all{4};
Tyk_75_Ce066_all = Tyk_75_Ce066_all(:);

Tyk_50_Ce066_all = S_raw_10_all{5};
Tyk_50_Ce066_all = Tyk_50_Ce066_all(:);

Tyk_25_Ce066_all = S_raw_10_all{6};
Tyk_25_Ce066_all = Tyk_25_Ce066_all(:);

%resistant data
Tykcpr_75_Ce0_all = R_raw_10_all{1};
Tykcpr_75_Ce0_all = Tykcpr_75_Ce0_all(:);

Tykcpr_50_Ce0_all = R_raw_10_all{2};
Tykcpr_50_Ce0_all = Tykcpr_50_Ce0_all(:);

Tykcpr_25_Ce0_all = R_raw_10_all{3};
Tykcpr_25_Ce0_all = Tykcpr_25_Ce0_all(:);

Tykcpr_75_Ce066_all = R_raw_10_all{4};
Tykcpr_75_Ce066_all = Tykcpr_75_Ce066_all(:);

Tykcpr_50_Ce066_all = R_raw_10_all{5};
Tykcpr_50_Ce066_all = Tykcpr_50_Ce066_all(:);

Tykcpr_25_Ce066_all = R_raw_10_all{6};
Tykcpr_25_Ce066_all = Tykcpr_25_Ce066_all(:);

%read the predicted data
%Sensitive, 75%, Ce=0.0
Tyk_75_Ce0_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_coculture/S_sim_i1_Ce0.00_sharedK_m14_10days.csv');
Tyk_75_Ce0_predicted = table2array(Tyk_75_Ce0_predicted(:,2)); %only the second column 
Tyk_75_Ce0_predicted = repmat(Tyk_75_Ce0_predicted,3,1); %repeat 3 times 

SS_res_Tyk_75_Ce0 = sum((Tyk_75_Ce0_all - Tyk_75_Ce0_predicted).^2);
SS_tot_Tyk_75_Ce0 = sum((Tyk_75_Ce0_all - mean(Tyk_75_Ce0_all)).^2);
R2_Tyk_75_Ce0 = 1 - SS_res_Tyk_75_Ce0 / SS_tot_Tyk_75_Ce0;

fprintf('R-squared for Tyk, Ce=0, 75: %.4f\n', R2_Tyk_75_Ce0);
%R-squared for Tyk, Ce=0, 75: 0.7200

%Sensitive, 50%, Ce=0.0
Tyk_50_Ce0_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_coculture/S_sim_i2_Ce0.00_sharedK_m14_10days.csv');
Tyk_50_Ce0_predicted = table2array(Tyk_50_Ce0_predicted(:,2)); %only the second column 
Tyk_50_Ce0_predicted = repmat(Tyk_50_Ce0_predicted,3,1); %repeat 3 times 

SS_res_Tyk_50_Ce0 = sum((Tyk_50_Ce0_all - Tyk_50_Ce0_predicted).^2);
SS_tot_Tyk_50_Ce0 = sum((Tyk_50_Ce0_all - mean(Tyk_50_Ce0_all)).^2);
R2_Tyk_50_Ce0 = 1 - SS_res_Tyk_50_Ce0 / SS_tot_Tyk_50_Ce0;

fprintf('R-squared for Tyk, Ce=0, 50: %.4f\n', R2_Tyk_50_Ce0);
%R-squared for Tyk, Ce=0, 50: 0.6457

%Sensitive, 25%, Ce=0.0
Tyk_25_Ce0_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_coculture/S_sim_i3_Ce0.00_sharedK_m14_10days.csv');
Tyk_25_Ce0_predicted = table2array(Tyk_25_Ce0_predicted(:,2)); %only the second column 
Tyk_25_Ce0_predicted = repmat(Tyk_25_Ce0_predicted,3,1); %repeat 3 times 

%delete certain rows
% rowsToRemove = [4,15,26];
% Tyk_25_Ce0_all(rowsToRemove) = [];
% Tyk_25_Ce0_predicted(rowsToRemove) = [];

SS_res_Tyk_25_Ce0 = sum((Tyk_25_Ce0_all - Tyk_25_Ce0_predicted).^2);
SS_tot_Tyk_25_Ce0 = sum((Tyk_25_Ce0_all - mean(Tyk_25_Ce0_all)).^2);
R2_Tyk_25_Ce0 = 1 - SS_res_Tyk_25_Ce0 / SS_tot_Tyk_25_Ce0;

fprintf('R-squared for Tyk, Ce=0, 25: %.4f\n', R2_Tyk_25_Ce0);
%R-squared for Tyk, Ce=0, 25: -0.7388


%Sensitive, 75%, Ce=0.66
Tyk_75_Ce066_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_coculture/S_sim_i4_Ce0.66_sharedK_m14_10days.csv');
Tyk_75_Ce066_predicted = table2array(Tyk_75_Ce066_predicted(:,2)); %only the second column 
Tyk_75_Ce066_predicted = repmat(Tyk_75_Ce066_predicted,3,1); %repeat 3 times 

SS_res_Tyk_75_Ce066 = sum((Tyk_75_Ce066_all - Tyk_75_Ce066_predicted).^2);
SS_tot_Tyk_75_Ce066 = sum((Tyk_75_Ce066_all - mean(Tyk_75_Ce066_all)).^2);
R2_Tyk_75_Ce066 = 1 - SS_res_Tyk_75_Ce066 / SS_tot_Tyk_75_Ce066;

fprintf('R-squared for Tyk, Ce=0.66, 75: %.4f\n', R2_Tyk_75_Ce066);
%R-squared for Tyk, Ce=0.66, 75: -3.7276

%Sensitive, 50%, Ce=0.66
Tyk_50_Ce066_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_coculture/S_sim_i5_Ce0.66_sharedK_m14_10days.csv');
Tyk_50_Ce066_predicted = table2array(Tyk_50_Ce066_predicted(:,2)); %only the second column 
Tyk_50_Ce066_predicted = repmat(Tyk_50_Ce066_predicted,3,1); %repeat 3 times 

SS_res_Tyk_50_Ce066 = sum((Tyk_50_Ce066_all - Tyk_50_Ce066_predicted).^2);
SS_tot_Tyk_50_Ce066 = sum((Tyk_50_Ce066_all - mean(Tyk_50_Ce066_all)).^2);
R2_Tyk_50_Ce066 = 1 - SS_res_Tyk_50_Ce066 / SS_tot_Tyk_50_Ce066;

fprintf('R-squared for Tyk, Ce=0.66, 50: %.4f\n', R2_Tyk_50_Ce066);
%R-squared for Tyk, Ce=0.66, 50: -3.4380

%Sensitive, 25%, Ce=0.66
Tyk_25_Ce066_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_coculture/S_sim_i6_Ce0.66_sharedK_m14_10days.csv');
Tyk_25_Ce066_predicted = table2array(Tyk_25_Ce066_predicted(:,2)); %only the second column 
Tyk_25_Ce066_predicted = repmat(Tyk_25_Ce066_predicted,3,1); %repeat 3 times 

SS_res_Tyk_25_Ce066 = sum((Tyk_25_Ce066_all - Tyk_25_Ce066_predicted).^2);
SS_tot_Tyk_25_Ce066 = sum((Tyk_25_Ce066_all - mean(Tyk_25_Ce066_all)).^2);
R2_Tyk_25_Ce066 = 1 - SS_res_Tyk_25_Ce066 / SS_tot_Tyk_25_Ce066;

fprintf('R-squared for Tyk, Ce=0.66, 25: %.4f\n', R2_Tyk_25_Ce066);
%R-squared for Tyk, Ce=0.66, 25: -3.3739


%Resistant, 75%, Ce=0.0
Tykcpr_75_Ce0_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_coculture/R_sim_i1_Ce0.00_sharedK_m14_10days.csv');
Tykcpr_75_Ce0_predicted = table2array(Tykcpr_75_Ce0_predicted(:,2)); %only the second column 
Tykcpr_75_Ce0_predicted = repmat(Tykcpr_75_Ce0_predicted,3,1); %repeat 3 times 

SS_res_Tykcpr_75_Ce0 = sum((Tykcpr_75_Ce0_all - Tykcpr_75_Ce0_predicted).^2);
SS_tot_Tykcpr_75_Ce0 = sum((Tykcpr_75_Ce0_all - mean(Tykcpr_75_Ce0_all)).^2);
R2_Tykcpr_75_Ce0 = 1 - SS_res_Tykcpr_75_Ce0 / SS_tot_Tykcpr_75_Ce0;

fprintf('R-squared for Tykcpr, Ce=0, 75: %.4f\n', R2_Tykcpr_75_Ce0);
%R-squared for Tykcpr, Ce=0, 75: 0.4224

%Resistant, 50%, Ce=0.0
Tykcpr_50_Ce0_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_coculture/R_sim_i2_Ce0.00_sharedK_m14_10days.csv');
Tykcpr_50_Ce0_predicted = table2array(Tykcpr_50_Ce0_predicted(:,2)); %only the second column 
Tykcpr_50_Ce0_predicted = repmat(Tykcpr_50_Ce0_predicted,3,1); %repeat 3 times 

SS_res_Tykcpr_50_Ce0 = sum((Tykcpr_50_Ce0_all - Tykcpr_50_Ce0_predicted).^2);
SS_tot_Tykcpr_50_Ce0 = sum((Tykcpr_50_Ce0_all - mean(Tykcpr_50_Ce0_all)).^2);
R2_Tykcpr_50_Ce0 = 1 - SS_res_Tykcpr_50_Ce0 / SS_tot_Tykcpr_50_Ce0;

fprintf('R-squared for Tykcpr, Ce=0, 50: %.4f\n', R2_Tykcpr_50_Ce0);
%R-squared for Tykcpr, Ce=0, 50: 0.6272

%Resistant, 25%, Ce=0.0
Tykcpr_25_Ce0_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_coculture/R_sim_i3_Ce0.00_sharedK_m14_10days.csv');
Tykcpr_25_Ce0_predicted = table2array(Tykcpr_25_Ce0_predicted(:,2)); %only the second column 
Tykcpr_25_Ce0_predicted = repmat(Tykcpr_25_Ce0_predicted,3,1); %repeat 3 times 

SS_res_Tykcpr_25_Ce0 = sum((Tykcpr_25_Ce0_all - Tykcpr_25_Ce0_predicted).^2);
SS_tot_Tykcpr_25_Ce0 = sum((Tykcpr_25_Ce0_all - mean(Tykcpr_25_Ce0_all)).^2);
R2_Tykcpr_25_Ce0 = 1 - SS_res_Tykcpr_25_Ce0 / SS_tot_Tykcpr_25_Ce0;

fprintf('R-squared for Tykcpr, Ce=0, 25: %.4f\n', R2_Tykcpr_25_Ce0);
%R-squared for Tykcpr, Ce=0, 25: 0.7537


%Resistant, 75%, Ce=0.66
Tykcpr_75_Ce066_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_coculture/R_sim_i4_Ce0.66_sharedK_m14_10days.csv');
Tykcpr_75_Ce066_predicted = table2array(Tykcpr_75_Ce066_predicted(:,2)); %only the second column 
Tykcpr_75_Ce066_predicted = repmat(Tykcpr_75_Ce066_predicted,3,1); %repeat 3 times 

SS_res_Tykcpr_75_Ce066 = sum((Tykcpr_75_Ce066_all - Tykcpr_75_Ce066_predicted).^2);
SS_tot_Tykcpr_75_Ce066 = sum((Tykcpr_75_Ce066_all - mean(Tykcpr_75_Ce066_all)).^2);
R2_Tykcpr_75_Ce066 = 1 - SS_res_Tykcpr_75_Ce066 / SS_tot_Tykcpr_75_Ce066;

fprintf('R-squared for Tykcpr, Ce=0.66, 75: %.4f\n', R2_Tykcpr_75_Ce066);
%R-squared for Tykcpr, Ce=0.66, 75: 0.1451

%Resistant, 50%, Ce=0.66
Tykcpr_50_Ce066_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_coculture/R_sim_i5_Ce0.66_sharedK_m14_10days.csv');
Tykcpr_50_Ce066_predicted = table2array(Tykcpr_50_Ce066_predicted(:,2)); %only the second column 
Tykcpr_50_Ce066_predicted = repmat(Tykcpr_50_Ce066_predicted,3,1); %repeat 3 times 

SS_res_Tykcpr_50_Ce066 = sum((Tykcpr_50_Ce066_all - Tykcpr_50_Ce066_predicted).^2);
SS_tot_Tykcpr_50_Ce066 = sum((Tykcpr_50_Ce066_all - mean(Tykcpr_50_Ce066_all)).^2);
R2_Tykcpr_50_Ce066 = 1 - SS_res_Tykcpr_50_Ce066 / SS_tot_Tykcpr_50_Ce066;

fprintf('R-squared for Tykcpr, Ce=0.66, 50: %.4f\n', R2_Tykcpr_50_Ce066);
%R-squared for Tykcpr, Ce=0.66, 50: 0.5846

%Resistant, 25%, Ce=0.66
Tykcpr_25_Ce066_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_coculture/R_sim_i6_Ce0.66_sharedK_m14_10days.csv');
Tykcpr_25_Ce066_predicted = table2array(Tykcpr_25_Ce066_predicted(:,2)); %only the second column 
Tykcpr_25_Ce066_predicted = repmat(Tykcpr_25_Ce066_predicted,3,1); %repeat 3 times 

SS_res_Tykcpr_25_Ce066 = sum((Tykcpr_25_Ce066_all - Tykcpr_25_Ce066_predicted).^2);
SS_tot_Tykcpr_25_Ce066 = sum((Tykcpr_25_Ce066_all - mean(Tykcpr_25_Ce066_all)).^2);
R2_Tykcpr_25_Ce066 = 1 - SS_res_Tykcpr_25_Ce066 / SS_tot_Tykcpr_25_Ce066;

fprintf('R-squared for Tykcpr, Ce=0.66, 25: %.4f\n', R2_Tykcpr_25_Ce066);
%R-squared for Tykcpr, Ce=0.66, 25: 0.5298
