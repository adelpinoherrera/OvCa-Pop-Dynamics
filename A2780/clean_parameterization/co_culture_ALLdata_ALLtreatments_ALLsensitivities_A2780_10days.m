%ALL DATA
%RE-PARAMETERIZE co-culture equations using all the data and
%re-parameterize only rS, K, alphaS and rR (most sensitive parameters)
%days, 4 parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Import co-culture data and format it 
clear all;
outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture_alldata';

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


%Import mono-culture data and format it 
%%%Sensitive 
monoS_files = {'/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_Untreated_30k_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC50_1um_30k_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC75_0.62um_30k_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_IC25_1.47um_30k_cellsinwell.csv'};
ce_mono_S = [0.00, 1.00, 0.62, 1.47];
monoS_labels = {'Ce=1.00', 'Ce=0.62', 'Ce=1.47'};

monoS_data = cell(4,1);
S_mono_all = cell(4,1);

for i = 1:4
    monoS_data{i} = readtable(monoS_files{i});
    S_mono_all{i} = [monoS_data{i}.Rep1, monoS_data{i}.Rep2, monoS_data{i}.Rep3];
end
t_mono_S = monoS_data{1}.Day;
idx_0_10_mono = (t_mono_S >= 0) & (t_mono_S <= 10);
t_mono_S_10 = t_mono_S(idx_0_10_mono);

S_mono_10_all = cell(4,1);
for i = 1:4
    S_mono_10_all{i} = S_mono_all{i}(idx_0_10_mono, :);
end

%%%Resistant 
monoR_files = {'/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780cis_Untreated_30k_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780cis_IC50_1um_30k_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780cis_IC75_0.62um_30k_cellsinwell.csv',
    '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780cis_IC25_1.47um_30k_cellsinwell.csv'};

ce_mono_R = [0.0, 1.00, 0.62, 1.47];

monoR_data = cell(4,1);
R_mono_all = cell(4,1);

for i = 1:4
    monoR_data{i} = readtable(monoR_files{i});
    R_mono_all{i} = [monoR_data{i}.Rep1, monoR_data{i}.Rep2, monoR_data{i}.Rep3];
end
t_mono_R = monoR_data{1}.Day;
idx_0_10_R = (t_mono_R >= 0) & (t_mono_R <= 10);
t_mono_R_10 = t_mono_R(idx_0_10_R);

R_mono_10_all = cell(4,1);
for i = 1:4
    R_mono_10_all{i} = R_mono_all{i}(idx_0_10_R, :);
end

%%%Fixed parameter from estimating with only shared K 
nS = 3.7779;
IC50S = 0.9960177;

alphaR = 1.7982;
nR = 2.5166;
IC50R = 8.338535;

%Estimate for parameters rS, K, alphaS, rR
p0 = [0.8, 1100000.0, 0.1, ... %rS, K, alphaS
      0.5]; %rR
lb = [0.001, 500000.0, 0.001, ...
      0.0001];
ub = [20.0, 3000000.0, 20.0, ...
      20.0];

res_fun_joint = @(p) joint_residuals_coculture_mono_sens2_m13( ...
    p, nS, IC50S, alphaR, nR, IC50R, ...
    t_coculture_10, S_raw_10_all, R_raw_10_all, ce_coculture, ...
    t_mono_S_10, S_mono_10_all, ce_mono_S, ...
    t_mono_R_10, R_mono_10_all, ce_mono_R);

problem_joint = createOptimProblem('lsqnonlin', ...
    'x0', p0, ...
    'objective', res_fun_joint, ...
    'lb', lb, ...
    'ub', ub);

ms_joint = MultiStart('UseParallel', true, 'Display', 'iter');
nStarts = 200;

[hat_joint, resnorm_joint] = run(ms_joint, problem_joint, nStarts);
fprintf('Joint fit (co-culture + mono-culture, 10 days)\n');
fprintf('Fitted rS = %.4f\n', hat_joint(1));
fprintf('Fitted K  = %.4f\n', hat_joint(2));
fprintf('Fitted alphaS = %.4f\n', hat_joint(3));
fprintf('Fitted rR = %.4f\n', hat_joint(4));
fprintf('Sum of squared residuals = %.2f\n', resnorm_joint);
% 71 out of 200 local solver runs converged with a positive local solver exitflag.
% Joint fit (co-culture + mono-culture, 10 days)
% Fitted rS = 0.7719
% Fitted K  = 1112765.1958
% Fitted alphaS = 0.0911
% Fitted rR = 0.6319

% 73 out of 200 local solver runs converged with a positive local solver exitflag.
% Joint fit (co-culture + mono-culture, 10 days)
% Fitted rS = 0.7719
% Fitted K  = 1112822.2016
% Fitted alphaS = 0.0911
% Fitted rR = 0.6319
% Sum of squared residuals = 7307519906411.63

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(res_fun_joint, hat_joint, lb, ub);

% Compute 95% confidence intervals
ci = nlparci(p_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95 Confidence Intervals:\n')
fprintf('rS: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
fprintf('K: [%.4f, %.4f]\n', ci(2,1), ci(2,2));
fprintf('alphaS: [%.4f, %.4f]\n', ci(3,1), ci(3,2));
fprintf('rR: [%.4f, %.4f]\n', ci(4,1), ci(4,2));
% Solver stopped prematurely.
% lsqnonlin stopped because it exceeded the function evaluation limit,
% options.MaxFunctionEvaluations = 4.000000e+02.
% 95 Confidence Intervals:
% rS: [0.7481, 0.7957] (0.7957-0.7481)/0.7719 = 0.0617
% K: [1075110.8286, 1150533.5746]
% (1150533.5746-1075110.8286)/1112822.2016 = 0.0678
% alphaS: [0.0838, 0.0983] (0.0983-0.0838)/0.0911 = 0.1592
% rR: [0.6152, 0.6487] (0.6487-0.6152)/0.6319 = 0.053

[~, BIC_info] = joint_residuals_coculture_mono_sens2_m13( ...
    hat_joint, nS, IC50S, alphaR, nR, IC50R, ...
    t_coculture_10, S_raw_10_all, R_raw_10_all, ce_coculture, ...
    t_mono_S_10, S_mono_10_all, ce_mono_S, ...
    t_mono_R_10, R_mono_10_all, ce_mono_R);

fprintf('nBIC summary:\n');
fprintf('  Co-culture = %.2f\n', BIC_info.co);
fprintf('  Mono S     = %.2f\n', BIC_info.monoS);
fprintf('  Mono R     = %.2f\n', BIC_info.monoR);
fprintf('  Total      = %.2f\n', BIC_info.total);

% nBIC summary:
%   Co-culture = 9231.39
%   Mono S     = 3076.51
%   Mono R     = 2999.10
%   Total      = 15290.24

%co-culture plot
figure;
for i = 1:6
    S_raw = S_raw_10_all{i};
    R_raw = R_raw_10_all{i};
    [nTimes, nReps] = size(S_raw);

    S_mod_all = zeros(nTimes, nReps);
    R_mod_all = zeros(nTimes, nReps);

    for j = 1:nReps
        S0 = S_raw(1,j);
        R0 = R_raw(1,j);

        [S_mod, R_mod] = simulate_coculture_all_sens2_m13( ...
            t_coculture_10, S0, R0, ...
            hat_joint(1), hat_joint(2), hat_joint(3), nS, IC50S, ...
            hat_joint(4), alphaR, nR, IC50R, ce_coculture(i));

        S_mod_all(:,j) = S_mod;
        R_mod_all(:,j) = R_mod;
    end

    S_mean = mean(S_raw, 2);
    S_sd   = std(S_raw, 0, 2);
    R_mean = mean(R_raw, 2);
    R_sd   = std(R_raw, 0, 2);

    S_mod_mean = mean(S_mod_all, 2);
    R_mod_mean = mean(R_mod_all, 2);
    Tot_mean = S_mean + R_mean;
    Tot_mod_mean = S_mod_mean + R_mod_mean;

    writematrix([t_coculture_10, S_mod_all], ...
    fullfile(outdir, sprintf('S_sim_i%d_Ce%.2f_top4param_m13_10days.csv', i, ce_coculture(i))));

writematrix([t_coculture_10, R_mod_all], ...
    fullfile(outdir, sprintf('R_sim_i%d_Ce%.2f_top4param_m13_10days.csv', i, ce_coculture(i))));

    subplot(2,3,i); hold on;

    Sdata = errorbar(t_coculture_10, S_mean, S_sd, 'r.', 'MarkerSize', 12, 'LineStyle', 'none');
    Rdata = errorbar(t_coculture_10, R_mean, R_sd, 'g.', 'MarkerSize', 12, 'LineStyle', 'none');
    Tdata = plot(t_coculture_10, Tot_mean, 'b.', 'MarkerSize', 12);

    Smodel = plot(t_coculture_10, S_mod_mean, 'r-', 'LineWidth', 2);
    Rmodel = plot(t_coculture_10, R_mod_mean, 'g-', 'LineWidth', 2);
    Tmodel = plot(t_coculture_10, Tot_mod_mean, 'b-', 'LineWidth', 2);

    title(sprintf('Ce=%.2f, %s', ce_coculture(i), ratio_labels{mod(i-1,3)+1}));
    xlabel('Time (days)');
    ylabel('Cell count');
    ylim([0 1400000]);
    grid on;
    legend([Sdata(1), Rdata(1), Tdata(1), Smodel(1), Rmodel(1), Tmodel(1)], ...
        {'S data','R data','Total data','S model','R model','Total model'}, ...
        'Location','best');
end
sgtitle(sprintf('rS = %.3f, K = %.3f, alphaS = %.3f, rR = %.3f', ...
    hat_joint(1), hat_joint(2), hat_joint(3), hat_joint(4)));

saveas(gcf, fullfile(outdir, 'co-culture_A2780_ALLdata_ALLtreat_ALLsensitivities_sens2_m13_constD_10days.png'));

%sensitive plot
figure;
for i = 1:4
    subplot(2,2,i); hold on;

    S_raw = S_mono_10_all{i};
    [nTimes, nReps] = size(S_raw);

    plot(t_mono_S_10, S_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);

    S_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        S0 = S_raw(1,j);
        S_fit_all(:,j) = simulate_mono_sens2_m13( ...
            t_mono_S_10, S0, hat_joint(1), hat_joint(2), hat_joint(3), nS, ...
            IC50S, ce_mono_S(i));
    end

    S_fit_mean = mean(S_fit_all, 2);
    plot(t_mono_S_10, S_fit_mean, 'r-', 'LineWidth', 3);

    title(sprintf('Sensitive mono, Ce = %.2f', ce_mono_S(i)));
    xlabel('Time (days)');
    ylabel('Sensitive cells');
    ylim([0 1400000]);
    grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Joint fit: rS = %.3f, K = %.3f, alphaS = %.3f, nS = %.3f', ...
    hat_joint(1), hat_joint(2), hat_joint(3), nS));

saveas(gcf, fullfile(outdir, 'co-culture_A2780_ALLdata_ALLtreat_ALLsens_ONLYSensitive_sens2_m13_constD_10days.png'));

%Resistant plot 
figure;
for i = 1:4
    subplot(2,2,i); hold on;

    R_raw = R_mono_10_all{i};
    [nTimes, nReps] = size(R_raw);

    plot(t_mono_R_10, R_raw, 'ko', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 6);

    R_fit_all = zeros(nTimes, nReps);
    for j = 1:nReps
        R0 = R_raw(1,j);
        R_fit_all(:,j) = simulate_mono_res_sens2_m13( ...
            t_mono_R_10, R0, hat_joint(4), hat_joint(2), alphaR, nR,...
            IC50R, ce_mono_R(i));
    end

    R_fit_mean = mean(R_fit_all, 2);
    plot(t_mono_R_10, R_fit_mean, 'g-', 'LineWidth', 3);

    title(sprintf('Resistant mono, Ce = %.2f', ce_mono_R(i)));
    xlabel('Time (days)');
    ylabel('Resistant cells');
    ylim([0 1400000]);
    grid on;
    legend('Replicates', 'Model fit', 'Location', 'best');
end
sgtitle(sprintf('Joint fit: rR = %.3f, K = %.3f, alphaR = %.3f, nR = %.3f', ...
    hat_joint(4), hat_joint(2), alphaR, nR));

saveas(gcf, fullfile(outdir, 'co-culture_A2780_ALLdata_ALLtreat_ALLsens_ONLYResistant_sens2_m13_constD_10days.png'));

%R squared values
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
A2780_75_Ce0_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture_alldata/S_sim_i1_Ce0.00_top4param_m13_10days.csv');
A2780_75_Ce0_predicted = table2array(A2780_75_Ce0_predicted(:,2)); %only the second column 
A2780_75_Ce0_predicted = repmat(A2780_75_Ce0_predicted,3,1); %repeat 3 times 

SS_res_A2780_75_Ce0 = sum((A2780_75_Ce0_all - A2780_75_Ce0_predicted).^2);
SS_tot_A2780_75_Ce0 = sum((A2780_75_Ce0_all - mean(A2780_75_Ce0_all)).^2);
R2_A2780_75_Ce0 = 1 - SS_res_A2780_75_Ce0 / SS_tot_A2780_75_Ce0;

fprintf('R-squared for A2780, Ce=0, 75: %.4f\n', R2_A2780_75_Ce0);
%R-squared for A2780, Ce=0, 75: 0.6586

%Sensitive, 50%, Ce=0.0
A2780_50_Ce0_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture_alldata/S_sim_i2_Ce0.00_top4param_m13_10days.csv');
A2780_50_Ce0_predicted = table2array(A2780_50_Ce0_predicted(:,2)); %only the second column 
A2780_50_Ce0_predicted = repmat(A2780_50_Ce0_predicted,3,1); %repeat 3 times 

SS_res_A2780_50_Ce0 = sum((A2780_50_Ce0_all - A2780_50_Ce0_predicted).^2);
SS_tot_A2780_50_Ce0 = sum((A2780_50_Ce0_all - mean(A2780_50_Ce0_all)).^2);
R2_A2780_50_Ce0 = 1 - SS_res_A2780_50_Ce0 / SS_tot_A2780_50_Ce0;

fprintf('R-squared for A2780, Ce=0, 50: %.4f\n', R2_A2780_50_Ce0);
%R-squared for A2780, Ce=0, 50: 0.9063

%Sensitive, 25%, Ce=0.0
A2780_25_Ce0_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture_alldata/S_sim_i3_Ce0.00_top4param_m13_10days.csv');
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
%R-squared for A2780, Ce=0, 25: 0.9025


%Sensitive, 75%, Ce=1.00
A2780_75_Ce1_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture_alldata/S_sim_i4_Ce1.00_top4param_m13_10days.csv');
A2780_75_Ce1_predicted = table2array(A2780_75_Ce1_predicted(:,2)); %only the second column 
A2780_75_Ce1_predicted = repmat(A2780_75_Ce1_predicted,3,1); %repeat 3 times 

SS_res_A2780_75_Ce1 = sum((A2780_75_Ce1_all - A2780_75_Ce1_predicted).^2);
SS_tot_A2780_75_Ce1 = sum((A2780_75_Ce1_all - mean(A2780_75_Ce1_all)).^2);
R2_A2780_75_Ce1 = 1 - SS_res_A2780_75_Ce1 / SS_tot_A2780_75_Ce1;

fprintf('R-squared for A2780, Ce=1.00, 75: %.4f\n', R2_A2780_75_Ce1);
%R-squared for A2780, Ce=0.66, 75: 0.4956

%Sensitive, 50%, Ce=1.00
A2780_50_Ce1_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture_alldata/S_sim_i5_Ce1.00_top4param_m13_10days.csv');
A2780_50_Ce1_predicted = table2array(A2780_50_Ce1_predicted(:,2)); %only the second column 
A2780_50_Ce1_predicted = repmat(A2780_50_Ce1_predicted,3,1); %repeat 3 times 

SS_res_A2780_50_Ce1 = sum((A2780_50_Ce1_all - A2780_50_Ce1_predicted).^2);
SS_tot_A2780_50_Ce1 = sum((A2780_50_Ce1_all - mean(A2780_50_Ce1_all)).^2);
R2_A2780_50_Ce1 = 1 - SS_res_A2780_50_Ce1 / SS_tot_A2780_50_Ce1;

fprintf('R-squared for A2780, Ce=1.00, 50: %.4f\n', R2_A2780_50_Ce1);
%R-squared for A2780, Ce=0.66, 50: 0.5665

%Sensitive, 25%, Ce=1.00
A2780_25_Ce1_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture_alldata/S_sim_i6_Ce1.00_top4param_m13_10days.csv');
A2780_25_Ce1_predicted = table2array(A2780_25_Ce1_predicted(:,2)); %only the second column 
A2780_25_Ce1_predicted = repmat(A2780_25_Ce1_predicted,3,1); %repeat 3 times 

SS_res_A2780_25_Ce1 = sum((A2780_25_Ce1_all - A2780_25_Ce1_predicted).^2);
SS_tot_A2780_25_Ce1 = sum((A2780_25_Ce1_all - mean(A2780_25_Ce1_all)).^2);
R2_A2780_25_Ce1 = 1 - SS_res_A2780_25_Ce1 / SS_tot_A2780_25_Ce1;

fprintf('R-squared for A2780, Ce=1.00, 25: %.4f\n', R2_A2780_25_Ce1);
%R-squared for Tyk, Ce=0.66, 25: 0.4521


%Resistant, 75%, Ce=0.0
A2780cis_75_Ce0_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture_alldata/R_sim_i1_Ce0.00_top4param_m13_10days.csv');
A2780cis_75_Ce0_predicted = table2array(A2780cis_75_Ce0_predicted(:,2)); %only the second column 
A2780cis_75_Ce0_predicted = repmat(A2780cis_75_Ce0_predicted,3,1); %repeat 3 times 

SS_res_A2780cis_75_Ce0 = sum((A2780cis_75_Ce0_all - A2780cis_75_Ce0_predicted).^2);
SS_tot_A2780cis_75_Ce0 = sum((A2780cis_75_Ce0_all - mean(A2780cis_75_Ce0_all)).^2);
R2_A2780cis_75_Ce0 = 1 - SS_res_A2780cis_75_Ce0 / SS_tot_A2780cis_75_Ce0;

fprintf('R-squared for A2780cis, Ce=0, 75: %.4f\n', R2_A2780cis_75_Ce0);
%R-squared for A2780cis, Ce=0, 75: -0.9905

%Resistant, 50%, Ce=0.0
A2780cis_50_Ce0_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture_alldata/R_sim_i2_Ce0.00_top4param_m13_10days.csv');
A2780cis_50_Ce0_predicted = table2array(A2780cis_50_Ce0_predicted(:,2)); %only the second column 
A2780cis_50_Ce0_predicted = repmat(A2780cis_50_Ce0_predicted,3,1); %repeat 3 times 

SS_res_A2780cis_50_Ce0 = sum((A2780cis_50_Ce0_all - A2780cis_50_Ce0_predicted).^2);
SS_tot_A2780cis_50_Ce0 = sum((A2780cis_50_Ce0_all - mean(A2780cis_50_Ce0_all)).^2);
R2_A2780cis_50_Ce0 = 1 - SS_res_A2780cis_50_Ce0 / SS_tot_A2780cis_50_Ce0;

fprintf('R-squared for A2780cis, Ce=0, 50: %.4f\n', R2_A2780cis_50_Ce0);
%R-squared for A2780cis, Ce=0, 50: -0.5023

%Resistant, 25%, Ce=0.0
A2780cis_25_Ce0_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture_alldata/R_sim_i3_Ce0.00_top4param_m13_10days.csv');
A2780cis_25_Ce0_predicted = table2array(A2780cis_25_Ce0_predicted(:,2)); %only the second column 
A2780cis_25_Ce0_predicted = repmat(A2780cis_25_Ce0_predicted,3,1); %repeat 3 times 

SS_res_A2780cis_25_Ce0 = sum((A2780cis_25_Ce0_all - A2780cis_25_Ce0_predicted).^2);
SS_tot_A2780cis_25_Ce0 = sum((A2780cis_25_Ce0_all - mean(A2780cis_25_Ce0_all)).^2);
R2_A2780cis_25_Ce0 = 1 - SS_res_A2780cis_25_Ce0 / SS_tot_A2780cis_25_Ce0;

fprintf('R-squared for A2780cis, Ce=0, 25: %.4f\n', R2_A2780cis_25_Ce0);
%R-squared for A2780cis, Ce=0, 25: 0.9044


%Resistant, 75%, Ce=1.00
A2780cis_75_Ce1_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture_alldata/R_sim_i4_Ce1.00_top4param_m13_10days.csv');
A2780cis_75_Ce1_predicted = table2array(A2780cis_75_Ce1_predicted(:,2)); %only the second column 
A2780cis_75_Ce1_predicted = repmat(A2780cis_75_Ce1_predicted,3,1); %repeat 3 times 

SS_res_A2780cis_75_Ce1 = sum((A2780cis_75_Ce1_all - A2780cis_75_Ce1_predicted).^2);
SS_tot_A2780cis_75_Ce1 = sum((A2780cis_75_Ce1_all - mean(A2780cis_75_Ce1_all)).^2);
R2_A2780cis_75_Ce1 = 1 - SS_res_A2780cis_75_Ce1 / SS_tot_A2780cis_75_Ce1;

fprintf('R-squared for A2780cis, Ce=1.00, 75: %.4f\n', R2_A2780cis_75_Ce1);
%R-squared for A2780cis, Ce=0.66, 75: 0.8064

%Resistant, 50%, Ce=1.00
A2780cis_50_Ce1_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture_alldata/R_sim_i5_Ce1.00_top4param_m13_10days.csv');
A2780cis_50_Ce1_predicted = table2array(A2780cis_50_Ce1_predicted(:,2)); %only the second column 
A2780cis_50_Ce1_predicted = repmat(A2780cis_50_Ce1_predicted,3,1); %repeat 3 times 

SS_res_A2780cis_50_Ce1 = sum((A2780cis_50_Ce1_all - A2780cis_50_Ce1_predicted).^2);
SS_tot_A2780cis_50_Ce1 = sum((A2780cis_50_Ce1_all - mean(A2780cis_50_Ce1_all)).^2);
R2_A2780cis_50_Ce1 = 1 - SS_res_A2780cis_50_Ce1 / SS_tot_A2780cis_50_Ce1;

fprintf('R-squared for A2780cis, Ce=1.00, 50: %.4f\n', R2_A2780cis_50_Ce1);
%R-squared for A2780cis, Ce=1.00, 50: 0.8352

%Resistant, 25%, Ce=1.00
A2780cis_25_Ce1_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_coculture_alldata/R_sim_i6_Ce1.00_top4param_m13_10days.csv');
A2780cis_25_Ce1_predicted = table2array(A2780cis_25_Ce1_predicted(:,2)); %only the second column 
A2780cis_25_Ce1_predicted = repmat(A2780cis_25_Ce1_predicted,3,1); %repeat 3 times 

SS_res_A2780cis_25_Ce1 = sum((A2780cis_25_Ce1_all - A2780cis_25_Ce1_predicted).^2);
SS_tot_A2780cis_25_Ce1 = sum((A2780cis_25_Ce1_all - mean(A2780cis_25_Ce1_all)).^2);
R2_A2780cis_25_Ce1 = 1 - SS_res_A2780cis_25_Ce1 / SS_tot_A2780cis_25_Ce1;

fprintf('R-squared for A2780cis, Ce=1.00, 25: %.4f\n', R2_A2780cis_25_Ce1);
%R-squared for Tykcpr, Ce=0.66, 25: 0.7638