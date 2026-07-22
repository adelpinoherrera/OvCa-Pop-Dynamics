%%Trying to estimate parameters with MATLAB for untreated cells
%Compare logistic and generalized logistic models 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Repeat everything for 10 days and lsqnonlin, and see how the parameters
%%%change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%need to change it a little bit to be able to use the lsqnonlin 
clear all;
outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_untreated';
%%%%Sensitive cells first 
Tyk_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_Untreated_40k_cellsinwell_new.csv');
Tyk_t_days = Tyk_data.Day;
Tyk_10t_days = Tyk_t_days(1:11);
Tyk_y_raw = [Tyk_data.Rep1, Tyk_data.Rep2, Tyk_data.Rep3];
Tyk_10y_raw = Tyk_y_raw(1:11,:);
[nDays, nReps] = size(Tyk_10y_raw);

%Vectorize time and data so all replicates are used, no need to do this to
%use lsqnonlin
Tyk_10t_all = repmat(Tyk_10t_days, nReps, 1); %should be 45 x 1 
Tyk_10y_all = Tyk_10y_raw(:); %should be 45 x 1 
Tyk_10N = numel(Tyk_10y_all); %total observations so 45


%Define objective function for both models, need to define a residual
%function to be able to use lsqnonlin
logistic_res = @(p, t_dummy) logistic_ode_model(p, Tyk_10t_days, Tyk_10y_raw)-Tyk_10y_all;

genlog_res   = @(p, t_dummy) genlog_ode_model(p, Tyk_10t_days, Tyk_10y_raw)-Tyk_10y_all;

%Define initial guesses for logistic
K0 = 1.4e6; %max(y_all); but it could be higher than the last day counts
r0 = 0.5; %pretty small

p0_log = [K0, r0];
lb_log = [0,  0];     % K >= 0, r >= 0
ub_log = [Inf, Inf]; %could define a differen upper bound

% Create lsqcurvefit problem for logistic
problem_log = createOptimProblem( ...
    'lsqnonlin', ... %eve though we are using a multistart, we still need to define the function to minimize
    'x0',        p0_log, ...
    'objective', logistic_res, ...
    'lb',        lb_log, ...
    'ub',        ub_log);

ms = MultiStart('UseParallel', true);  % set false if no Parallel Toolbox
nStarts = 200; %could do more starts if we wanted to 

[p_log, resnorm_log] = run(ms, problem_log, nStarts); %solve the problem with multistart 
fprintf('Fitted logistic parameters for Tyk (10 days)\n')
fprintf('Fitted growth rate r = %.4f\n', p_log(2));
fprintf('Fitted carrying capacity k = %.4f\n', p_log(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_log);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted logistic parameters for Tyk (10 days)
% Fitted growth rate r = 0.6433
% Fitted carrying capacity k = 545681.7362
% Sum of squared residuals = 18119492660.41

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_log_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(logistic_res, p_log, lb_log, ub_log);

% Compute 95% confidence intervals
ci = nlparci(p_log_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95% Confidence Intervals:\n');
fprintf('r: [%.4f, %.4f]\n', ci(2,1), ci(2,2));
fprintf('K: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
% Local minimum possible.
% lsqnonlin stopped because the final change in the sum of squares relative to 
% its initial value is less than the value of the function tolerance.
% <stopping criteria details>
% 95r: [0.6178, 0.6688] (0.6688-0.6178)/0.6433 = 0.079279
% K: [512235.2248, 579129.1008] (579129.1008-512235.2248)/545681.7362 =
% 0.122588

%Define initial guesses for generalized logistic
K0  = 1.4e6;         %max(y_all);
r0  = 0.5;
nu0 = 1;            % nu=1 reduces to logistic

p0_gen = [K0, r0, nu0];
lb_gen = [0,  0,   0.1];   % nu > 0
ub_gen = [Inf, Inf, 30];

problem_gen = createOptimProblem( ...
    'lsqnonlin', ...
    'x0',        p0_gen, ...
    'objective', genlog_res, ...
    'lb',        lb_gen, ...
    'ub',        ub_gen);

[p_gen, resnorm_gen] = run(ms, problem_gen, nStarts);
fprintf('Fitted generalized logistic parameters for Tyk (10 days)\n')
fprintf('Fitted growth rate r = %.4f\n', p_gen(2));
fprintf('Fitted carrying capacity k = %.4f\n', p_gen(1));
fprintf('Fitted n parameter n = %.4f\n', p_gen(3));
fprintf('Sum of squared residuals = %.2f\n', resnorm_gen);
% 199 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted generalized logistic parameters for Tyk (10 days)
% Fitted growth rate r = 0.5938
% Fitted carrying capacity k = 516423.7758
% Fitted n parameter n = 1.4037
% Sum of squared residuals = 17332788874.74

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_gen_final, resnorm_gen_final, residual_gen_final, exitflag_gen, output_gen, lambda_gen, J_gen] = ...
lsqnonlin(genlog_res, p_gen, lb_gen, ub_gen);

% Compute 95% confidence intervals
ci_gen = nlparci(p_gen_final, residual_gen_final, 'jacobian', J_gen);

% Print confidence intervals
fprintf('95% Confidence Intervals:\n');
fprintf('r: [%.4f, %.4f]\n', ci_gen(2,1), ci_gen(2,2));
fprintf('K: [%.4f, %.4f]\n', ci_gen(1,1), ci_gen(1,2));
fprintf('nu: [%.4f, %.4f]\n', ci_gen(3,1), ci_gen(3,2)); 
% Local minimum possible.
% lsqnonlin stopped because the final change in the sum of squares relative to 
% its initial value is less than the value of the function tolerance.
% <stopping criteria details>
% 95r: [0.5285, 0.6591] (0.6591-0.5285)/0.5938 = 0.2199
% K: [467927.6554, 564919.8963] (564919.8963-467927.6554)/516423.7758 =
% 0.1878
% nu: [0.6465, 2.1607] (2.1607-0.6465)/1.4037 = 1.0787

%%Plot BIC for the models 
N = numel(Tyk_10y_all);  
k_log = length(p_log);
k_gen = length(p_gen);
SSE_log = resnorm_log; %here the resnorm is already squared
SSE_gen = resnorm_gen;

% BIC formula: N*ln(SSE/N) + k*ln(N)
BIC_log = N * log(SSE_log / N) + k_log * log(N);
BIC_gen = N * log(SSE_gen / N) + k_gen * log(N);

% Display results
fprintf('Logistic:     BIC = %.3f\n', BIC_log);
fprintf('Gen Logistic: BIC = %.3f\n', BIC_gen);

if BIC_log < BIC_gen
    fprintf('Logistic preferred (lower BIC)\n');
else
    fprintf('Gen Logistic preferred (lower BIC)\n');
end

% Logistic:     BIC = 671.077
% Gen Logistic: BIC = 673.108
% Logistic preferred (lower BIC)

%%Plot the curves
% Dense time grid for plotting continuous curves
t_fine = linspace(0, 14, 200)';

% Predict for each replicate separately at t_fine, using fitted params
% Example: use the same ODE solvers but with t_fine instead of 1–14
% (for simplicity, here we interpolate from 1–14)

y_log_fit_full = logistic_ode_model(p_log, Tyk_10t_days, Tyk_10y_raw);
y_gen_fit_full = genlog_ode_model(p_gen, Tyk_10t_days, Tyk_10y_raw);

% Reshape back to 14×3 for plotting
y_log_fit = reshape(y_log_fit_full, [nDays, nReps]);
y_gen_fit = reshape(y_gen_fit_full, [nDays, nReps]);

figure;
hold on;
% Plot data
plot(Tyk_10t_days, Tyk_10y_raw, 'ko', 'MarkerFaceColor',[0.8 0.8 0.8]);

% Plot fitted curves per replicate (logistic)
for j = 1:nReps
    plot(Tyk_10t_days, y_log_fit(:,j), 'b-', 'LineWidth', 1.5);
end

% Plot generalized logistic as dashed red (optional)
for j = 1:nReps
    plot(Tyk_10t_days, y_gen_fit(:,j), 'r--', 'LineWidth', 1.5);
end

xlabel('Time (days)');
ylabel('Cell count (or density)');
legend('Data','Logistic fit','Gen. logistic fit','Location','best');
grid on;
hold off;

saveas(gcf, fullfile(outdir,'Tyk_untreated_logandgenlog_10days_new.png'));

writematrix(y_log_fit, fullfile(outdir,'Tyk_untreated_log_10days_new.csv'));
writematrix(y_gen_fit, fullfile(outdir,'Tyk_untreated_genlog_10days_new.csv'));

%calculate Rsquared values 
clear all;
%%%%Sensitive cells first 
Tyk_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tyk_Untreated_40k_cellsinwell_new.csv');
Tyk_t_days = Tyk_data.Day;
Tyk_y_raw = [Tyk_data.Rep1, Tyk_data.Rep2, Tyk_data.Rep3];
Tyk_10y_raw = Tyk_y_raw(1:11,:);
Tyk_10y_all = Tyk_10y_raw(:); %should be 33 x 1 


Tyk_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_untreated/Tyk_untreated_log_10days_new.csv');
Tyk_predicted = table2array(Tyk_predicted);
Tyk_predicted = Tyk_predicted(:); %33x1

SS_res = sum((Tyk_10y_all - Tyk_predicted).^2);
SS_tot = sum((Tyk_10y_all - mean(Tyk_10y_all)).^2);
R2 = 1 - SS_res / SS_tot;

fprintf('R-squared for logistic fit: %.4f\n', R2);
%R-squared for logistic fit: 0.9833

Tyk_predicted_gen = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tyk_untreated/Tyk_untreated_genlog_10days_new.csv');
Tyk_predicted_gen = table2array(Tyk_predicted_gen);
Tyk_predicted_gen = Tyk_predicted_gen(:);

SS_res = sum((Tyk_10y_all - Tyk_predicted_gen).^2);
SS_tot = sum((Tyk_10y_all - mean(Tyk_10y_all)).^2);
R2 = 1 - SS_res / SS_tot;

fprintf('R-squared for generalized fit: %.4f\n', R2);
%R-squared for generalized fit: 0.9840

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Repeat the 10 day process for Tykcpr with lsqnonlin
clear all;
outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tykcpr_untreated';

%%%%Resistant cells cells  
Tykcpr_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tykcpr_Untreated_40k_cellsinwell_new.csv');
Tykcpr_t_days = Tykcpr_data.Day;
Tykcpr_10t_days = Tykcpr_t_days(1:11);
Tykcpr_y_raw = [Tykcpr_data.Rep1, Tykcpr_data.Rep2, Tykcpr_data.Rep3];
Tykcpr_10y_raw = Tykcpr_y_raw(1:11,:);
[nDays, nReps] = size(Tykcpr_10y_raw);

%Vectorize time and data so all replicates are used, no need to do this to
%use lsqnonlin
Tykcpr_10t_all = repmat(Tykcpr_10t_days, nReps, 1); %should be 45 x 1 
Tykcpr_10y_all = Tykcpr_10y_raw(:); %should be 45 x 1 
Tykcpr_10N = numel(Tykcpr_10y_all); %total observations so 45

%Define objective function for both models, need to define a residual
%function to be able to use lsqnonlin
logistic_res = @(p, t_dummy) logistic_ode_model(p, Tykcpr_10t_days, Tykcpr_10y_raw)-Tykcpr_10y_all;

genlog_res   = @(p, t_dummy) genlog_ode_model(p, Tykcpr_10t_days, Tykcpr_10y_raw)-Tykcpr_10y_all;

%Define initial guesses for logistic
K0 = 1.4e6; %max(y_all); but it could be higher than the last day counts
r0 = 0.5; %pretty small

p0_log = [K0, r0];
lb_log = [0,  0];     % K >= 0, r >= 0
ub_log = [Inf, Inf]; %could define a differen upper bound

% Create lsqcurvefit problem for logistic
problem_log = createOptimProblem( ...
    'lsqnonlin', ... %eve though we are using a multistart, we still need to define the function to minimize
    'x0',        p0_log, ...
    'objective', logistic_res, ...
    'lb',        lb_log, ...
    'ub',        ub_log);

ms = MultiStart('UseParallel', true);  % set false if no Parallel Toolbox
nStarts = 200; %could do more starts if we wanted to 

[p_log, resnorm_log] = run(ms, problem_log, nStarts); %solve the problem with multistart 
fprintf('Fitted logistic parameters for Tykcpr (10 days)\n')
fprintf('Fitted growth rate r = %.4f\n', p_log(2));
fprintf('Fitted carrying capacity k = %.4f\n', p_log(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_log);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted logistic parameters for Tykcpr (10 days)
% Fitted growth rate r = 0.5339
% Fitted carrying capacity k = 421668.8619
% Sum of squared residuals = 13684837143.90

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_log_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(logistic_res, p_log, lb_log, ub_log);

% Compute 95% confidence intervals
ci = nlparci(p_log_final, residual_final, 'jacobian', J);

% Print confidence intervals
fprintf('95% Confidence Intervals:\n');
fprintf('r: [%.4f, %.4f]\n', ci(2,1), ci(2,2));
fprintf('K: [%.4f, %.4f]\n', ci(1,1), ci(1,2));
% Local minimum possible.
% lsqnonlin stopped because the final change in the sum of squares relative to 
% its initial value is less than the value of the function tolerance.
% <stopping criteria details>
% 95r: [0.5039, 0.5640] (0.5640-0.5039)/0.5339 = 0.112568
% K: [373990.7229, 469340.7711] (469340.7711-373990.7229)/421668.8619 =
% 0.226125

%Define initial guesses for generalized logistic
K0  = 1.4e6;         %max(y_all);
r0  = 0.5;
nu0 = 1;            % nu=1 reduces to logistic

p0_gen = [K0, r0, nu0];
lb_gen = [0,  0,   0.1];   % nu > 0
ub_gen = [Inf, Inf, 30];

problem_gen = createOptimProblem( ...
    'lsqnonlin', ...
    'x0',        p0_gen, ...
    'objective', genlog_res, ...
    'lb',        lb_gen, ...
    'ub',        ub_gen);

[p_gen, resnorm_gen] = run(ms, problem_gen, nStarts);
fprintf('Fitted generalized logistic parameters for Tykcpr (10 days)\n')
fprintf('Fitted growth rate r = %.4f\n', p_gen(2));
fprintf('Fitted carrying capacity k = %.4f\n', p_gen(1));
fprintf('Fitted n term n = %.4f\n', p_gen(3));
fprintf('Sum of squared residuals = %.2f\n', resnorm_gen);
% 197 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted generalized logistic parameters for Tykcpr (10 days)
% Fitted growth rate r = 0.4482
% Fitted carrying capacity k = 344144.3826
% Fitted n term n = 3.2848
% Sum of squared residuals = 9871504074.81

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_gen_final, resnorm_gen_final, residual_gen_final, exitflag_gen, output_gen, lambda_gen, J_gen] = ...
lsqnonlin(genlog_res, p_gen, lb_gen, ub_gen);

% Compute 95% confidence intervals
ci_gen = nlparci(p_gen_final, residual_gen_final, 'jacobian', J_gen);

% Print confidence intervals
fprintf('95% Confidence Intervals:\n');
fprintf('r: [%.4f, %.4f]\n', ci_gen(2,1), ci_gen(2,2));
fprintf('K: [%.4f, %.4f]\n', ci_gen(1,1), ci_gen(1,2));
fprintf('nu: [%.4f, %.4f]\n', ci_gen(3,1), ci_gen(3,2));
% Local minimum possible.
% lsqnonlin stopped because the final change in the sum of squares relative to 
% its initial value is less than the value of the function tolerance.
% <stopping criteria details>
% 95r: [0.4252, 0.4711] (0.4711-0.4252)/0.4482 = 0.1024
% K: [316061.3452, 372238.2150] (372238.2150-316061.3452)/344144.3826 =
% 0.1632
% nu: [0.9562, 5.6140] (5.6140-0.9562)/3.2848 = 1.418

%%Plot BIC for the models 
N = numel(Tykcpr_10y_all);  
k_log = length(p_log);
k_gen = length(p_gen);
SSE_log = resnorm_log; %here the resnorm is already squared
SSE_gen = resnorm_gen;

% BIC formula: N*ln(SSE/N) + k*ln(N)
BIC_log = N * log(SSE_log / N) + k_log * log(N);
BIC_gen = N * log(SSE_gen / N) + k_gen * log(N);

% Display results
fprintf('Logistic:     BIC = %.3f\n', BIC_log);
fprintf('Gen Logistic: BIC = %.3f\n', BIC_gen);

if BIC_log < BIC_gen
    fprintf('Logistic preferred (lower BIC)\n');
else
    fprintf('Gen Logistic preferred (lower BIC)\n');
end

% Logistic:     BIC = 661.814
% Gen Logistic: BIC = 654.531
% Gen Logistic preferred (lower BIC)

%%Plot the curves
% Dense time grid for plotting continuous curves
t_fine = linspace(0, 14, 200)';

% Predict for each replicate separately at t_fine, using fitted params
% Example: use the same ODE solvers but with t_fine instead of 1–14
% (for simplicity, here we interpolate from 1–14)

y_log_fit_full = logistic_ode_model(p_log, Tykcpr_10t_days, Tykcpr_10y_raw);
y_gen_fit_full = genlog_ode_model(p_gen, Tykcpr_10t_days, Tykcpr_10y_raw);

% Reshape back to 14×3 for plotting
y_log_fit = reshape(y_log_fit_full, [nDays, nReps]);
y_gen_fit = reshape(y_gen_fit_full, [nDays, nReps]);

figure;
hold on;
% Plot data
plot(Tykcpr_10t_days, Tykcpr_10y_raw, 'ko', 'MarkerFaceColor',[0.8 0.8 0.8]);

% Plot fitted curves per replicate (logistic)
for j = 1:nReps
    plot(Tykcpr_10t_days, y_log_fit(:,j), 'b-', 'LineWidth', 1.5);
end

% Plot generalized logistic as dashed red (optional)
for j = 1:nReps
    plot(Tykcpr_10t_days, y_gen_fit(:,j), 'r--', 'LineWidth', 1.5);
end

xlabel('Time (days)');
ylabel('Cell count (or density)');
legend('Data','Logistic fit','Gen. logistic fit','Location','best');
grid on;
hold off;

saveas(gcf, fullfile(outdir,'Tykcpr_untreated_logandgenlog_10days_new.png'));

writematrix(y_log_fit, fullfile(outdir,'Tykcpr_untreated_log_10days_new.csv'));
writematrix(y_gen_fit, fullfile(outdir,'Tykcpr_untreated_genlog_10days_new.csv'));

%calculate Rsquared values
clear all;
%%%%Sensitive cells first 
Tykcpr_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/data/Tykcpr_Untreated_40k_cellsinwell_new.csv');
Tykcpr_t_days = Tykcpr_data.Day;
Tykcpr_y_raw = [Tykcpr_data.Rep1, Tykcpr_data.Rep2, Tykcpr_data.Rep3];
Tykcpr_10y_raw = Tykcpr_y_raw(1:11,:);
Tykcpr_10y_all = Tykcpr_10y_raw(:); %should be 45 x 1 


Tykcpr_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tykcpr_untreated/Tykcpr_untreated_log_10days_new.csv');
Tykcpr_predicted = table2array(Tykcpr_predicted);
Tykcpr_predicted = Tykcpr_predicted(:);

SS_res = sum((Tykcpr_10y_all - Tykcpr_predicted).^2);
SS_tot = sum((Tykcpr_10y_all - mean(Tykcpr_10y_all)).^2);
R2 = 1 - SS_res / SS_tot;

fprintf('R-squared for logistic fit: %.4f\n', R2);
%R-squared for logistic fit: 0.9729

Tykcpr_predicted_gen = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/Tyk/clean_parameterization/results/Tykcpr_untreated/Tykcpr_untreated_genlog_10days_new.csv');
Tykcpr_predicted_gen = table2array(Tykcpr_predicted_gen);
Tykcpr_predicted_gen = Tykcpr_predicted_gen(:);

SS_res = sum((Tykcpr_10y_all - Tykcpr_predicted_gen).^2);
SS_tot = sum((Tykcpr_10y_all - mean(Tykcpr_10y_all)).^2);
R2 = 1 - SS_res / SS_tot;

fprintf('R-squared for generalized fit: %.4f\n', R2);
%R-squared for generalized fit: 0.9805