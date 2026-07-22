%%Trying to estimate parameters with MATLAB for untreated cells
%Compare logistic and generalized logistic models 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Repeat everything for 10 days and lsqnonlin, and see how the parameters
%%%change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%need to change it a little bit to be able to use the lsqnonlin 
clear all;
outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_untreated';
%%%%Sensitive cells first 
A2780_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_Untreated_30k_cellsinwell.csv');
A2780_t_days = A2780_data.Day;
A2780_10t_days = A2780_t_days(1:11);
A2780_y_raw = [A2780_data.Rep1, A2780_data.Rep2, A2780_data.Rep3];
A2780_10y_raw = A2780_y_raw(1:11,:);
[nDays, nReps] = size(A2780_10y_raw);

%Vectorize time and data so all replicates are used, no need to do this to
%use lsqnonlin
A2780_10t_all = repmat(A2780_10t_days, nReps, 1); %should be 45 x 1 
A2780_10y_all = A2780_10y_raw(:); %should be 45 x 1 
A2780_10N = numel(A2780_10y_all); %total observations so 45


%Define objective function for both models, need to define a residual
%function to be able to use lsqnonlin
logistic_res = @(p, t_dummy) logistic_ode_model(p, A2780_10t_days, A2780_10y_raw)-A2780_10y_all;

genlog_res   = @(p, t_dummy) genlog_ode_model(p, A2780_10t_days, A2780_10y_raw)-A2780_10y_all;

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
fprintf('Fitted logistic parameters for A2780 (10 days)\n')
fprintf('Fitted growth rate r = %.4f\n', p_log(2));
fprintf('Fitted carrying capacity k = %.4f\n', p_log(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_log);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted logistic parameters for A2780 (10 days)
% Fitted growth rate r = 0.8547
% Fitted carrying capacity k = 1114404.9924
% Sum of squared residuals = 295866465375.16

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_log_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(logistic_res, p0_log, lb_log, ub_log);

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
% 95r: [0.7903, 0.9190] (0.9190 - 0.7903)/0.8547 = 0.1506
% K: [1046016.8381, 1182793.1467] (1182793.1467 -
% 1046016.8381)/1114404.9924 = 0.1227

%Define initial guesses for generalized logistic
K0  = 1.4e6;         %max(y_all);
r0  = 0.5;
nu0 = 1;            % nu=1 reduces to logistic

p0_gen = [K0, r0, nu0];
lb_gen = [0,  0,   0.1];   % nu > 0
ub_gen = [Inf, Inf, 10];

problem_gen = createOptimProblem( ...
    'lsqnonlin', ...
    'x0',        p0_gen, ...
    'objective', genlog_res, ...
    'lb',        lb_gen, ...
    'ub',        ub_gen);

[p_gen, resnorm_gen] = run(ms, problem_gen, nStarts);
fprintf('Fitted generalized logistic parameters for A2780 (10 days)\n')
fprintf('Fitted growth rate r = %.4f\n', p_gen(2));
fprintf('Fitted carrying capacity k = %.4f\n', p_gen(1));
fprintf('Fitted n parameter n = %.4f\n', p_gen(3));
fprintf('Sum of squared residuals = %.2f\n', resnorm_gen);
% 197 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted generalized logistic parameters for A2780 (10 days)
% Fitted growth rate r = 0.7813
% Fitted carrying capacity k = 1095502.8261
% Fitted n parameter n = 1.3332
% Sum of squared residuals = 291704500733.76

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
% lsqnonlin stopped because the size of the current step is less than
% the value of the step size tolerance.
% <stopping criteria details>
% 95r: [0.5942, 0.9683] (0.9683-0.5942)/0.7813 = 0.4788
% K: [1015315.3794, 1175690.2728] (1175690.2728-1015315.3794)/1095502.8261
% = 0.1464
% nu: [0.2116, 2.4549] (2.4549-0.2116)/ 1.3332 = 1.6826


%%Plot BIC for the models 
N = numel(A2780_10y_all);  
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

% Logistic:     BIC = 763.243
% Gen Logistic: BIC = 766.272
% Logistic preferred (lower BIC)

%%Plot the curves
% Dense time grid for plotting continuous curves
t_fine = linspace(0, 14, 200)';

% Predict for each replicate separately at t_fine, using fitted params
% Example: use the same ODE solvers but with t_fine instead of 1–14
% (for simplicity, here we interpolate from 1–14)

y_log_fit_full = logistic_ode_model(p_log, A2780_10t_days, A2780_10y_raw);
y_gen_fit_full = genlog_ode_model(p_gen, A2780_10t_days, A2780_10y_raw);

% Reshape back to 14×3 for plotting
y_log_fit = reshape(y_log_fit_full, [nDays, nReps]);
y_gen_fit = reshape(y_gen_fit_full, [nDays, nReps]);

figure;
hold on;
% Plot data
plot(A2780_10t_days, A2780_10y_raw, 'ko', 'MarkerFaceColor',[0.8 0.8 0.8]);

% Plot fitted curves per replicate (logistic)
for j = 1:nReps
    plot(A2780_10t_days, y_log_fit(:,j), 'b-', 'LineWidth', 1.5);
end

% Plot generalized logistic as dashed red (optional)
for j = 1:nReps
    plot(A2780_10t_days, y_gen_fit(:,j), 'r--', 'LineWidth', 1.5);
end

xlabel('Time (days)');
ylabel('Cell count (or density)');
legend('Data','Logistic fit','Gen. logistic fit','Location','best');
grid on;
hold off;

saveas(gcf, fullfile(outdir,'A2780_untreated_logandgenlog_10days.png'));

writematrix(y_log_fit, fullfile(outdir,'A2780_untreated_log_10days.csv'));
writematrix(y_gen_fit, fullfile(outdir,'A2780_untreated_genlog_10days.csv'));

%calculate Rsquared values 
clear all;
%%%%Sensitive cells first 
A2780_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780_Untreated_30k_cellsinwell.csv');
A2780_t_days = A2780_data.Day;
A2780_y_raw = [A2780_data.Rep1, A2780_data.Rep2, A2780_data.Rep3];
A2780_10y_raw = A2780_y_raw(1:11,:);
A2780_10y_all = A2780_10y_raw(:); %should be 33 x 1 


A2780_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_untreated/A2780_untreated_log_10days.csv');
A2780_predicted = table2array(A2780_predicted);
A2780_predicted = A2780_predicted(:); %33x1

SS_res = sum((A2780_10y_all - A2780_predicted).^2);
SS_tot = sum((A2780_10y_all - mean(A2780_10y_all)).^2);
R2 = 1 - SS_res / SS_tot;

fprintf('R-squared for logistic fit: %.4f\n', R2);
%R-squared for logistic fit: 0.9522

A2780_predicted_gen = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780_untreated/A2780_untreated_genlog_10days.csv');
A2780_predicted_gen = table2array(A2780_predicted_gen);
A2780_predicted_gen = A2780_predicted_gen(:);

SS_res = sum((A2780_10y_all - A2780_predicted_gen).^2);
SS_tot = sum((A2780_10y_all - mean(A2780_10y_all)).^2);
R2 = 1 - SS_res / SS_tot;

fprintf('R-squared for generalized fit: %.4f\n', R2);
%R-squared for generalized fit: 0.9529


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Repeat the 10 day process for A2780cis with lsqnonlin
clear all;
outdir = '/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780cis_untreated';

%%%%Resistant cells cells  
A2780cis_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780cis_Untreated_30k_cellsinwell.csv');
A2780cis_t_days = A2780cis_data.Day;
A2780cis_10t_days = A2780cis_t_days(1:11);
A2780cis_y_raw = [A2780cis_data.Rep1, A2780cis_data.Rep2, A2780cis_data.Rep3];
A2780cis_10y_raw = A2780cis_y_raw(1:11,:);
[nDays, nReps] = size(A2780cis_10y_raw);

%Vectorize time and data so all replicates are used, no need to do this to
%use lsqnonlin
A2780cis_10t_all = repmat(A2780cis_10t_days, nReps, 1); %should be 45 x 1 
A2780cis_10y_all = A2780cis_10y_raw(:); %should be 45 x 1 
A2780cis_10N = numel(A2780cis_10y_all); %total observations so 45

%Define objective function for both models, need to define a residual
%function to be able to use lsqnonlin
logistic_res = @(p, t_dummy) logistic_ode_model(p, A2780cis_10t_days, A2780cis_10y_raw)-A2780cis_10y_all;

genlog_res   = @(p, t_dummy) genlog_ode_model(p, A2780cis_10t_days, A2780cis_10y_raw)-A2780cis_10y_all;

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
fprintf('Fitted logistic parameters for A2780cis (10 days)\n')
fprintf('Fitted growth rate r = %.4f\n', p_log(2));
fprintf('Fitted carrying capacity k = %.4f\n', p_log(1));
fprintf('Sum of squared residuals = %.2f\n', resnorm_log);
% All 200 local solver runs converged with a positive local solver exitflag.
% Fitted logistic parameters for A2780cis (10 days)
% Fitted growth rate r = 0.5553
% Fitted carrying capacity k = 1356851.0591
% Sum of squared residuals = 166036884856.39

% Run lsqnonlin separately to get Jacobian using MultiStart solution as initial guess
[p_log_final, resnorm_final, residual_final, exitflag, output_final, lambda_final, J] = ...
lsqnonlin(logistic_res, p0_log, lb_log, ub_log);

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
% 95r: [0.5233, 0.5873] (0.5873 - 0.5233)/0.5553 = 0.1153
% K: [1206767.4755, 1506827.3486] (1506827.3486 -
% 1206767.4755)/1356851.0591 = 0.2211

%Define initial guesses for generalized logistic
K0  = 1.4e6;         %max(y_all);
r0  = 0.5;
nu0 = 1;            % nu=1 reduces to logistic

p0_gen = [K0, r0, nu0];
lb_gen = [0,  0,   0.1];   % nu > 0
ub_gen = [Inf, Inf, 10];

problem_gen = createOptimProblem( ...
    'lsqnonlin', ...
    'x0',        p0_gen, ...
    'objective', genlog_res, ...
    'lb',        lb_gen, ...
    'ub',        ub_gen);

[p_gen, resnorm_gen] = run(ms, problem_gen, nStarts);
fprintf('Fitted generalized logistic parameters for A2780cis (10 days)\n')
fprintf('Fitted growth rate r = %.4f\n', p_gen(2));
fprintf('Fitted carrying capacity k = %.4f\n', p_gen(1));
fprintf('Fitted n term n = %.4f\n', p_gen(3));
fprintf('Sum of squared residuals = %.2f\n', resnorm_gen);
% 190 out of 200 local solver runs converged with a positive local solver exitflag.
% Fitted generalized logistic parameters for A2780cis (10 days)
% Fitted growth rate r = 0.4621
% Fitted carrying capacity k = 1093175.1044
% Fitted n term n = 4.5931
% Sum of squared residuals = 66862349478.50

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
% 95r: [0.4478, 0.4764] (0.4764-0.4478)/0.4621 = 0.0619
% K: [1042705.8743, 1143644.3188] (1143644.3188-1042705.8743)/1093175.1044
% = 0.0923
% nu: [1.6989, 7.4873] (7.4873-1.6989)/4.5931 = 1.26

%%Plot BIC for the models 
N = numel(A2780cis_10y_all);  
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

% Logistic:     BIC = 744.179
% Gen Logistic: BIC = 717.660
% Gen Logistic preferred (lower BIC)

%%Plot the curves
% Dense time grid for plotting continuous curves
t_fine = linspace(0, 14, 200)';

% Predict for each replicate separately at t_fine, using fitted params
% Example: use the same ODE solvers but with t_fine instead of 1–14
% (for simplicity, here we interpolate from 1–14)

y_log_fit_full = logistic_ode_model(p_log, A2780cis_10t_days, A2780cis_10y_raw);
y_gen_fit_full = genlog_ode_model(p_gen, A2780cis_10t_days, A2780cis_10y_raw);

% Reshape back to 14×3 for plotting
y_log_fit = reshape(y_log_fit_full, [nDays, nReps]);
y_gen_fit = reshape(y_gen_fit_full, [nDays, nReps]);

figure;
hold on;
% Plot data
plot(A2780cis_10t_days, A2780cis_10y_raw, 'ko', 'MarkerFaceColor',[0.8 0.8 0.8]);

% Plot fitted curves per replicate (logistic)
for j = 1:nReps
    plot(A2780cis_10t_days, y_log_fit(:,j), 'b-', 'LineWidth', 1.5);
end

% Plot generalized logistic as dashed red (optional)
for j = 1:nReps
    plot(A2780cis_10t_days, y_gen_fit(:,j), 'r--', 'LineWidth', 1.5);
end

xlabel('Time (days)');
ylabel('Cell count (or density)');
legend('Data','Logistic fit','Gen. logistic fit','Location','best');
grid on;
hold off;

saveas(gcf, fullfile(outdir,'A2780cis_untreated_logandgenlog_10days.png'));

writematrix(y_log_fit, fullfile(outdir,'A2780cis_untreated_log_10days.csv'));
writematrix(y_gen_fit, fullfile(outdir,'A2780cis_untreated_genlog_10days.csv'));

%calculate Rsquared values
clear all;
%%%%Sensitive cells first 
A2780cis_data = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/data/A2780cis_Untreated_30k_cellsinwell.csv');
A2780cis_t_days = A2780cis_data.Day;
A2780cis_y_raw = [A2780cis_data.Rep1, A2780cis_data.Rep2, A2780cis_data.Rep3];
A2780cis_10y_raw = A2780cis_y_raw(1:11,:);
A2780cis_10y_all = A2780cis_10y_raw(:); %should be 45 x 1 


A2780cis_predicted = readtable('/blue/ferrallm/01_analysis/adelpinoherrera_OvCa-Aim2-MathModelingMATLAB/clean_parameterization/results/A2780cis_untreated/A2780cis_untreated_log_10days.csv');
A2780cis_predicted = table2array(A2780cis_predicted);
A2780cis_predicted = A2780cis_predicted(:);

SS_res = sum((A2780cis_10y_all - A2780cis_predicted).^2);
SS_tot = sum((A2780cis_10y_all - mean(A2780cis_10y_all)).^2);
R2 = 1 - SS_res / SS_tot;

fprintf('R-squared for logistic fit: %.4f\n', R2);
%R-squared for logistic fit: 0.9719