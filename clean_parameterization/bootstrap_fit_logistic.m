function p_boot = bootstrap_fit_logistic(idx, data)
    % Resample data using indices
    t_boot = data{1}(idx);
    y_boot = data{2}(idx);
    
    % Create objective function for bootstrap sample
    t_boot_vec = reshape(t_boot, [], 1);
    y_boot_vec = reshape(y_boot, [], 1);
    
    obj_fun = @(p) logistic_ode_model(p, t_boot_vec, y_boot_vec) - y_boot_vec;
    
    % Known parameters
    lb_log = [0, 0];
    ub_log = [Inf, Inf];
    p0 = [1.4e6, 0.5];
    
    % Fit to bootstrap sample
    p_boot = lsqnonlin(obj_fun, p0, lb_log, ub_log);
end