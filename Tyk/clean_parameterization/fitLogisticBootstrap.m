function p_boot = fitLogisticBootstrap(t_boot, y_boot, p0, lb, ub)
    % Reshape to vectors
    t_boot_vec = reshape(t_boot, [], 1);
    y_boot_vec = reshape(y_boot, [], 1);
    
    % Create objective function for bootstrap sample
    obj_fun = @(p) logistic_ode_model(p, t_boot_vec, y_boot_vec) - y_boot_vec;
    
    % Fit to bootstrap sample using the provided initial guess
    p_boot = lsqnonlin(obj_fun, p0, lb, ub);
end