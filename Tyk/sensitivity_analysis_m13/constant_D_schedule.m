function D_t = constant_D_schedule(~, ce)
    % Constant drug concentration at ALL time points
    D_t = ce * ones(size(~));
end