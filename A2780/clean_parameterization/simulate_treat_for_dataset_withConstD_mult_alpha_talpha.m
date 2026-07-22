function S_model = simulate_treat_for_dataset_withConstD_mult_alpha_talpha(n, alpha, t_data, S0, ce, rS, KS, IC50)
    tspan = [min(t_data), max(t_data)];
    ode_fun = @(t, S) treated_ode_withConstD_mult_alpha_talpha(t, S, n, alpha, ce, rS, KS, IC50);
    opts = odeset('RelTol',1e-8, 'AbsTol',1e-10);
    [t_sol, S_sol] = ode15s(ode_fun, tspan, S0, opts);
    S_model = interp1(t_sol, S_sol, t_data, 'pchip');
end