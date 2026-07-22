function R_model = simulate_mono_res_sens2_m13(t_data, R0, rR, K, alphaR, nR, IC50R, ce)
    ode_fun = @(t,R) mono_ode_res_sens2_m13(t, R, rR, K, alphaR, nR, IC50R, ce);
    opts = odeset('RelTol',1e-8, 'AbsTol',1e-10);
    [t_sol, R_sol] = ode15s(ode_fun, [min(t_data), max(t_data)], R0, opts);
    R_model = interp1(t_sol, R_sol, t_data, 'pchip');
end