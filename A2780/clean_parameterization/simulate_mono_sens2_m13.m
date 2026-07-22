function S_model = simulate_mono_sens2_m13(t_data, S0, rS, K, alphaS, nS, IC50S, ce)
    ode_fun = @(t,S) mono_ode_sens2_m13(t, S, rS, K, alphaS, nS, IC50S, ce);
    opts = odeset('RelTol',1e-8, 'AbsTol',1e-10);
    [t_sol, S_sol] = ode15s(ode_fun, [min(t_data), max(t_data)], S0, opts);
    S_model = interp1(t_sol, S_sol, t_data, 'pchip');
end